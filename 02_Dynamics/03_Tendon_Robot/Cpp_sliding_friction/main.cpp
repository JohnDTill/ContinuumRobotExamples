#include <fstream>
#include <iostream>
#include <stdio.h>
#include <commonmath.h>
#include <convexoptimization.h>
#include <numericalpde.h>
#include <timemanagement.h>
using namespace ContinuumRobotLibrary;

#ifdef QT_CORE_LIB
#include <simpleplotting.h>
#endif

//Independent Parameters
const double nu = 0.1;
const double E = 207e9;
const double rad = 0.00135/2;
const double total_mass = 0.034;
const Vector3d g = -9.81*Vector3d::UnitX();
const double L = 0.7144;
const double base_to_motor = 0.0518; //m
const double dt = 5e-2;
const double alpha = -0.4;
const double T = 25;
const int N = 200;
const int num_tendons = 1;
typedef Array<double, num_tendons, 1> ArrayNd;
const ArrayNd compliance = ArrayNd::Constant(1.6e-3);
const DiagonalMatrix<double, 3> Bse (0, 0, 0);
const DiagonalMatrix<double, 3> Bbt (5e-4, 5e-4, 5e-4);
const DiagonalMatrix<double, 3> C (1e-4, 1e-4, 1e-4);
const double tendon_offset = 0.01506;
const Vector3d r[num_tendons] = { tendon_offset*Vector3d::UnitX() };

const double zA = 0;          //   z(m)  Ramp Motor Input:
const double zB = -0.01619;   // zA ^ ____           ______
const double t1 = 3.652;      //    |     \         /
const double t2 = 3.967;      // zB |      \_______/
const double t3 = 17.63;      //    -------------------------> t(s)
const double t4 = 17.94;      //         t1 t2    t3 t4
double z(double t){
    if(t>t1 && t<=t2)       return zA + (zB-zA)*(t-t1)/(t2-t1); //Ramp lower
    else if(t>t2 && t<=t3)  return zB;                          //Low state
    else if(t>t3 && t<=t4)  return zB + (zA-zB)*(t-t3)/(t4-t3); //Ramp higher
    else                    return zA;                          //High state
}

//Dependent parameter calculations
const double c0 = TimeManagerBdfAlpha::getC0(dt,alpha);
const double G = E/(2*1.3);
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double rho = total_mass/(L*A);
const DiagonalMatrix<double, 3> Kse (G*A,G*A,E*A);
const DiagonalMatrix<double, 3> Kbt (E*I,E*I,G*2*I);
const DiagonalMatrix<double, 3> J (I, I, 2*I);
const Matrix3d Kse_dense = Kse.toDenseMatrix();
const Matrix3d Kbt_dense = Kbt.toDenseMatrix();
const Matrix3d Kse_c0Bse = Kse.toDenseMatrix() + c0*Bse.toDenseMatrix();
const Matrix3d Kbt_c0Bbt = Kbt.toDenseMatrix() + c0*Bbt.toDenseMatrix();
const Vector3d rhoAg = rho*A*g;
const double rhoA = rho*A;
const DiagonalMatrix<double, 3> rhoJ = rho*J;
static double t = 0;

double smoothSaturation(double val){
    const double eps = 1e-3;
    return val/sqrt(val*val + eps*eps);
}

static Matrix3d R;
static Vector3d v;
static Vector3d u;
static Vector3d q;
static Vector3d w;
static Vector3d q_t;
static Vector3d w_t;
static Vector3d nb;
static Vector3d mb;
static ArrayNd tau;
static ArrayNd tau_s;
static Vector3d v_sh;
static Vector3d u_sh;
static ArrayNd si_t;
static ArrayNd si_s;
VectorXd odeObjFunc(VectorXd guess){
    //Probably better conditioned to guess nb_s and mb_s
    //- In fact, that might avoid your conditioning problem
    Vector3d v_s = guess.head(3);
    Vector3d u_s = guess.tail(3);

    Vector3d nb_s = rhoA*(w.cross(q) + q_t) + C*q.cwiseProduct(q.cwiseAbs()) - transposeMultiply(R,rhoAg) - u.cross(nb);
    Vector3d mb_s = w.cross(rhoJ*w) + rhoJ*w_t - v.cross(nb) - u.cross(mb);

    for(int i = 0; i < num_tendons; i++){
        Vector3d pbsi = u.cross(r[i]) + v;
        Vector3d pbssi = u.cross(pbsi) + u_s.cross(r[i]) + v_s;
        double pbsi_norm = pbsi.norm();
        Vector3d pbsi_normalized = pbsi/pbsi_norm;
        Vector3d pbsi_normalized_s = (Matrix3d::Identity() - pbsi_normalized*pbsi_normalized.transpose()) / pbsi_norm * pbssi;

        Vector3d fi_b_star = -tau(i) * pbsi_normalized_s;
        double fi_b_norm = (fi_b_star - fi_b_star.dot(pbsi_normalized) * pbsi_normalized).norm();

        si_s(i) = pbsi_norm / (1 + compliance(i)*tau(i));
        tau_s(i) = -nu*smoothSaturation(si_t(i))*si_s(i)*fi_b_norm;

        Vector3d fi_b = -tau_s(i) * pbsi.normalized() + fi_b_star;
        Vector3d li_b = r[i].cross(fi_b);

        nb_s -= (-fi_b);
        mb_s -= (-li_b);
    }

    Vector3d v_s_calculated = Kse_c0Bse.inverse() * (nb_s - Bse*v_sh);
    Vector3d u_s_calculated = Kbt_c0Bbt.inverse() * (mb_s - Bbt*u_sh);

    Vector6d err;
    err << v_s_calculated - v_s, u_s_calculated - u_s;

    return err;
}

//PDE semi-discretized in time describing tendon backbone dynamics
void tendonBackbonePDE(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y, Map<VectorXd> z_h){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d v = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);
    Vector3d q = Map<Vector3d>(&y[18]);
    Vector3d w = Map<Vector3d>(&y[21]);
    ArrayNd si = Map<ArrayNd>(&y[24]);
    ArrayNd tau = Map<ArrayNd>(&y[24+num_tendons]);

    //Map state vector derivative
    Map<Vector3d> p_s = Map<Vector3d>(&y_s_out[0]);
    Map<Matrix3d> R_s = Map<Matrix3d>(&y_s_out[3]);
    Map<Vector6d> vs_and_us(&y_s_out[12]);
    Map<Vector3d> q_s = Map<Vector3d>(&y_s_out[18]);
    Map<Vector3d> w_s = Map<Vector3d>(&y_s_out[21]);
    Map<ArrayNd> si_s(&y_s_out[24]);
    Map<ArrayNd> tau_s(&y_s_out[24+num_tendons]);

    Map<Vector3d> v_h(&z_h[0]);
    Map<Vector3d> u_h(&z_h[3]);
    Map<Vector3d> q_h(&z_h[6]);
    Map<Vector3d> w_h(&z_h[9]);
    Map<Vector3d> v_sh(&z_h[12]);
    Map<Vector3d> u_sh(&z_h[15]);
    Map<ArrayNd> si_h(&z_h[18]);

    Vector3d v_t = c0*v + v_h;
    Vector3d u_t = c0*u + u_h;
    Vector3d q_t = c0*q + q_h;
    Vector3d w_t = c0*w + w_h;
    ArrayNd si_t = c0*si + si_h;

    Vector3d nb = Kse*(v - Vector3d::UnitZ()) + Bse*v_t;
    Vector3d mb = Kbt*u + Bbt*u_t;

    ::R = R;
    ::v = v;
    ::u = u;
    ::q = q;
    ::w = w;
    ::q_t = q_t;
    ::w_t = w_t;
    ::tau = tau;
    ::nb = nb;
    ::mb = mb;
    ::v_sh = v_sh;
    ::u_sh = u_sh;
    ::tau_s = ArrayNd::Zero();
    ::si_t = si_t;

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    //std::cout << "solving..." << std::endl;
    vs_and_us = solveLevenbergMarquardt<odeObjFunc>(Vector6d::Zero(),1e-20,500,1e-2,0.5,1e-5,1e-6);
    //std::cout << "solved" << std::endl;
    q_s = v_t - u.cross(q) + w.cross(v);
    w_s = u_t - u.cross(w);
    tau_s = ::tau_s;
    si_s = ::si_s;

    //std::cout << tau_s.transpose() << std::endl;

    //Output argument for variables with time derivatives
    z_out << v, u, q, w, vs_and_us, si;
}

//ODE describing tendon robot statics
void cosseratTendonRobotOde(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d v = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);
    ArrayNd si = Map<ArrayNd>(&y[24]);
    ArrayNd tau = Map<ArrayNd>(&y[24+num_tendons]);

    //Map state vector derivative
    Map<Vector3d> p_s(&y_s_out[0]);
    Map<Matrix3d> R_s(&y_s_out[3]);
    Map<Vector6d> vs_and_us(&y_s_out[12]);
    Map<Vector3d> q_s = Map<Vector3d>(&y_s_out[18]);
    Map<Vector3d> w_s = Map<Vector3d>(&y_s_out[21]);
    Map<ArrayNd> si_s(&y_s_out[24]);
    Map<ArrayNd> tau_s(&y_s_out[24+num_tendons]);

    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();
    Matrix3d A_plus_Kse = Kse_dense;
    Matrix3d G = Matrix3d::Zero();
    Matrix3d H_plus_Kbt = Kbt_dense;

    for(int i = 0; i < num_tendons; i++){
        Vector3d pb_si = u.cross(r[i]) + v;
        double pb_s_norm = pb_si.norm();
        si_s(i) = pb_s_norm / (1 + compliance(i)*tau(i));
        Matrix3d A_i = -hat_squared(pb_si)*(tau(i)/pow(pb_s_norm,3));
        Matrix3d G_i = -hat_postmultiply(A_i,r[i]);
        Vector3d a_i = A_i*(u.cross(pb_si));

        a += a_i;
        b += r[i].cross(a_i);
        A_plus_Kse += A_i;
        G += G_i;
        H_plus_Kbt += hat_premultiply(r[i],G_i);
    }

    Matrix6d K;
    K << A_plus_Kse, G, G.transpose(), H_plus_Kbt;

    Vector3d nb = Kse*(v - Vector3d::UnitZ());
    Vector3d mb = Kbt*u;

    Vector6d rhs;
    rhs << -u.cross(nb) - transposeMultiply(R,rhoAg) - a,
           -u.cross(mb) - v.cross(nb) - b;

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    vs_and_us = K.selfadjointView<Eigen::Upper>().llt().solve(rhs);
    q_s = Vector3d::Zero();
    w_s = Vector3d::Zero();
    tau_s = ArrayNd::Zero();

    //Output argument for variables with time derivatives
    z_out << v, u, Vector6d::Zero(), vs_and_us, si;
}

static MatrixXd Y(24+2*num_tendons,N), Z(18+num_tendons,N), Z_h(18+num_tendons,N);
template<bool is_dynamic>
VectorXd obj(VectorXd guess){
    Vector3d v0 = Kse.inverse()*guess.head(3) + Vector3d::UnitZ(); //not exactly guessing n0 due to viscoelastic constitutive equation
    Vector3d u0 = guess.segment<3>(3);
    //ArrayNd tau0 = guess.tail<num_tendons>().cwiseMax(0);
    //ArrayNd slack = -(guess.tail<num_tendons>().cwiseMin(0));
    ArrayNd tau0 = guess.tail<num_tendons>();

    VectorXd y0(24+2*num_tendons);
    y0 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, v0, u0, Vector6d::Zero(), base_to_motor / (1+compliance(0)*tau0(0)), tau0;
    Y.col(0) = y0;

    //Numerically integrate the Cosserat rod equations
    if(is_dynamic) TimeBdfAlpha_SpaceEuler<tendonBackbonePDE,24+2*num_tendons,18,N>(Y,Z,y0,L,Z_h);
    else euler<cosseratTendonRobotOde,24+2*num_tendons,18,N>(Y,Z,y0,L);

    //Find the internal forces in the backbone prior to the final plate
    Vector3d vL = Y.block<3,1>(12,N-1);
    Vector3d uL = Y.block<3,1>(15,N-1);
    ArrayNd tauL = Y.bottomRightCorner<num_tendons,1>();

    Vector3d vL_t = is_dynamic ? (c0*vL + Z_h.block<3,1>(0,N-1)).eval() : Vector3d::Zero();
    Vector3d uL_t = is_dynamic ? (c0*uL + Z_h.block<3,1>(3,N-1)).eval() : Vector3d::Zero();

    Vector3d nbL = Kse*(vL - Vector3d::UnitZ()) + Bse*vL_t;
    Vector3d mbL = Kbt*uL + Bbt*uL_t;

    //Find the equilibrium error at the tip, considering tendon forces
    Vector3d force_error = -nbL;
    Vector3d moment_error = -mbL;
    for(int i = 0; i < num_tendons; i++){
        Vector3d pb_si = uL.cross(r[i]) + vL;
        Vector3d Fb_i = -tauL(i)*pb_si.normalized();
        force_error += Fb_i;
        moment_error += r[i].cross(Fb_i);
    }

    //Find the length violation error
    ArrayNd si_L = Y.block<num_tendons,1>(24,N-1);
    ArrayNd l_star = ArrayNd::Constant(L+base_to_motor+z(t));
    ArrayNd length_error = si_L - l_star;

    VectorXd distal_error(6 + num_tendons);
    distal_error << force_error, moment_error, length_error;

    return distal_error;
}

int main(int, char**){
    //Solve static
    VectorXd guess = VectorXd::Zero(6 + num_tendons);
    guess = solveLevenbergMarquardt<obj<false> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);

    //Solve dynamic
    TimeManagerBdfAlpha time_scheme(Z, Z_h, dt, alpha);
    int M = static_cast<int>(T/dt);
    MatrixXd tendon(3*M,N+1);
    Matrix3Xd tip(3,M);
    MatrixXd centerline(3*M,N);
    MatrixXd disks(3*M,4*7);

    for(int i = 0; i < M; i++){
        guess = solveLevenbergMarquardt<obj<true> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);
        std::cout << t << std::endl;

        t += dt;
        time_scheme.advanceTime();

        //Store results
        tip.col(i) = Y.block<3,1>(0,N-1);
        centerline.block<3,N>(3*i,0) = Y.block<3,N>(0,0);
        tendon.block<3,1>(3*i,0) = r[0] + z(t)*Vector3d::UnitZ();
        for(int j = 0; j < N; j++){
            Matrix3d R = Map<Matrix3d>(&Y(3,j));
            tendon.block<3,1>(3*i,j+1) = centerline.block<3,1>(3*i,j) + R*r[0];
        }
        for(int j = 1; j <= 7; j++){
            int k = (N-1)*j/7;
            disks.block<3,1>(3*i,4*(j-1)+3) = centerline.block<3,1>(3*i,k);
            disks.block<3,3>(3*i,4*(j-1)) = Map<Matrix3d>(&Y(3,k));
        }
    }

    //Save results for Blender visualization
    std::fstream file("tendon.dat", std::fstream::out);
    file << tendon;
    file.close();

    file = std::fstream("disks.dat", std::fstream::out);
    file << disks;
    file.close();

    file = std::fstream("centerline.dat", std::fstream::out);
    file << centerline;
    file.close();

    //Show the end-effector trajectory
    #ifdef QT_CORE_LIB
    plot(VectorXd::LinSpaced(M,0,t), tip.row(0), "Tip Trajectory", "t (s)", "x (m)");
    #endif

    return 0;
}
