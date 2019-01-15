#include <fstream>
#include <stdio.h>
#include <commonmath.h>
#include <numericalintegration.h>
#include <convexoptimization.h>
using namespace ContinuumRobotLibrary;

//Independent Parameters
const double E = 200e9;
const double G = 80e9;
const double rad = 0.001;
const double rho = 8000;
const Vector3d g = 9.81*Vector3d::UnitX();
const double L = 0.5;
const int num_tendons = 4;

//Helical routing
const double offset = 0.01506;
void getRouting(double arclength, Vector3d* r, Vector3d* r_s, Vector3d* r_ss){
    const double f = 2*pi/L; //tendons make a full revolution

    double c = offset*cos(f*arclength);
    double s = offset*sin(f*arclength);
    double fc = f*c;
    double fs = f*s;
    double ffc = f*fc;
    double ffs = f*fs;

    r[0] = Vector3d(c, s, 0);
    r[1] = Vector3d(-s, c, 0);
    r[2] = Vector3d(-c, -s, 0);
    r[3] = Vector3d(s, -c, 0);

    r_s[0] = Vector3d(-fs, fc, 0);
    r_s[1] = Vector3d(-fc, -fs, 0);
    r_s[2] = Vector3d(fs, -fc, 0);
    r_s[3] = Vector3d(fc, fs, 0);

    r_ss[0] = Vector3d(-ffc, -ffs, 0);
    r_ss[1] = Vector3d(ffs, -ffc, 0);
    r_ss[2] = Vector3d(ffc, ffs, 0);
    r_ss[3] = Vector3d(-ffs, ffc, 0);
}

const VectorXd tau = Vector4d(15,0,0,0);

//Dependent parameter calculations
const double area = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
const DiagonalMatrix<double, 3> Kse (G*area,G*area,E*area);
const DiagonalMatrix<double, 3> Kbt (E*I,E*I,G*J);
const Matrix3d Kse_dense = Kse.toDenseMatrix();
const Matrix3d Kbt_dense = Kbt.toDenseMatrix();

void cosseratTendonRobotOde(VectorXd& y_s_out, double s, VectorXd& y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d v = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);

    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();
    Matrix3d A_plus_Kse = Kse_dense;
    Matrix3d G = Matrix3d::Zero();
    Matrix3d H_plus_Kbt = Kbt_dense;

    Vector3d r[num_tendons], r_s[num_tendons], r_ss[num_tendons];
    getRouting(s, r, r_s, r_ss);

    for(int i = 0; i < num_tendons; i++){
        Vector3d ri = r[i];
        Vector3d ri_s = r_s[i];
        Vector3d ri_ss = r_ss[i];

        Vector3d pb_si = u.cross(ri) + ri_s + v;
        double pb_s_norm = pb_si.norm();
        Matrix3d A_i = -hat_squared(pb_si)*(tau(i)/pow(pb_s_norm,3));
        Matrix3d G_i = -hat_postmultiply(A_i,ri);
        Vector3d a_i = A_i*(u.cross(pb_si + ri_s) + ri_ss);

        a += a_i;
        b += ri.cross(a_i);
        A_plus_Kse += A_i;
        G += G_i;
        H_plus_Kbt += hat_premultiply(ri,G_i);
    }

    Matrix6d K;
    K << A_plus_Kse, G, G.transpose(), H_plus_Kbt;

    Vector3d nb = Kse*(v - Vector3d::UnitZ());
    Vector3d mb = Kbt*u;

    Vector6d rhs;
    rhs << -u.cross(nb) - transposeMultiply(R,rho*area*g) - a,
           -u.cross(mb) - v.cross(nb) - b;

    //Pack state vector derivative
    Map<Vector3d> p_s(&y_s_out[0]);
    Map<Matrix3d> R_s(&y_s_out[3]);
    Map<Vector6d> vs_and_us(&y_s_out[12]);

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    vs_and_us = K.selfadjointView<Eigen::Upper>().llt().solve(rhs);
}

//Boundary conditions
const Vector3d p0 = Vector3d::Zero();
const Matrix3d R0 = Matrix3d::Identity();

static MatrixXd Y; //Declare Y global to save results
VectorXd shootingFunction(VectorXd guess){
    Vector3d n0 = guess.segment<3>(0);
    Vector3d v0 = Kse.inverse()*n0 + Vector3d::UnitZ();
    Vector3d u0 = guess.segment<3>(3);

    VectorXd y0(18);
    y0 << p0, Map<VectorXd>(Matrix3d(R0).data(), 9), v0, u0;

    //Numerically integrate the Cosserat rod equations
    Y = ode4<cosseratTendonRobotOde>(y0, 0, L);

    //Find the internal forces in the backbone prior to the final plate
    Vector3d vL = Y.block<3,1>(12,Y.cols()-1);
    Vector3d uL = Y.block<3,1>(15,Y.cols()-1);

    Vector3d nb = Kse*(vL - Vector3d::UnitZ());
    Vector3d mb = Kbt*uL;

    //Find the equilibrium error at the tip, considering tendon forces
    Vector3d force_error = -nb;
    Vector3d moment_error = -mb;
    Vector3d r[num_tendons], r_s[num_tendons], r_ss[num_tendons];
    getRouting(L, r, r_s, r_ss);
    for(int i = 0; i < num_tendons; i++){
        Vector3d pb_si = uL.cross(r[i]) + r_s[i] + vL;
        Vector3d Fb_i = -tau(i)*pb_si.normalized();
        force_error += Fb_i;
        moment_error += r[i].cross(Fb_i);
    }

    Vector6d distal_error;
    distal_error << force_error, moment_error;

    return distal_error;
}

int main(int, char**){
    Vector6d init_guess = Vector6d::Zero(); //nb and u

    //Solve with shooting method
    VectorXd wrench_soln = solveLevenbergMarquardt<shootingFunction>(init_guess, 1e-12, 500, 1e-2, 0.5, 1e-7, 1e-9);

    //Save results for Blender visualization
    std::fstream file("centerline.dat", std::fstream::out);
    file << Y.block(0,0,3,Y.cols());
    file.close();

    MatrixXd tendonlines(3*num_tendons, Y.cols());
    for(int i = 0; i < Y.cols(); i++){
        Vector3d p = Y.block<3,1>(0,i);
        Matrix3d R = Map<Matrix3d>(&Y(3,i));
        Vector3d r[num_tendons], r_s[num_tendons], r_ss[num_tendons];
        getRouting(i*L/(Y.cols()-1), r, r_s, r_ss);
        for(int j = 0; j < num_tendons; j++)
            tendonlines.block<3,1>(3*j,i) = p + R*r[j];
    }
    file = std::fstream("tendonlines.dat", std::fstream::out);
    file << tendonlines;
    file.close();

    const int num_disks = 9;
    Matrix3Xd disks(3,4*num_disks);
    for(int i = 1; i <= num_disks; i++){
        int j = ((Y.cols()-1) * i) / num_disks;
        Vector3d p = Y.block<3,1>(0,j);
        Matrix3d R = Map<Matrix3d>(&Y(3,j));
        disks.block<3,3>(0,4*(i-1)) = R;
        disks.block<3,1>(0,4*(i-1)+3) = p;
    }
    file = std::fstream("disks.dat", std::fstream::out);
    file << disks;
    file.close();

    return 0;
}
