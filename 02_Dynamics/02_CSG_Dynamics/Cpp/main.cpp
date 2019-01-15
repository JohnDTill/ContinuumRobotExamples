#include <fstream>
#include <iostream>
#include <stdio.h>
#include <ctime>
#include <commonmath.h>
#include <convexoptimization.h>
#include <numericalpde.h>
#include <timemanagement.h>
using namespace ContinuumRobotLibrary;

#ifdef QT_CORE_LIB
#include <simpleplotting.h>
#endif

//Independent Parameters
const double E = 207e9;
const double rad = 0.001;
const double rho = 8000;
const Vector3d g = -9.81*Vector3d::UnitZ();
const double scrib_R = 0.087; //radius of the hole pattern
const double alpha1 = 100*pi/180; //major angle of the hole pattern
const double plate_mass = 0.0921;
const DiagonalMatrix<double, 3> plate_mass_moment(1.91e-4, 1.91e-4, 3.82e-4);
const int N = 200; //Number of spatial points
const double dt = 1.0/120.0; //Time step
const double alpha = -0.2; //Parameter for implicit time discretization
const Vector3d p_start(0, 0, 0.5);
const Vector3d p_final(-0.2, 0, 0.5);
const int num_actuated = 300;
const int num_swaying = 300;

//Dependent parameter calculations
const double c0 = TimeManagerBdfAlpha::getC0(dt,alpha);
const double G = E/(2*(1+0.3));
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
static Vector3d pE = p_start;
static Matrix3d RE = Matrix3d::Identity();
const double alpha2 = 120*pi/180 - alpha1; //minor angle of the hole pattern
const DiagonalMatrix<double, 3> Kse_inv = DiagonalMatrix<double, 3>(G*A,G*A,E*A).inverse();
const DiagonalMatrix<double, 3> Kbt_inv = DiagonalMatrix<double, 3>(E*I,E*I,G*2*I).inverse();
const DiagonalMatrix<double, 3> rhoJ = rho*DiagonalMatrix<double, 3>(I, I, 2*I);
const double rhoA = rho*A;
const Vector3d rhoAg = rho*A*g;
static Vector6d L;

//Hexagonal Stewart-Gough hole pattern
static Vector3d p0[6], r[6];
void initialize_StewartGough_pattern(){
    for(int i = 0; i < 6; i++){
        double theta_b = (-alpha2 + ((i+1) - (i+1)%2)*alpha2 + (i-i%2)*alpha1)/2;
        p0[i] = scrib_R * Vector3d(cos(theta_b), sin(theta_b), 0);
        double theta_L = (-alpha1 + ((i+1) - (i+1)%2)*alpha1 + (i-i%2)*alpha2)/2;
        r[i] = scrib_R * Vector3d(cos(theta_L), sin(theta_L), 0);
    }
}

//ODE describing elastic rod statics
void cosseratRodOde(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Map<Vector3d> n(&y[12]);
    Map<Vector3d> m(&y[15]);

    //Hard-coded constitutive law w/ no precurvature
    Vector3d v = Kse_inv*transposeMultiply(R,n); v(2) += 1;
    Vector3d u = Kbt_inv*transposeMultiply(R,m);

    //Refer to the state vector derivative by its components
    Map<Vector3d> p_s (&y_s_out[0]);
    Map<Matrix3d> R_s (&y_s_out[3]);
    Map<Vector3d> n_s (&y_s_out[12]);
    Map<Vector3d> m_s (&y_s_out[15]);
    Map<Vector3d> q_s (&y_s_out[18]);
    Map<Vector3d> w_s (&y_s_out[21]);

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    n_s = -rhoAg;
    m_s = -p_s.cross(n);
    q_s = Vector3d::Zero();
    w_s = Vector3d::Zero();

    //Output argument for variables with time derivatives
    z_out << Vector6d::Zero(), v, u;
}

//PDE semi-discretized in time describing elastic rod dynamics
void cosseratRodPDE(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y, Map<VectorXd> z_h){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d n = Map<Vector3d>(&y[12]);
    Vector3d m = Map<Vector3d>(&y[15]);
    Vector3d q = Map<Vector3d>(&y[18]);
    Vector3d w = Map<Vector3d>(&y[21]);

    //Unpack history terms
    Map<Vector3d> q_h (&z_h[0]);
    Map<Vector3d> w_h (&z_h[3]);
    Map<Vector3d> v_h (&z_h[6]);
    Map<Vector3d> u_h (&z_h[9]);

    //Material constitutive equation
    Vector3d v = Kse_inv*transposeMultiply(R,n); v(2) += 1;
    Vector3d u = Kbt_inv*transposeMultiply(R,m);
    Vector3d q_t = c0*q + q_h;
    Vector3d w_t = c0*w + w_h;
    Vector3d v_t = c0*v + v_h;
    Vector3d u_t = c0*u + u_h;

    //Pack state vector derivative
    Map<Vector3d> p_s (&y_s_out[0]);
    Map<Matrix3d> R_s (&y_s_out[3]);
    Map<Vector3d> n_s (&y_s_out[12]);
    Map<Vector3d> m_s (&y_s_out[15]);
    Map<Vector3d> q_s (&y_s_out[18]);
    Map<Vector3d> w_s (&y_s_out[21]);

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    n_s = rhoA*R*(q_t + w.cross(q)) - rhoAg;
    m_s = R*(rhoJ*w_t + w.cross(rhoJ*w)) - p_s.cross(n);
    q_s = v_t - u.cross(q) + w.cross(v);
    w_s = u_t - u.cross(w);

    //Output argument for variables with time derivatives
    z_out << q, w, v, u;
}

static MatrixXd Y[6], Z[6], Z_h[6];
static MatrixXd X(24,1), X_h(24,1);
const Map<Vector3d> pE_h(&X_h(0,0));
const Map<Vector3d> vE_h(&X_h(3,0));
const Map<Matrix3d> RE_h(&X_h(6,0));
const Map<Vector3d> wE_h(&X_h(15,0));
const Map<Vector6d> L_h(&X_h(18,0));

template<bool is_forward, bool is_dynamic>
VectorXd obj(VectorXd guess){
    VectorXd residual(42);

    if(is_forward){
        pE = guess.segment<3>(36);
        RE = generateRotation(guess.segment<3>(39));
    }else{
        L = guess.segment<6>(36);
    }

    Matrix3d Jg = RE*plate_mass_moment*RE.transpose();
    Vector3d vE, wE, EF, EM;
    Vector6d L_t;
    if(is_dynamic){
        vE = c0*pE + pE_h;
        Vector3d aE = c0*vE + vE_h;
        wE = inv_hat((c0*RE+RE_h)*RE.transpose());
        Vector3d wE_t = c0*wE + wE_h;
        L_t = c0*L + L_h;

        EF = plate_mass*(g - aE);
        EM = -Jg*wE_t - wE.cross(Jg*wE);
    }else{
        vE = wE = Vector3d::Zero();
        L_t = Vector6d::Zero();

        EF = plate_mass*g;
        EM = Vector3d::Zero();
    }
    X << pE, vE, Map<VectorXd>(RE.data(), 9), wE, L;

    for(int i = 0; i < 6; i++){
        Vector3d n0 = guess.segment<3>(6*i);
        Vector3d m0 = guess.segment<3>(6*i+3);
        double wz0 = 0; //Assume negligible torsional velocity through baseplate
        VectorXd y0(24);
        y0 << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0, 0, 0, L_t(i), 0, 0, wz0;
        Y[i].col(0) = y0;

        //Numerically integrate the Cosserat rod equations
        if(is_dynamic) TimeBdfAlpha_SpaceEuler<cosseratRodPDE,24,12,N>(Y[i],Z[i],y0,L(i),Z_h[i]);
        else euler<cosseratRodOde,24,12,N>(Y[i],Z[i],y0,L(i));

        Vector3d pL_shot = Y[i].block<3,1>(0, N-1);
        Matrix3d RL_shot = Map<Matrix3d>(Y[i].block<9,1>(3, N-1).data());
        Vector3d nL = Y[i].block<3,1>(12, N-1);
        Vector3d mL = Y[i].block<3,1>(15, N-1);

        residual.segment<3>(6*i) = pL_shot - (pE + RE*r[i]);
        residual.segment<3>(6*i+3) = linear_rotation_error(RL_shot,RE);

        EF -= nL;
        EM -= (mL + (RE*r[i]).cross(nL));
    }
    residual.segment<3>(36) = EF;
    residual.segment<3>(39) = EM;

    return residual;
}

int main(int, char**){
    initialize_StewartGough_pattern();

    const bool Kinematics = false;
    const bool Dynamics = true;
    const bool Forward = true;
    const bool Inverse = false;

    //Setup the grid
    for(int i = 0; i < 6; i++){
        Y[i].resize(24,N);
        Z[i].resize(12,N);
        Z_h[i].resize(12,N);
    }

    //Solve the "home" configuration
    VectorXd guess(42);
    guess << VectorXd::Zero(36), pE(2)*VectorXd::Ones(6);
    guess = solveLevenbergMarquardt<obj<Inverse,Kinematics> >(guess);

    //Declare actuation profile and set first entry
    Vector6d leg_lengths[num_actuated+1];
    leg_lengths[0] = guess.segment<6>(36);

    //Inverse kinematics to find the input trajectory
    for(int i = 1; i <= num_actuated; i++){
        std::cout << "IK: " << i << std::endl;

        pE = p_start + (p_final - p_start)*i/static_cast<double>(num_actuated);
        guess = solveLevenbergMarquardt<obj<Inverse,Kinematics> >(guess);
        leg_lengths[i] = guess.segment<6>(36);
    }

    //Reset the robot state and solve ICs
    pE = p_start;
    L = leg_lengths[0];
    guess << VectorXd::Zero(36), p_start, Vector3d::Zero();
    guess = solveLevenbergMarquardt<obj<Forward,Kinematics> >(guess);

    //Prepare for the dynamic problem
    int total_frames = 1 + num_actuated + num_swaying;
    MatrixXd centerlines(6*3*total_frames,N), g_frames(3*total_frames,4);
    VectorXd pEx(total_frames), pEz(total_frames);
    TimeManagerBdfAlpha* time_scheme[7];
    for(int i = 0; i < 6; i++) time_scheme[i] = new TimeManagerBdfAlpha(Z[i], Z_h[i], dt, alpha);
    time_scheme[6] = new TimeManagerBdfAlpha(X, X_h, dt, alpha);

    //Simulate the dynamic behavior
    for(int i = 0; i < total_frames; i++){
        std::cout << "Dynamic: " << i << std::endl;

        if(i <= num_actuated) L = leg_lengths[i];
        guess = solveLevenbergMarquardt<obj<Forward,Dynamics> >(guess);
        for(int j = 0; j < 7; j++) time_scheme[j]->advanceTime();

        pEx(i) = pE.x();
        pEz(i) = pE.z();
        g_frames.block<3,3>(3*i,0) = RE;
        g_frames.block<3,1>(3*i,3) = pE;
        for(int j = 0; j < 6; j++) centerlines.block<3,N>(18*i+3*j,0) = Y[j].block<3,N>(0,0);
    }

    //Save results for Blender visualization
    std::fstream file("centerlines.dat", std::fstream::out);
    file << centerlines;
    file.close();

    file = std::fstream("ee_transformations.dat", std::fstream::out);
    file << g_frames;
    file.close();

    //Show the end-effector trajectory
    #ifdef QT_CORE_LIB
    plot(pEx, pEz, "CSG Dynamic Trajectory", "x (m)", "z (m)");
    #endif

    return 0;
}
