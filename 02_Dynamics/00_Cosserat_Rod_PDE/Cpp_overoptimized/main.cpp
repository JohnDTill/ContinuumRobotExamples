#include <iostream>
#include <stdio.h>
#include <ctime>
#include <commonmath.h>
#include <convexoptimization.h>
#include "timemanagement.h"
#include "numericalpde.h"
using namespace ContinuumRobotLibrary;

#ifdef QT_CORE_LIB
#include <simpleanimation.h>
#endif

//Independent Parameters
const double E = 200e9;
const double G = 80e9;
const double rad = 0.001;
const double rho = 8000;
const Vector3d g = -9.81*Vector3d::UnitX();
const double L = 0.5;
const int N = 40; //Number of spatial points
const double dt = 2e-3; //Time step
const double alpha = -0.2; //Parameter for implicit time discretization
const double T = 25; //Total simulation time
const double initial_tip_mass = 0.02; //Weight for ICs, released at t0+
const double Bs = 0; //Shearing material damping
const double Be = 0; //Extension material damping
const double Bb = 0; //Bending material damping
const double Bt = 0; //Torsional material damping

//Dependent parameter calculations
const int Nm1 = N-1;
const double ds = L/Nm1;
static Vector3d F = initial_tip_mass*g;
const double c0 = TimeManagerBdfAlpha::getC0(dt,alpha);
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const DiagonalMatrix<double, 3> Kse = DiagonalMatrix<double, 3>(G*A,G*A,E*A);
const Vector3d Kse_e3(0, 0, E*A);
const double EA = E*A;
const DiagonalMatrix<double, 3> Kbt = DiagonalMatrix<double, 3>(E*I,E*I,G*2*I);
const DiagonalMatrix<double, 3> Kse_inv = DiagonalMatrix<double, 3>(G*A,G*A,E*A).inverse();
const DiagonalMatrix<double, 3> Kbt_inv = DiagonalMatrix<double, 3>(E*I,E*I,G*2*I).inverse();
const DiagonalMatrix<double, 3> rhoJ = rho*DiagonalMatrix<double, 3>(I, I, 2*I);
const DiagonalMatrix<double, 3> Bse(Bs, Bs, Be);
const DiagonalMatrix<double, 3> Bbt(Bb, Bb, Bt);
const double Ks_plus_c0_Bs_inv = 1/(G*A + c0*Bs);
const double Ke_plus_c0_Be_inv = 1/(E*A + c0*Be);
const double Kb_plus_c0_Bb_inv = 1/(E*I + c0*Bb);
const double Kt_plus_c0_Bt_inv = 1/(G*2*I + c0*Bt);
const double rhoA = rho*A;
const Vector3d rhoAg = rho*A*g;
const double rhoI = rho*I;
const double rho2I = rho*2*I;

static MatrixXd Y(24,N); //Declare Y global for shooting and plotting
static MatrixXd Z(12,N);
static MatrixXd Z_h(12,N);

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

//Integration routine for the PDE
void cosseratRodPDE_Euler(){
    for(int i = 0; i < Nm1; i++){
        int state_offset = 24*i;
        int deriv_offset = 12*i;
        int state_next = 24*(i+1);

        #define p0 Y(state_offset)
        #define p1 Y(state_offset+1)
        #define p2 Y(state_offset+2)

        #define R00 Y(state_offset+3)
        #define R10 Y(state_offset+4)
        #define R20 Y(state_offset+5)
        #define R01 Y(state_offset+6)
        #define R11 Y(state_offset+7)
        #define R21 Y(state_offset+8)
        #define R02 Y(state_offset+9)
        #define R12 Y(state_offset+10)
        #define R22 Y(state_offset+11)

        #define nl0 Y(state_offset+12)
        #define nl1 Y(state_offset+13)
        #define nl2 Y(state_offset+14)

        #define ml0 Y(state_offset+15)
        #define ml1 Y(state_offset+16)
        #define ml2 Y(state_offset+17)

        #define q0 Y(state_offset+18)
        #define q1 Y(state_offset+19)
        #define q2 Y(state_offset+20)

        #define w0 Y(state_offset+21)
        #define w1 Y(state_offset+22)
        #define w2 Y(state_offset+23)

        #define v0 Z(deriv_offset+6)
        #define v1 Z(deriv_offset+7)
        #define v2 Z(deriv_offset+8)

        #define u0 Z(deriv_offset+9)
        #define u1 Z(deriv_offset+10)
        #define u2 Z(deriv_offset+11)

        #define q0_h Z_h(deriv_offset)
        #define q1_h Z_h(deriv_offset+1)
        #define q2_h Z_h(deriv_offset+2)

        #define w0_h Z_h(deriv_offset+3)
        #define w1_h Z_h(deriv_offset+4)
        #define w2_h Z_h(deriv_offset+5)

        #define v0_h Z_h(deriv_offset+6)
        #define v1_h Z_h(deriv_offset+7)
        #define v2_h Z_h(deriv_offset+8)

        #define u0_h Z_h(deriv_offset+9)
        #define u1_h Z_h(deriv_offset+10)
        #define u2_h Z_h(deriv_offset+11)

        v0 = Ks_plus_c0_Bs_inv * (nl0 - Bs*v0_h);
        v1 = Ks_plus_c0_Bs_inv * (nl1 - Bs*v1_h);
        v2 = Ke_plus_c0_Be_inv * (nl2 + EA - Be*v2_h);

        u0 = Kb_plus_c0_Bb_inv * (ml0 - Bb*u0_h);
        u1 = Kb_plus_c0_Bb_inv * (ml1 - Bb*u1_h);
        u2 = Kt_plus_c0_Bt_inv * (ml2 - Bt*u2_h);

        double wt0 = c0*w0 + w0_h;
        double wt1 = c0*w1 + w1_h;
        double wt2 = c0*w2 + w2_h;

        double p0s = R00*v0 + R01*v1 + R02*v2;
        double p1s = R10*v0 + R11*v1 + R12*v2;
        double p2s = R20*v0 + R21*v1 + R22*v2;

        double R00s = R01*u2 - R02*u1;
        double R10s = R11*u2 - R12*u1;
        double R20s = R21*u2 - R22*u1;
        double R01s = R02*u0 - R00*u2;
        double R11s = R12*u0 - R10*u2;
        double R21s = R22*u0 - R20*u2;

        double nl0s = u2*nl1 - u1*nl2 + rhoA*( c0*q0 + q0_h - q1*w2 + q2*w1 ) - R00*rhoAg(0) - R01*rhoAg(1) - R02*rhoAg(2);
        double nl1s = u0*nl2 - u2*nl0 + rhoA*( c0*q1 + q1_h + q0*w2 - q2*w0 ) - R10*rhoAg(0) - R11*rhoAg(1) - R12*rhoAg(2);
        double nl2s = u1*nl0 - u0*nl1 + rhoA*( c0*q2 + q2_h - q0*w1 + q1*w0 ) - R20*rhoAg(0) - R21*rhoAg(1) - R22*rhoAg(2);

        double ml0s = u2*ml1 - u1*ml2 + rhoI*(wt0 + w1*w2) + v2*nl1 - v1*nl2;
        double ml1s = u0*ml2 - u2*ml0 + rhoI*(wt1 - w0*w2) + v0*nl2 - v2*nl0;
        double ml2s = u1*ml0 - u0*ml1 + rho2I*wt2          + v1*nl0 - v0*nl1;

        double qs0 = c0*v0 + v0_h - q2*u1 + q1*u2 + v2*w1 - v1*w2;
        double qs1 = c0*v1 + v1_h + q2*u0 - q0*u2 + v0*w2 - v2*w0;
        double qs2 = c0*v2 + v2_h + q0*u1 - q1*u0 + v1*w0 - v0*w1;

        double ws0 = c0*u0 + u0_h - u1*w2 + u2*w1;
        double ws1 = c0*u1 + u1_h - u2*w0 + u0*w2;
        double ws2 = c0*u2 + u2_h - u0*w1 + u1*w0;

        #define p0_next Y(state_next)
        #define p1_next Y(state_next+1)
        #define p2_next Y(state_next+2)

        #define R00_next Y(state_next+3)
        #define R10_next Y(state_next+4)
        #define R20_next Y(state_next+5)
        #define R01_next Y(state_next+6)
        #define R11_next Y(state_next+7)
        #define R21_next Y(state_next+8)
        #define R02_next Y(state_next+9)
        #define R12_next Y(state_next+10)
        #define R22_next Y(state_next+11)

        #define nl0_next Y(state_next+12)
        #define nl1_next Y(state_next+13)
        #define nl2_next Y(state_next+14)

        #define ml0_next Y(state_next+15)
        #define ml1_next Y(state_next+16)
        #define ml2_next Y(state_next+17)

        #define q0_next Y(state_next+18)
        #define q1_next Y(state_next+19)
        #define q2_next Y(state_next+20)

        #define w0_next Y(state_next+21)
        #define w1_next Y(state_next+22)
        #define w2_next Y(state_next+23)

        p0_next = p0 + ds*p0s;
        p1_next = p1 + ds*p1s;
        p2_next = p2 + ds*p2s;

        R00_next = R00 + ds*R00s;
        R10_next = R10 + ds*R10s;
        R20_next = R20 + ds*R20s;
        R01_next = R01 + ds*R01s;
        R11_next = R11 + ds*R11s;
        R21_next = R21 + ds*R21s;

        R02_next = R10_next*R21_next - R11_next*R20_next;
        R12_next = R20_next*R01_next - R21_next*R00_next;
        R22_next = R00_next*R11_next - R01_next*R10_next;

        nl0_next = nl0 + ds*nl0s;
        nl1_next = nl1 + ds*nl1s;
        nl2_next = nl2 + ds*nl2s;

        ml0_next = ml0 + ds*ml0s;
        ml1_next = ml1 + ds*ml1s;
        ml2_next = ml2 + ds*ml2s;

        q0_next = q0 + ds*qs0;
        q1_next = q1 + ds*qs1;
        q2_next = q2 + ds*qs2;

        w0_next = w0 + ds*ws0;
        w1_next = w1 + ds*ws1;
        w2_next = w2 + ds*ws2;

        #undef p0
        #undef p1
        #undef p2

        #undef R00
        #undef R10
        #undef R20
        #undef R01
        #undef R11
        #undef R21
        #undef R02
        #undef R12
        #undef R22

        #undef nl0
        #undef nl1
        #undef nl2

        #undef ml0
        #undef ml1
        #undef ml2

        #undef q0
        #undef q1
        #undef q2

        #undef w0
        #undef w1
        #undef w2

        #undef p0_next
        #undef p1_next
        #undef p2_next

        #undef R00_next
        #undef R10_next
        #undef R20_next
        #undef R01_next
        #undef R11_next
        #undef R21_next
        #undef R02_next
        #undef R12_next
        #undef R22_next

        #undef nl0_next
        #undef nl1_next
        #undef nl2_next

        #undef ml0_next
        #undef ml1_next
        #undef ml2_next

        #undef q0_next
        #undef q1_next
        #undef q2_next

        #undef w0_next
        #undef w1_next
        #undef w2_next

        #undef v0
        #undef v1
        #undef v2
        #undef u0
        #undef u1
        #undef u2

        #undef q0_h
        #undef q1_h
        #undef q2_h
        #undef w0_h
        #undef w1_h
        #undef w2_h

        #undef v0_h
        #undef v1_h
        #undef v2_h
        #undef u0_h
        #undef u1_h
        #undef u2_h
    }
}

//Shooting method objective function
template<bool is_dynamic>
VectorXd objFunc(VectorXd guess){
    VectorXd y0(24);
    y0 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, guess, Vector6d::Zero();
    Y.col(0) = y0;

    //Numerically integrate the semi-discretized rod PDEs
    if(is_dynamic) cosseratRodPDE_Euler();
    else           euler<cosseratRodOde,24,12,N>(Y,Z,y0,L);

    Vector3d nL_shooting = Y.block<3,1>(12, N-1);
    Vector3d mL_shooting = Y.block<3,1>(15, N-1);

    Vector6d distal_error;
    if(is_dynamic) distal_error << nL_shooting, mL_shooting;
    else           distal_error << nL_shooting - F, mL_shooting;

    return distal_error;
}

int main(){
    //Solve initial conditions
    Vector6d guess = Vector6d::Zero();
    Vector6d reaction_wrench = solveLevenbergMarquardt<objFunc<false> >(guess);
    TimeManagerBdfAlpha time_scheme(Z, Z_h, dt, alpha);
    F = Vector3d::Zero(); //The tip load is removed

    int M = static_cast<int>(T/dt);
    MatrixXd px(M,N), pz(M,N);
    std::clock_t start = std::clock();
    for(int i = 0; i < M; i++){
        px.row(i) = Y.row(0);
        pz.row(i) = Y.row(2);
        time_scheme.advanceTime(); //Update Z_h from history of Z
        //Solve dynamic problem
        reaction_wrench = solveLevenbergMarquardt<objFunc<true> >(reaction_wrench);
        Z.block<6,N>(0,0) = Y.block<6,N>(18,0); //This part of Z isn't set during integration
    }
    double duration = ( std::clock() - start ) / static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "Time elasped: " << duration << "s" << std::endl;
    std::cout << "Real-time ratio: " << T/duration << std::endl;

    #ifdef QT_CORE_LIB
    playAnimation(pz, px, dt, "Dynamic Cantilever", "z (m)", "x (m)");
    #endif

    return 0;
}
