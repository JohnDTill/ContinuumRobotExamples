#include <stdio.h>
#include <ctime>
#include <iostream>
#include <commonmath.h>
#include <convexoptimization.h>
using namespace ContinuumRobotLibrary;

//Independent Parameters
const double E = 207e9;
const double rad = 0.0013/2;
const double scrib_R = 0.087; //radius of the hole pattern
const double alpha1 = 100*pi/180; //major angle of the hole pattern

//Dependent parameter calculations
const double G = E/(2*(1+0.305));
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
const double GA_inv = 1/(G*A);
const double EA_inv = 1/(E*A);
const double EI_inv = 1/(E*I);
const DiagonalMatrix<double, 3> Kse_inv = DiagonalMatrix<double, 3>(G*A,G*A,E*A).inverse();
const DiagonalMatrix<double, 3> Kbt_inv = DiagonalMatrix<double, 3>(E*I,E*I,G*J).inverse();
static Vector3d pE = 0.4*Vector3d::UnitZ();
const Matrix3d RE = Matrix3d::Identity();
const double alpha2 = 120*pi/180 - alpha1; //minor angle of the hole pattern

//Hexagonal Stewart-Gough hole pattern
static Vector3d p0[6];
static Vector3d r[6];
void initialize_StewartGough_pattern(){
    for(int i = 0; i < 6; i++){
        double theta_b = (-alpha2 + ((i+1) - (i+1)%2)*alpha2 + (i-i%2)*alpha1)/2;
        p0[i] = scrib_R * Vector3d(cos(theta_b), sin(theta_b), 0);

        double theta_L = (-alpha1 + ((i+1) - (i+1)%2)*alpha1 + (i-i%2)*alpha2)/2;
        r[i] = scrib_R * Vector3d(cos(theta_L), sin(theta_L), 0);
    }
}

//Ordinary differential equation describing elastic rod
void cosseratRodOde(VectorXd& y_s_out, VectorXd& y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Map<Vector3d> n(&y[12]);
    Map<Vector3d> m(&y[15]);

    //Hard-coded material constitutive equation w/ no precurvature
    Vector3d v = Kse_inv*transposeMultiply(R,n); v(2) += 1;
    Vector3d u = Kbt_inv*transposeMultiply(R,m);

    //Refer to the state vector derivative by its components
    Map<Vector3d> p_s (&y_s_out[0]);
    Map<Matrix3d> R_s (&y_s_out[3]);
    Map<Vector3d> n_s (&y_s_out[12]);
    Map<Vector3d> m_s (&y_s_out[15]);

    //ODEs
    p_s = R*v;
    R_s = hat_postmultiply(R,u);
    n_s = Vector3d::Zero();
    m_s = -p_s.cross(n);
}

//Hard-coded numerical integration routine for a Cosserat rod
const int N = 40;
const int Nm1 = N-1;
void ode4_CosseratRodOde(VectorXd& y_out, double L){
    double ds = L/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    for(int i = Nm1; i > 0; i--){
        #define p0 y_out[0]
        #define p1 y_out[1]
        #define p2 y_out[2]

        #define R00 y_out[3]
        #define R10 y_out[4]
        #define R20 y_out[5]
        #define R01 y_out[6]
        #define R11 y_out[7]
        #define R21 y_out[8]
        #define R02 y_out[9]
        #define R12 y_out[10]
        #define R22 y_out[11]

        #define n0 y_out[12]
        #define n1 y_out[13]
        #define n2 y_out[14]

        #define m0 y_out[15]
        #define m1 y_out[16]
        #define m2 y_out[17]

        double v0 = GA_inv*(R00*n0 + R10*n1 + R20*n2);
        double v1 = GA_inv*(R01*n0 + R11*n1 + R21*n2);
        double v2 = EA_inv*(R02*n0 + R12*n1 + R22*n2) + 1;

        double u0 = EI_inv*(R00*m0 + R10*m1 + R20*m2);
        double u1 = EI_inv*(R01*m0 + R11*m1 + R21*m2);

        double ps0 = R00*v0 + R01*v1 + R02*v2;
        double ps1 = R10*v0 + R11*v1 + R12*v2;
        double ps2 = R20*v0 + R21*v1 + R22*v2;

        double Rs00 = - R02*u1;
        double Rs10 = - R12*u1;
        double Rs20 = - R22*u1;
        double Rs01 = R02*u0;
        double Rs11 = R12*u0;
        double Rs21 = R22*u0;

        double ms0 = ps2*n1 - ps1*n2;
        double ms1 = ps0*n2 - ps2*n0;
        double ms2 = ps1*n0 - ps0*n1;

        double y2R00 = R00 + half_ds*Rs00;
        double y2R10 = R10 + half_ds*Rs10;
        double y2R20 = R20 + half_ds*Rs20;
        double y2R01 = R01 + half_ds*Rs01;
        double y2R11 = R11 + half_ds*Rs11;
        double y2R21 = R21 + half_ds*Rs21;
        double y2R02 = y2R10*y2R21 - y2R11*y2R20;
        double y2R12 = y2R20*y2R01 - y2R21*y2R00;
        double y2R22 = y2R00*y2R11 - y2R01*y2R10;

        double y2m0 = m0 + half_ds*ms0;
        double y2m1 = m1 + half_ds*ms1;
        double y2m2 = m2 + half_ds*ms2;

        double y2v0 = GA_inv*(y2R00*n0 + y2R10*n1 + y2R20*n2);
        double y2v1 = GA_inv*(y2R01*n0 + y2R11*n1 + y2R21*n2);
        double y2v2 = EA_inv*(y2R02*n0 + y2R12*n1 + y2R22*n2) + 1;

        double y2u0 = EI_inv*(y2R00*y2m0 + y2R10*y2m1 + y2R20*y2m2);
        double y2u1 = EI_inv*(y2R01*y2m0 + y2R11*y2m1 + y2R21*y2m2);

        double k2ps0 = y2R00*y2v0 + y2R01*y2v1 + y2R02*y2v2;
        double k2ps1 = y2R10*y2v0 + y2R11*y2v1 + y2R12*y2v2;
        double k2ps2 = y2R20*y2v0 + y2R21*y2v1 + y2R22*y2v2;

        double k2Rs00 = - y2R02*y2u1;
        double k2Rs10 = - y2R12*y2u1;
        double k2Rs20 = - y2R22*y2u1;
        double k2Rs01 = y2R02*y2u0;
        double k2Rs11 = y2R12*y2u0;
        double k2Rs21 = y2R22*y2u0;

        double k2ms0 = k2ps2*n1 - k2ps1*n2;
        double k2ms1 = k2ps0*n2 - k2ps2*n0;
        double k2ms2 = k2ps1*n0 - k2ps0*n1;

        double y3R00 = R00 + half_ds*k2Rs00;
        double y3R10 = R10 + half_ds*k2Rs10;
        double y3R20 = R20 + half_ds*k2Rs20;
        double y3R01 = R01 + half_ds*k2Rs01;
        double y3R11 = R11 + half_ds*k2Rs11;
        double y3R21 = R21 + half_ds*k2Rs21;
        double y3R02 = y3R10*y3R21 - y3R11*y3R20;
        double y3R12 = y3R20*y3R01 - y3R21*y3R00;
        double y3R22 = y3R00*y3R11 - y3R01*y3R10;

        double y3m0 = m0 + half_ds*k2ms0;
        double y3m1 = m1 + half_ds*k2ms1;
        double y3m2 = m2 + half_ds*k2ms2;

        double y3v0 = GA_inv*(y3R00*n0 + y3R10*n1 + y3R20*n2);
        double y3v1 = GA_inv*(y3R01*n0 + y3R11*n1 + y3R21*n2);
        double y3v2 = EA_inv*(y3R02*n0 + y3R12*n1 + y3R22*n2) + 1;

        double y3u0 = EI_inv*(y3R00*y3m0 + y3R10*y3m1 + y3R20*y3m2);
        double y3u1 = EI_inv*(y3R01*y3m0 + y3R11*y3m1 + y3R21*y3m2);

        double k3ps0 = y3R00*y3v0 + y3R01*y3v1 + y3R02*y3v2;
        double k3ps1 = y3R10*y3v0 + y3R11*y3v1 + y3R12*y3v2;
        double k3ps2 = y3R20*y3v0 + y3R21*y3v1 + y3R22*y3v2;

        double k3Rs00 = - y3R02*y3u1;
        double k3Rs10 = - y3R12*y3u1;
        double k3Rs20 = - y3R22*y3u1;
        double k3Rs01 = y3R02*y3u0;
        double k3Rs11 = y3R12*y3u0;
        double k3Rs21 = y3R22*y3u0;

        double k3ms0 = k3ps2*n1 - k3ps1*n2;
        double k3ms1 = k3ps0*n2 - k3ps2*n0;
        double k3ms2 = k3ps1*n0 - k3ps0*n1;

        double y4R00 = R00 + ds*k3Rs00;
        double y4R10 = R10 + ds*k3Rs10;
        double y4R20 = R20 + ds*k3Rs20;
        double y4R01 = R01 + ds*k3Rs01;
        double y4R11 = R11 + ds*k3Rs11;
        double y4R21 = R21 + ds*k3Rs21;
        double y4R02 = y4R10*y4R21 - y4R11*y4R20;
        double y4R12 = y4R20*y4R01 - y4R21*y4R00;
        double y4R22 = y4R00*y4R11 - y4R01*y4R10;

        double y4m0 = m0 + ds*k3ms0;
        double y4m1 = m1 + ds*k3ms1;
        double y4m2 = m2 + ds*k3ms2;

        double y4v0 = GA_inv*(y4R00*n0 + y4R10*n1 + y4R20*n2);
        double y4v1 = GA_inv*(y4R01*n0 + y4R11*n1 + y4R21*n2);
        double y4v2 = EA_inv*(y4R02*n0 + y4R12*n1 + y4R22*n2) + 1;

        double y4u0 = EI_inv*(y4R00*y4m0 + y4R10*y4m1 + y4R20*y4m2);
        double y4u1 = EI_inv*(y4R01*y4m0 + y4R11*y4m1 + y4R21*y4m2);

        double k4ps0 = y4R00*y4v0 + y4R01*y4v1 + y4R02*y4v2;
        double k4ps1 = y4R10*y4v0 + y4R11*y4v1 + y4R12*y4v2;
        double k4ps2 = y4R20*y4v0 + y4R21*y4v1 + y4R22*y4v2;

        double k4Rs00 = - y4R02*y4u1;
        double k4Rs10 = - y4R12*y4u1;
        double k4Rs20 = - y4R22*y4u1;
        double k4Rs01 = y4R02*y4u0;
        double k4Rs11 = y4R12*y4u0;
        double k4Rs21 = y4R22*y4u0;

        double k4ms0 = k4ps2*n1 - k4ps1*n2;
        double k4ms1 = k4ps0*n2 - k4ps2*n0;
        double k4ms2 = k4ps1*n0 - k4ps0*n1;

        p0 += sixth_ds*( ps0 + 2*( k2ps0 + k3ps0 ) + k4ps0 );
        p1 += sixth_ds*( ps1 + 2*( k2ps1 + k3ps1 ) + k4ps1 );
        p2 += sixth_ds*( ps2 + 2*( k2ps2 + k3ps2 ) + k4ps2 );

        R00 += sixth_ds*( Rs00 + 2*( k2Rs00 + k3Rs00 ) + k4Rs00 );
        R10 += sixth_ds*( Rs10 + 2*( k2Rs10 + k3Rs10 ) + k4Rs10 );
        R20 += sixth_ds*( Rs20 + 2*( k2Rs20 + k3Rs20 ) + k4Rs20 );
        R01 += sixth_ds*( Rs01 + 2*( k2Rs01 + k3Rs01 ) + k4Rs01 );
        R11 += sixth_ds*( Rs11 + 2*( k2Rs11 + k3Rs11 ) + k4Rs11 );
        R21 += sixth_ds*( Rs21 + 2*( k2Rs21 + k3Rs21 ) + k4Rs21 );

        m0 += sixth_ds*( ms0 + 2*( k2ms0 + k3ms0 ) + k4ms0 );
        m1 += sixth_ds*( ms1 + 2*( k2ms1 + k3ms1 ) + k4ms1 );
        m2 += sixth_ds*( ms2 + 2*( k2ms2 + k3ms2 ) + k4ms2 );

        R02 = R10*R21 - R11*R20;
        R12 = R20*R01 - R21*R00;
        R22 = R00*R11 - R01*R10;

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

        #undef n0
        #undef n1
        #undef n2

        #undef m0
        #undef m1
        #undef m2
    }
}

//Shooting method objective function
static bool integrate = true;
static Vector3d pL_shot[6];
static Matrix3d RL_shot[6];
static Vector3d nL[6];
static Vector3d mL[6];
static VectorXd residual(36);
static VectorXd y(18);
VectorXd shootingFunction(VectorXd guess){

    Vector3d EF = Vector3d::Zero();
    Vector3d EM = Vector3d::Zero();

    for(int i = 0; i < 6; i++){
        if(integrate){
            Vector3d n0 = guess.segment<3>(5*i);
            Vector2d m0xy = guess.segment<2>(5*i+3);
            double L = guess(30+i);

            y << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;

            //Numerically integrate the Cosserat rod equations
            ode4_CosseratRodOde(y, L);

            pL_shot[i] = y.segment<3>(0);
            RL_shot[i] = Map<Matrix3d>(y.segment<9>(3).data());
            nL[i] = y.segment<3>(12);
            mL[i] = y.segment<3>(15);
        }

        residual.segment<3>(5*i) = pL_shot[i] - (pE + RE*r[i]);
        residual.segment<2>(5*i+3) =
                linear_rotation_error(RL_shot[i],RE).segment<2>(0);

        EF -= nL[i];
        EM -= (mL[i] + (RE*r[i]).cross(nL[i]));
    }
    integrate = true;

    residual.segment<3>(30) = EF;
    residual.segment<3>(33) = EM;

    return residual;
}

static VectorXd y_s(18);
const double incr = 1e-8;
void jacobianFunction(MatrixXd& J_out, VectorXd& guess, VectorXd&){
    for(int i = 0; i < 6; i++){
        double L = guess(30+i);
        for(int j = 0; j < 5; j++){
            //A guessed variable is incremented for the finite difference
            double temp = guess(5*i+j);
            guess(5*i+j) += incr;
            Vector3d n0 = Map<Vector3d>(&guess[5*i]);
            Vector2d m0xy = Map<Vector2d>(&guess[5*i+3]);
            guess(5*i+j) = temp;

            y << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;

            //Numerically integrate the Cosserat rod equations
            ode4_CosseratRodOde(y, L);
            Vector3d pL_incr = y.segment<3>(0);
            Matrix3d RL_incr = Map<Matrix3d>(&y(3));
            Vector3d nL_incr = y.segment<3>(12);
            Vector3d mL_incr = y.segment<3>(15);

            //Jacobian blocks are partial derivatives of objective function equations
            J_out.block<3,1>(5*i, 5*i+j) = (pL_incr - pL_shot[i]) / incr;
            J_out.block<2,1>(5*i+3, 5*i+j) = linear_rotation_error( (RL_incr - RL_shot[i])/incr, RE).segment<2>(0);
            Vector3d nL_partial = (nL_incr - nL[i])/incr;
            J_out.block<3,1>(30, 5*i+j) = -nL_partial;
            J_out.block<3,1>(33, 5*i+j) = -( (mL_incr - mL[i]) / incr  + (RE*r[i]).cross(nL_partial) );
        }

        //Partial derivatives w.r.t. arc length do not require integration
        y << pL_shot[i], Map<VectorXd>(RL_shot[i].data(), 9), nL[i], mL[i];
        cosseratRodOde(y_s, y);

        J_out.block<3,1>(5*i, 30+i) = y_s.segment<3>(0);
        J_out.block<2,1>(5*i+3, 30+i) = linear_rotation_error( Map<Matrix3d>(&y_s[3]), RE).segment<2>(0);

        J_out.block<3,1>(30, 30+i) = -y_s.segment<3>(12);
        J_out.block<3,1>(33, 30+i) = -( y_s.segment<3>(15)  + (RE*r[i]).cross(y_s.segment<3>(12)));
    }
}

int main(){
    initialize_StewartGough_pattern();

    VectorXd guess(36);
    guess << VectorXd::Zero(30), pE(2)*VectorXd::Ones(6);

    //Solve with shooting method
    guess = solveLevenbergMarquardt<shootingFunction, jacobianFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    pE = Vector3d(0, 0.02, 0.48);
    guess = solveLevenbergMarquardt<shootingFunction, jacobianFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    //Benchmark speed test
    std::clock_t start = std::clock();
    const int M = 50000;
    for(int i = 0; i < M; i++){
        integrate = false;
        if(i%200 < 100){
            pE(1) += 0.001;
            pE(2) += 0.001;
        }else{
            pE(1) -= 0.001;
            pE(2) -= 0.001;
        }
        guess = solveLevenbergMarquardt<shootingFunction, jacobianFunction>(guess, 1e-6, 500, 1e-4, 0.5);
    }
    double duration = ( std::clock() - start ) / double(CLOCKS_PER_SEC);
    std::cout << "Time elasped: " << duration << "s" << std::endl;
    std::cout << "Solution rate: " << M/duration << "Hz" << std::endl;

    return 0;
}
