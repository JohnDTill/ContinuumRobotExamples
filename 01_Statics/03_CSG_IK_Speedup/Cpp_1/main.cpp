#include <stdio.h>
#include <ctime>
#include <iostream>
#include <commonmath.h>
#include <numericalintegration.h>
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
const DiagonalMatrix<double, 3> Kse = DiagonalMatrix<double, 3>(G*A,G*A,E*A);
const DiagonalMatrix<double, 3> Kbt = DiagonalMatrix<double, 3>(E*I,E*I,G*J);
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
VectorXd cosseratRodOde(VectorXd y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(y.segment<9>(3).data());
    Vector3d n = y.segment<3>(12);
    Vector3d m = y.segment<3>(15);

    //Hard-coded material constitutive equation w/ no precurvature
    Vector3d v = Kse.inverse()*R.transpose()*n + Vector3d::UnitZ();
    Vector3d u = Kbt.inverse()*R.transpose()*m;

    //ODEs
    Vector3d p_s = R*v;
    Matrix3d R_s = R*hat(u);
    Vector3d n_s = Vector3d::Zero();
    Vector3d m_s = -p_s.cross(n);

    //Pack state vector derivative
    VectorXd y_s(18);
    y_s << p_s, Map<VectorXd>(R_s.data(), 9), n_s, m_s;

    return y_s;
}

//Shooting method objective function
const int N = 40;
VectorXd shootingFunction(VectorXd guess){
    VectorXd residual(36);

    Vector3d EF = Vector3d::Zero();
    Vector3d EM = Vector3d::Zero();

    for(int i = 0; i < 6; i++){
        Vector3d n0 = guess.segment<3>(5*i);
        Vector2d m0xy = guess.segment<2>(5*i+3);
        double L = guess(30+i);

        VectorXd y0(18);
        y0 << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;

        //Numerically integrate the Cosserat rod equations
        MatrixXd Y = ode4<cosseratRodOde,N>(y0, L);

        Vector3d pL_shot = Y.block<3,1>(0, N-1);
        Matrix3d RL_shot = Map<Matrix3d>(Y.block<9,1>(3, N-1).data());
        Vector3d nL = Y.block<3,1>(12, N-1);
        Vector3d mL = Y.block<3,1>(15, N-1);

        residual.segment<3>(5*i) = pL_shot - (pE + RE*r[i]);
        residual.segment<2>(5*i+3) =
                linear_rotation_error(RL_shot,RE).segment<2>(0);

        EF -= nL;
        EM -= (mL + (RE*r[i]).cross(nL));
    }

    residual.segment<3>(30) = EF;
    residual.segment<3>(33) = EM;

    return residual;
}

int main(){
    initialize_StewartGough_pattern();

    VectorXd guess(36);
    guess << VectorXd::Zero(30), pE(2)*VectorXd::Ones(6);

    //Solve with shooting method
    guess = solveLevenbergMarquardt<shootingFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    pE = Vector3d(0, 0.02, 0.48);
    guess = solveLevenbergMarquardt<shootingFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    //Benchmark speed test
    std::clock_t start = std::clock();
    const int M = 300;
    for(int i = 0; i < M; i++){
        if(i%200 < 100){
            pE(1) += 0.001;
            pE(2) += 0.001;
        }else{
            pE(1) -= 0.001;
            pE(2) -= 0.001;
        }
        guess = solveLevenbergMarquardt<shootingFunction>(guess, 1e-6, 500, 1e-4, 0.5);
    }
    double duration = ( std::clock() - start ) / double(CLOCKS_PER_SEC);
    std::cout << "Time elasped: " << duration << "s" << std::endl;
    std::cout << "Solution rate: " << M/duration << "Hz" << std::endl;

    return 0;
}
