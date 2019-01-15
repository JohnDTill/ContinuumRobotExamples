#include <iostream>
#include <commonmath.h>
#include <numericalintegration.h>
#include <convexoptimization.h>
using namespace ContinuumRobotLibrary;

#ifdef QT_CORE_LIB
#include <simpleplotting.h>
#endif

//Independent Parameters
const double E = 200e9;
const double G = 80e9;
const double rad = 0.001;
const double rho = 8000;
const Vector3d g = 9.81*Vector3d::UnitX();
const double L = 0.5;

//Dependent parameter calculations
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
const DiagonalMatrix<double, 3> Kse = DiagonalMatrix<double, 3>(G*A,G*A,E*A);
const DiagonalMatrix<double, 3> Kbt = DiagonalMatrix<double, 3>(E*I,E*I,G*J);

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
    Vector3d n_s = -rho*A*g;
    Vector3d m_s = -p_s.cross(n);

    //Pack state vector derivative
    VectorXd y_s(18);
    y_s << p_s, Map<VectorXd>(R_s.data(), 9), n_s, m_s;

    return y_s;
}

//Boundary conditions
const Vector3d p0 = Vector3d::Zero();
const Matrix3d R0 = Matrix3d::Identity();
const Vector3d pL = Vector3d(0, -0.1*L, 0.8*L);
const Matrix3d RL = Matrix3d::Identity();

static MatrixXd Y; //Declare Y global for shooting and visualization
VectorXd shootingFunction(VectorXd guess){
    VectorXd y0(18);
    y0 << p0, Map<VectorXd>(Matrix3d(R0).data(), 9), guess;

    //Numerically integrate the Cosserat rod equations
    Y = ode4<cosseratRodOde>(y0, L);

    //Calculate the pose (position and orientation) error at the far end
    Vector3d pL_shooting = Y.block<3,1>(0, Y.cols()-1);
    Matrix3d RL_shooting = Map<Matrix3d>(Y.block<9,1>(3, Y.cols()-1).data());

    Vector6d distal_error;
    distal_error << pL - pL_shooting, rotation_error(RL,RL_shooting);

    return distal_error;
}

int main(int, char**){
    Vector6d init_guess = Vector6d::Zero();

    //Solve with shooting method
    VectorXd wrench_soln = solveLevenbergMarquardt<shootingFunction>(init_guess);

    #ifdef QT_CORE_LIB
    plot(Y.row(1), Y.row(2), "Cosserat Rod BVP Solution", "y (m)", "z (m)");
    #endif

    std::cout << "Solved rod centerline:\n" << Y.block(0,0,3,Y.cols()) << std::endl;

    return 0;
}
