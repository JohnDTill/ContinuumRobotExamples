#include <iostream>
#include <commonmath.h>
#include <numericalintegration.h>
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

int main(int, char**){
    //Set initial conditions
    Vector3d p0 = Vector3d::Zero();
    Matrix3d R0 = Matrix3d::Identity();
    Vector3d n0 = Vector3d::UnitY();
    Vector3d m0 = Vector3d::Zero();

    VectorXd y0(18);
    y0 << p0, Map<VectorXd>(R0.data(), 9), n0, m0;

    //Numerically integrate the Cosserat rod equations
    MatrixXd Y = ode4<cosseratRodOde>(y0, L);

    #ifdef QT_CORE_LIB //Plot the solution if Qt is used
    plot(Y.row(1), Y.row(2), "Cosserat Rod IVP Solution", "y (m)", "z (m)");
    #endif

    //Output the data
    std::cout << "Solved rod centerline:\n" << Y.block(0,0,3,Y.cols()) << std::endl;

    return 0;
}
