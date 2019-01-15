//Note: Run this with compiler optimizations enabled. It is quite slow otherwise.

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
const Vector3d g = -9.81*Vector3d::UnitZ();
const double ee_mass = 0.1;
const double scrib_R = 0.087; //radius of the hole pattern
const double alpha1 = 100*pi/180; //major angle of the hole pattern

//Dependent parameter calculations
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
const DiagonalMatrix<double, 3> Kse = DiagonalMatrix<double, 3>(G*A,G*A,E*A);
const DiagonalMatrix<double, 3> Kbt = DiagonalMatrix<double, 3>(E*I,E*I,G*J);
const Vector3d F = ee_mass*g;
const Vector3d M = Vector3d::Zero();
const Vector3d pE = 0.4*Vector3d::UnitZ();
const Matrix3d RE = Ry(10*pi/180); //Bend 10 degrees above the y-axis
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
    Vector3d n_s = -rho*A*g;
    Vector3d m_s = -p_s.cross(n);

    //Pack state vector derivative
    VectorXd y_s(18);
    y_s << p_s, Map<VectorXd>(R_s.data(), 9), n_s, m_s;

    return y_s;
}

//Shooting method objective function
const int N = 100;
static std::vector<VectorXd> px(12);
static std::vector<VectorXd> pz(12);
static MatrixXd p(3*6,N);
VectorXd shootingFunction(VectorXd guess){
    VectorXd residual(36);

    Vector3d EF = F;
    Vector3d EM = M;

    for(int i = 0; i < 6; i++){
        Vector3d n0 = guess.segment<3>(5*i);
        Vector2d m0xy = guess.segment<2>(5*i+3);
        double L = guess(30+i);

        VectorXd y0(18);
        y0 << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;

        //Numerically integrate the Cosserat rod equations
        MatrixXd Y = ode4<cosseratRodOde,N>(y0, L);
        px[i] = Y.row(0);
        pz[i] = Y.row(2);
        p.block<3,N>(3*i,0) = Y.block<3,N>(0,0);

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

int main(int, char**){
    initialize_StewartGough_pattern();

    VectorXd init_guess(36);
    init_guess << VectorXd::Zero(30), pE(2)*VectorXd::Ones(6);

    //Solve with shooting method
    solveLevenbergMarquardt<shootingFunction>(init_guess);

    #ifdef QT_CORE_LIB
    //Include lines for the end-effector by connecting the rod ends
    for(auto i = 0; i < 5; i++){
        px[6+i] = Vector2d(px[i](N-1), px[i+1](N-1));
        pz[6+i] = Vector2d(pz[i](N-1), pz[i+1](N-1));
    }
    px[11] = Vector2d(px[5](N-1), px[0](N-1));
    pz[11] = Vector2d(pz[5](N-1), pz[0](N-1));
    //Build the color argument and show the solution
    const Vector3d b = Vector3d(0, 0, 255);
    const Vector3d r = Vector3d(255, 0, 0);
    std::vector<Vector3d> colors = {b,b,b,b,b,b, r,r,r,r,r,r};
    plot(px, pz, colors, "Continuum Stewart Gough BVP Solution", "x (m)", "z (m)");
    #endif

    for(int i = 0; i < 6; i++)
        std::cout << "Rod " << i << " centerline:\n" << p.block<3,N>(3*i,0) << '\n' << std::endl;

    return 0;
}
