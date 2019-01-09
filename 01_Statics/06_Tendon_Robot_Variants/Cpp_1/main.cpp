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

const double tendon_offset = 0.01506;
#define R(theta) (tendon_offset*Vector3d(cos(theta), sin(theta), 0))
const Vector3d r[num_tendons] = {R(0), R(pi/2), R(pi), R(pi*3/2)};
#undef R

const VectorXd tau = Vector4d(15,0,0,0);

//Dependent parameter calculations
const double area = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
const DiagonalMatrix<double, 3> Kbt (E*I,E*I,G*J);
const Matrix3d Kbt_dense = Kbt.toDenseMatrix();

void kirchhoffTendonRobotOde(VectorXd& y_s_out, VectorXd& y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d nb = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);

    Vector3d a(0,0,0);
    Vector3d b(0,0,0);
    Matrix3d G; G.setZero();
    Matrix3d H_plus_Kbt = Kbt_dense;

    for(int i = 0; i < num_tendons; i++){
        Vector3d pb_si = u.cross(r[i]);  pb_si.z() += 1;
        double pb_s_norm = pb_si.norm();
        Matrix3d A_i = -hat_squared(pb_si)*(tau(i)/pow(pb_s_norm,3));
        Matrix3d G_i = -hat_postmultiply(A_i,r[i]);
        Vector3d a_i = A_i*(u.cross(pb_si));

        a += a_i;
        b += r[i].cross(a_i);
        G += G_i;
        H_plus_Kbt += hat_premultiply(r[i],G_i);
    }

    //Pack state vector derivative
    Map<Vector3d> p_s(&y_s_out[0]);
    Map<Matrix3d> R_s(&y_s_out[3]);
    Map<Vector3d> nb_s(&y_s_out[12]);
    Map<Vector3d> u_s(&y_s_out[15]);

    //ODEs
    p_s = R.col(2);
    R_s = hat_postmultiply(R,u);
    u_s = H_plus_Kbt.inverse()*(-u.cross(Kbt*u) - Vector3d::UnitZ().cross(nb) - b);
    nb_s = -u.cross(nb) - G*u_s - a - transposeMultiply(R,rho*area*g);
}

//Boundary conditions
const Vector3d p0 = Vector3d::Zero();
const Matrix3d R0 = Matrix3d::Identity();

static MatrixXd Y; //Declare Y global to save results
VectorXd shootingFunction(VectorXd guess){
    VectorXd y0(18);
    y0 << p0, Map<VectorXd>(Matrix3d(R0).data(), 9), guess;

    //Numerically integrate the Cosserat rod equations
    Y = ode4<kirchhoffTendonRobotOde>(y0, L);

    //Find the internal forces in the backbone prior to the final plate
    Vector3d nb = Y.block<3,1>(12,Y.cols()-1);
    Vector3d uL = Y.block<3,1>(15,Y.cols()-1);
    Vector3d mb = Kbt*uL;

    //Find the equilibrium error at the tip, considering tendon forces
    Vector3d force_error = -nb;
    Vector3d moment_error = -mb;
    for(int i = 0; i < num_tendons; i++){
        Vector3d pb_si = uL.cross(r[i]) + Vector3d::UnitZ();
        Vector3d Fb_i = -tau(i)*pb_si.normalized();
        force_error += Fb_i;
        moment_error += r[i].cross(Fb_i);
    }

    Vector6d distal_error;
    distal_error << force_error, moment_error;

    return distal_error;
}

int main(){
    Vector6d init_guess = Vector6d::Zero(); //nb and u

    //Solve with shooting method
    VectorXd wrench_soln = solveLevenbergMarquardt<shootingFunction>(init_guess);

    //Save results for Blender visualization
    std::fstream file("centerline.dat", std::fstream::out);
    file << Y.block(0,0,3,Y.cols());
    file.close();

    MatrixXd tendonlines(3*num_tendons, Y.cols());
    for(int i = 0; i < Y.cols(); i++){
        Vector3d p = Y.block<3,1>(0,i);
        Matrix3d R = Map<Matrix3d>(&Y(3,i));
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
