#include <stdio.h>
#include <ctime>
#include <iostream>
#include <commonmath.h>
#include <numericalintegration.h>
#include <convexoptimization.h>
using namespace ContinuumRobotLibrary;
#include <thread>
#include <atomic>

//Independent Parameters
const int num_threads = 2;
const double E = 207e9;
const double rad = 0.0013/2;
const double scrib_R = 0.087; //radius of the hole pattern
const double alpha1 = 100*pi/180; //major angle of the hole pattern

//Dependent parameter calculations
const double G = E/(2*(1+0.305));
const double A = pi*pow(rad,2);
const double I = pi*pow(rad,4)/4;
const double J = 2*I;
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

enum WORKER_STATE{
    CALCULATING_OBJFUNC,
    CALCULATING_JACOBIAN,
    RESTING,
    TERMINATED
};
static std::atomic<WORKER_STATE> worker_states[num_threads-1];

static void waitOnWorkers(){
    for(int i = 0; i < num_threads-1; i++)
        while( worker_states[i]!=RESTING ){ /* Execution is blocked */ }
}

static void setWorkerState(WORKER_STATE command){
    for(int i = 0; i < num_threads-1; i++) worker_states[i] = command;
}

static MatrixXd J_global = MatrixXd::Zero(36,36);
static VectorXd guess_global;
static Vector3d pL_shot[6];
static Matrix3d RL_shot[6];
static Vector3d nL[6];
static Vector3d mL[6];
const int N = 40;

static void integrateRod(int i){
    Vector3d n0 = Map<Vector3d>(&guess_global[5*i]);
    Vector2d m0xy = Map<Vector2d>(&guess_global[5*i+3]);
    double L = guess_global(30+i);

    VectorXd y0(18);
    y0 << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;

    //Numerically integrate the Cosserat rod equations
    MatrixXd Y = ode4<cosseratRodOde,N>(y0, L);

    pL_shot[i] = Y.block<3,1>(0, N-1);
    RL_shot[i] = Map<Matrix3d>(Y.block<9,1>(3, N-1).data());
    nL[i] = Y.block<3,1>(12, N-1);
    mL[i] = Y.block<3,1>(15, N-1);
}

//Shooting method objective function
static bool integrate = true;
static VectorXd shootingFunction(VectorXd guess){
    guess_global = guess;
    VectorXd residual(36);

    if(integrate){
        setWorkerState(CALCULATING_OBJFUNC);
        switch(num_threads){
            case 1: //Do all the work
                for(int i = 0; i < 6; i++) integrateRod(i);
                break;
            case 2: //Do half the work
                for(int i = 3; i < 6; i++) integrateRod(i);
                break;
            case 3: //Do a third of the work
                for(int i = 4; i < 6; i++) integrateRod(i);
                break;
            case 6: //Do a sixth of the work
                integrateRod(5);
                break;
            default:
                std::cout << "Invalid number of threads: must be divisor of 6" << std::endl;
                throw(1);
        }
        waitOnWorkers();
    }
    integrate = true;

    Vector3d EF = Vector3d::Zero();
    Vector3d EM = Vector3d::Zero();
    for(int i = 0; i < 6; i++){
        residual.segment<3>(5*i) = pL_shot[i] - (pE + RE*r[i]);
        residual.segment<2>(5*i+3) = linear_rotation_error(RL_shot[i],RE).segment<2>(0);
        EF -= nL[i];
        EM -= (mL[i] + (RE*r[i]).cross(nL[i]));
    }
    residual.segment<3>(30) = EF;
    residual.segment<3>(33) = EM;

    return residual;
}

const double incr = 1e-8;
static void setJacobianOfLeg(int i){
    for(int j = 0; j < 5; j++){
        double temp = guess_global(5*i+j);
        guess_global(5*i+j) += incr;
        Vector3d n0 = Map<Vector3d>(&guess_global[5*i]);
        Vector2d m0xy = Map<Vector2d>(&guess_global[5*i+3]);
        guess_global(5*i+j) = temp;
        double L = guess_global(30+i);
        VectorXd y0(18);
        y0 << p0[i], 1, 0, 0, 0, 1, 0, 0, 0, 1, n0, m0xy, 0;
        //Numerically integrate the Cosserat rod equations
        VectorXd yf = ode4_endpoint<cosseratRodOde,N>(y0, L);
        Vector3d pL_incr = yf.segment<3>(0);
        Matrix3d RL_incr = Map<Matrix3d>(&yf(3));
        Vector3d nL_incr = yf.segment<3>(12);
        Vector3d mL_incr = yf.segment<3>(15);

        J_global.block<3,1>(5*i, 5*i+j) = (pL_incr - pL_shot[i]) / incr;
        J_global.block<2,1>(5*i+3, 5*i+j) = linear_rotation_error( (RL_incr - RL_shot[i])/incr, RE).segment<2>(0);
        Vector3d nL_partial = (nL_incr - nL[i])/incr;
        J_global.block<3,1>(30, 5*i+j) = -nL_partial;
        J_global.block<3,1>(33, 5*i+j) = -( (mL_incr - mL[i]) / incr  + (RE*r[i]).cross(nL_partial) );
    }

    VectorXd y_s(18), yL(18);
    yL << pL_shot[i], Map<VectorXd>(RL_shot[i].data(), 9), nL[i], mL[i];
    cosseratRodOde(y_s, yL);

    J_global.block<3,1>(5*i, 30+i) = y_s.segment<3>(0);
    J_global.block<2,1>(5*i+3, 30+i) = linear_rotation_error( Map<Matrix3d>(&y_s[3]), RE).segment<2>(0);

    J_global.block<3,1>(30, 30+i) = -y_s.segment<3>(12);
    J_global.block<3,1>(33, 30+i) = -( y_s.segment<3>(15)  + (RE*r[i]).cross(y_s.segment<3>(12)));
}

static void jacobianFunction(MatrixXd& J, VectorXd&, VectorXd&){
    setWorkerState(CALCULATING_JACOBIAN);
    switch(num_threads){
        case 1:
            for(int i = 0; i < 6; i++) setJacobianOfLeg(i);
            break;
        case 2:
            for(int i = 3; i < 6; i++) setJacobianOfLeg(i);
            break;
        case 3:
            for(int i = 4; i < 6; i++) setJacobianOfLeg(i);
            break;
        case 6:
            setJacobianOfLeg(5);
    }
    waitOnWorkers();

    J = J_global;
}

static void workerFunction(int id){
    while( worker_states[id]!=TERMINATED ){
        if( worker_states[id]==CALCULATING_OBJFUNC ){
            switch(num_threads){
                case 2:
                    integrateRod(0);
                    integrateRod(1);
                    integrateRod(2);
                    break;
                case 3:
                    integrateRod(2*id);
                    integrateRod(2*id+1);
                    break;
                case 6:
                    integrateRod(id);
                    break;
            }
            worker_states[id] = RESTING;
        }

        if( worker_states[id]==CALCULATING_JACOBIAN ){
            switch(num_threads){
                case 2:
                    setJacobianOfLeg(0);
                    setJacobianOfLeg(1);
                    setJacobianOfLeg(2);
                    break;
                case 3:
                    setJacobianOfLeg(2*id);
                    setJacobianOfLeg(2*id+1);
                    break;
                case 6:
                    setJacobianOfLeg(id);
                    break;
            }
            worker_states[id] = RESTING;
        }
    }
}

int main(){
    std::thread workers[num_threads-1];
    setWorkerState(RESTING);
    for(int i = 0; i < num_threads-1; i++)
        workers[i] = std::thread(workerFunction, i);

    initialize_StewartGough_pattern();

    VectorXd guess(36);
    guess << VectorXd::Zero(30), pE(2)*VectorXd::Ones(6);

    //Solve with shooting method
    guess = solveLevenbergMarquardt<shootingFunction, jacobianFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    pE = Vector3d(0, 0.02, 0.48);
    guess = solveLevenbergMarquardt<shootingFunction, jacobianFunction>(guess, 1e-6, 500, 1e-4, 0.5);

    //Benchmark speed test
    std::clock_t start = std::clock();
    const int M = 5000;
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

    setWorkerState(TERMINATED);
    for(int i = 0; i < num_threads-1; i++) workers[i].join();

    return 0;
}
