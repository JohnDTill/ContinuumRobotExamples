#ifndef NUMERICALPDE_H
#define NUMERICALPDE_H

#include "eigenmatrixtypes.h"

namespace ContinuumRobotLibrary{

//A function pointer typedef for [y_s,z] = f(y)
typedef void(AutonomousOdeZFunc)(VectorXd&,Map<VectorXd>,Map<VectorXd>); //[y_s,z] = f(y)

/*! Integrate an ODE using Euler's method, while setting a z vector to solve initial conditions. */
template<AutonomousOdeZFunc ODE, int statesize, int derivsize, int N = 100>
inline void euler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        ODE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    ODE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize));
    #undef y_s
}

//A function pointer typedef for [y_s,z] = f(y,z_h)
typedef void(AutonomousPdeFunc)(VectorXd&,Map<VectorXd>,Map<VectorXd>,Map<VectorXd>);

/*! Integrate a semi-discretized PDE using Euler's method.
    The PDE function must set time-differientiated variables in the z vector,
    and the history terms z_h are calculated elsewhere. */
template<AutonomousPdeFunc PDE, int statesize, int derivsize, int N = 100>
inline void TimeBdfAlpha_SpaceEuler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L,
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        PDE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    PDE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize), Map<VectorXd>(&Z_h(0,N-1),derivsize));
    #undef y_s
}

}

#endif // NUMERICALPDE_H
