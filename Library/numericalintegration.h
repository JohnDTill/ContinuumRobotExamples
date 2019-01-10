#ifndef NUMERICALINTEGRATION_H
#define NUMERICALINTEGRATION_H

#include "eigenmatrixtypes.h"

namespace ContinuumRobotLibrary{

//A function pointer typedef for y_s = f(s,y)
typedef VectorXd(OdeFunc)(double,VectorXd);

/*! Integrate an ODE using the classic 4th-order Runge-Kutta algorithm. */
template<OdeFunc ODE, int N = 100>
inline MatrixXd ode4(VectorXd y0, double s0, double sf){
    MatrixXd Y( y0.size(), N );
    Y.col(0) = y0;

    double ds = (sf-s0)/(N-1);
    double half_ds = ds/2;
    double sixth_ds = ds/6;
    double L = sf-s0;
    int Nm1 = N-1;

    //Classic 4th-order Runge-Kutta method
    VectorXd k0, y1, k1, y2, k2, y3, k3;
    double s = s0;
    for(int i = 0; i < N-1; i++){
        y0 = Y.col(i);
        k0 = ODE(s, y0);
        s += half_ds;
        y1 = Y.col(i) + k0*half_ds;
        k1 = ODE(s, y1);
        y2 = Y.col(i) + k1*half_ds;
        k2 = ODE(s, y2);
        s = s0 + (L*(i+1))/Nm1;
        y3 = Y.col(i) + k2*ds;
        k3 = ODE(s, y3);

        Y.col(i+1) = Y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
    }

    return Y;
}

//A function pointer typedef for y_s = f(s,y), where y_s is an output argument
typedef void(OdeOutFunc)(VectorXd&, double s, VectorXd&);

/*! Integrate an ODE using the classic 4th-order Runge-Kutta algorithm. */
template<OdeOutFunc ODE, int N = 100>
inline MatrixXd ode4(VectorXd y0, double s0, double sf){
    EigenBase<VectorXd>::Index sze = y0.size();

    MatrixXd Y( sze, N );
    Y.col(0) = y0;

    double ds = (sf-s0)/(N-1);
    double half_ds = ds/2;
    double sixth_ds = ds/6;
    double L = sf-s0;
    int Nm1 = N-1;

    //Classic 4th-order Runge-Kutta method
    VectorXd k0(sze), k1(sze), k2(sze), k3(sze);
    double s = s0;
    for(int i = 0; i < Nm1; i++){
        ODE(k0,s,y0);
        y0 += k0*half_ds;
        s += half_ds;
        ODE(k1,s,y0);
        y0 = Y.col(i) + k1*half_ds;
        ODE(k2,s,y0);
        s = s0 + (L*(i+1))/Nm1;
        y0 = Y.col(i) + k2*ds;
        ODE(k3,s,y0);

        y0 = Y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
        Y.col(i+1) = y0;
    }

    return Y;
}

//A function pointer typedef for y_s = f(y), where there is no explicit dependence on s
typedef VectorXd(AutonomousOdeFunc)(VectorXd);

/*! Integrate an autonomous ODE (no explicit dependence on the variable of integration)
 *  using the classic 4th-order Runge-Kutta algorithm. */
template<AutonomousOdeFunc ODE, int N = 100>
inline MatrixXd ode4(VectorXd y0, double L){
    MatrixXd Y( y0.size(), N );
    Y.col(0) = y0;

    double ds = L/(N-1);
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    //Classic 4th-order Runge-Kutta method
    VectorXd k0, y1, k1, y2, k2, y3, k3;
    for(int i = 0; i < N-1; i++){
        y0 = Y.col(i);
        k0 = ODE(y0);
        y1 = Y.col(i) + k0*half_ds;
        k1 = ODE(y1);
        y2 = Y.col(i) + k1*half_ds;
        k2 = ODE(y2);
        y3 = Y.col(i) + k2*ds;
        k3 = ODE(y3);

        Y.col(i+1) = Y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
    }

    return Y;
}

//A function pointer typedef for y_s = f(y), where y_s is an output argument
typedef void(AutonomousOdeOutFunc)(VectorXd&, VectorXd&);

/*! Integrate an autonomous ODE (no explicit dependence on the variable of integration)
 *  using the classic 4th-order Runge-Kutta algorithm. */
template<AutonomousOdeOutFunc ODE, int N = 100>
inline MatrixXd ode4(VectorXd y0, double L){
    EigenBase<VectorXd>::Index sze = y0.size();

    MatrixXd Y( sze, N );
    Y.col(0) = y0;

    double ds = L/(N-1);
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    //Classic 4th-order Runge-Kutta method
    VectorXd k0(sze), k1(sze), k2(sze), k3(sze);
    for(int i = 0; i < N-1; i++){
        ODE(k0,y0);
        y0 += k0*half_ds;
        ODE(k1,y0);
        y0 = Y.col(i) + k1*half_ds;
        ODE(k2,y0);
        y0 = Y.col(i) + k2*ds;
        ODE(k3,y0);

        y0 = Y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
        Y.col(i+1) = y0;
    }

    return Y;
}

/*! Integrate an autonomous ODE (no explicit dependence on the variable of integration)
 *  using the classic 4th-order Runge-Kutta algorithm. */
template<AutonomousOdeOutFunc ODE, int N = 100>
inline VectorXd ode4_endpoint(VectorXd y0, double L){
    EigenBase<VectorXd>::Index sze = y0.size();

    double ds = L/(N-1);
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    //Classic 4th-order Runge-Kutta method
    VectorXd k0(sze), k1(sze), k2(sze), k3(sze), y(sze);
    for(int i = N-1; i > 0; i--){
        ODE(k0,y0);
        y = y0 + k0*half_ds;
        ODE(k1,y);
        y = y0 + k1*half_ds;
        ODE(k2,y);
        y = y0 + k2*ds;
        ODE(k3,y);

        y0 += (k0 + 2*(k1 + k2) + k3) * sixth_ds;
    }

    return y0;
}

}

#endif // NUMERICALINTEGRATION_H
