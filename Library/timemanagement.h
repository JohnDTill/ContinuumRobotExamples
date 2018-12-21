#ifndef TIMEMANAGEMENT_H
#define TIMEMANAGEMENT_H

#include "eigenmatrixtypes.h"

namespace ContinuumRobotLibrary{

/*! Manages a grid of time-differentiated variables.
    Z_h is set statically upon construction and updated
    dynamically with each call to "advanceTime()". */
class TimeManagerBdfAlpha{
private:
    MatrixXd* curr;
    MatrixXd old;
    MatrixXd* hist;

    double c2;
    double d1;
    double c1_plus_c0_times_d1;

public:
    TimeManagerBdfAlpha(MatrixXd& z, MatrixXd& zh, double dt, double alpha){
        curr = &z;
        hist = &zh;
        old = z;

        const double c0 = (1.5 + alpha) / (dt * ( 1 + alpha ));
        const double c1 = -2.0/dt;
        c2 = (0.5 + alpha) / (dt * ( 1 + alpha ));
        d1 = alpha / ( 1 + alpha );
        const double c1_plus_c2 = (-1.5 - alpha) / (dt * ( 1 + alpha ));
        c1_plus_c0_times_d1 = c1 + d1*c0;

        *hist = c1_plus_c2*z;
    }

    void advanceTime(){
        *hist = c1_plus_c0_times_d1*(*curr) + c2*old + d1*(*hist);
        old = *curr;
    }

    static double getC0(double dt, double alpha){
        return (1.5 + alpha) / (dt * ( 1 + alpha ));
    }
};

}

#endif // TIMEMANAGEMENT_H
