#ifndef CONVEXOPTIMIZATION_H
#define CONVEXOPTIMIZATION_H

#include <iostream>
#include "eigenmatrixtypes.h"
#include "numericaldifferentiation.h"

namespace ContinuumRobotLibrary{

//Error handling helper functions
const int SOLVER_FAILURE = 15;

static void checkInitialResidualForNaN(VectorXd& r){
    if( r.hasNaN() ){
        std::cout << "ERROR: The objective function is undefined at the initial guess:\n"
                     "f = " << r.transpose() << '\n' << std::endl;
        throw(SOLVER_FAILURE);
    }
}

static void checkJacobianForNaN(MatrixXd& J){
    if( J.hasNaN() ){
        std::cout << "ERROR: The Jacobian has large elements beyond machine precision.\n";
        if(J.rows()*J.cols() < 5000) std::cout << "J =\n" << J << "\n\n";
        std::cout.flush();
        throw(SOLVER_FAILURE);
    }
}

static void incrementIterationCounter(int& counter, VectorXd& y, VectorXd& r, MatrixXd& J, int& max_iter){
    if( ++counter > max_iter ){
        std::cout << "ERROR: Max iter of " << max_iter << " exceeded.\n"
                     "Final sum of squares: E = " << r.dot(r) << "\n"
                     "Final residual: f = " << r.transpose() << "\n"
                     "Final stalled guess: y = " << y.transpose() << "\n\n";
        FullPivLU<MatrixXd> lu_decomp(J);
        if( lu_decomp.rank() < r.size() && y.size() >= r.size() ){
            std::cout << "The jacobian is rank deficient.\n";
            if(J.rows()*J.cols() < 5000) std::cout << "J =\n" << J << "\n\n";
        }
        std::cout.flush();
        throw(SOLVER_FAILURE);
    }
}

/*! Solve a convex system with a Levenberg algorithm.
    The last call to the objective function f is guaranteed to to use the final solved y. */
template<ObjFunc f>
inline VectorXd solveLevenbergMarquardt(VectorXd y0,
                        double sos_tol = 1e-12,
                        int max_iter = 500,
                        double damp = 1e-2,
                        double alm_adptv_coeff = 0.5,
                        double incr_scale = 1e-8,
                        double incr_floor = 1e-12){
    VectorXd r = f(y0);
    checkInitialResidualForNaN(r);
    double E = r.squaredNorm();
    if( E <= sos_tol ) return y0;

    VectorXd y  = y0;
    VectorXd ones = VectorXd::Ones(y.size());
    MatrixXd J = MatrixXd::Zero(r.size(), y.size());
    MatrixXd B;
    int counter = 0;

    while( E > sos_tol ){
        getJacobian<f>(J,y,r,incr_scale,incr_floor);
        checkJacobianForNaN(J);
        incrementIterationCounter(counter, y, r, J, max_iter);
        B.noalias() = J.transpose()*J;
        VectorXd rhs = -J.transpose()*r;

        damp *= alm_adptv_coeff;
        //lhs is symmetric, positive semi-definite, positive definite for damp > 0
        B.diagonal() += damp*ones;
        VectorXd pc = B.selfadjointView<Eigen::Lower>().llt().solve(rhs);
        r = f(y + pc);
        double Ec = r.squaredNorm();

        if( !r.hasNaN() && Ec < E ){
            y += pc;
            E = Ec;
        }else{
            while( r.hasNaN() || Ec > E ){
                incrementIterationCounter(counter, y, r, J, max_iter);

                double old_damp = damp;
                damp /= alm_adptv_coeff;

                B.diagonal() += (damp - old_damp)*ones;
                pc.noalias() = B.selfadjointView<Eigen::Lower>().llt().solve(rhs);
                r.noalias() = f(y + pc);
                Ec = r.squaredNorm();
            }

            y += pc;
            E = Ec;
        }
    }

    return y;
}

//Function pointer to calculate J = f(y,r) using an output argument for J.
//J is the Jacobian, y the current guess, and r the residual.
typedef void(JacobFunc)(MatrixXd&,VectorXd&,VectorXd&);

/*! Solve a convex system with a Levenberg algorithm with a user-supplied Jacobian.
    The last call to the objective function f is guaranteed to to use the final solved y.
    The Jacobian is initialized to zeroes, so only non-zero elements need to be assigned. */
template<ObjFunc f, JacobFunc suppliedJacobianFunction>
inline VectorXd solveLevenbergMarquardt(VectorXd y0,
                        double sos_tol = 1e-12,
                        int max_iter = 500,
                        double damp = 1e-2,
                        double alm_adptv_coeff = 0.5){
    VectorXd r = f(y0);
    checkInitialResidualForNaN(r);
    double E = r.squaredNorm();
    if( E <= sos_tol ) return y0;

    VectorXd y  = y0;
    VectorXd ones = VectorXd::Ones(y.size());
    MatrixXd J = MatrixXd::Zero(r.size(),y.size());
    MatrixXd B;
    int counter = 0;

    while( E > sos_tol ){
        suppliedJacobianFunction(J,y,r);
        checkJacobianForNaN(J);
        incrementIterationCounter(counter, y, r, J, max_iter);
        B.noalias() = J.transpose()*J;
        VectorXd rhs = -J.transpose()*r;

        damp *= alm_adptv_coeff;
        //lhs is symmetric, positive semi-definite, positive definite for damp > 0
        B.diagonal() += damp*ones;
        VectorXd pc = B.selfadjointView<Eigen::Lower>().llt().solve(rhs);
        r = f(y + pc);
        double Ec = r.squaredNorm();

        if( !r.hasNaN() && Ec < E ){
            y += pc;
            E = Ec;
        }else{
            while( r.hasNaN() || Ec > E ){
                incrementIterationCounter(counter, y, r, J, max_iter);

                double old_damp = damp;
                damp /= alm_adptv_coeff;

                B.diagonal() += (damp - old_damp)*ones;
                pc.noalias() = B.selfadjointView<Eigen::Lower>().llt().solve(rhs);
                r.noalias() = f(y + pc);
                Ec = r.squaredNorm();
            }

            y += pc;
            E = Ec;
        }
    }

    return y;
}

}

#endif // CONVEXOPTIMIATION_H
