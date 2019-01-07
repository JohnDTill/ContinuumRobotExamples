//Eigen is a popular and feature-rich C++ matrix math library
//Project Website: http://eigen.tuxfamily.org

#ifndef EIGENMATRIXTYPES_H
#define EIGENMATRIXTYPES_H

#include "eigen/Eigen/Dense"
using namespace Eigen;

//Extra types
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

#endif // EIGENMATRIXTYPES_H
