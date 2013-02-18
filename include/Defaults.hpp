#ifndef ___GFTOOLS_COMMON_H___
#define ___GFTOOLS_COMMON_H___

#include <complex>
#include <iostream>
#include<Eigen/Core>
#include<Eigen/StdVector>
//#include "EigenIterator.h"

#define REALTYPE_DOUBLE

namespace GFTools {

typedef double RealType;
typedef std::complex<RealType> ComplexType;

template <typename T>
using VectorType = Eigen::Matrix<T,Eigen::Dynamic, 1>;

template <typename T>
using MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor>;

/** A short name for imaginary unit. */
static const ComplexType I = ComplexType(0.0,1.0);    // 'static' to prevent linking problems

const RealType PI = std::atan(1.0)*4;
/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;

template <bool Fermion> inline ComplexType Matsubara(int n, RealType beta){return PI*I/beta*ComplexType(2*n+Fermion);};
template <bool Fermion> inline int MatsubaraIndex(ComplexType in, RealType beta){return std::round((beta*imag(in)/PI-Fermion)/2.0);};

inline ComplexType FMatsubara(int n, RealType beta){return Matsubara<1>(n,beta);};
inline ComplexType BMatsubara(int n, RealType beta){return Matsubara<0>(n,beta);};
inline int FMatsubaraIndex(ComplexType in, RealType beta){return MatsubaraIndex<1>(in,beta);};
inline int BMatsubaraIndex(ComplexType in, RealType beta){return MatsubaraIndex<0>(in,beta);};

} // end namespace GFTools

#endif // endif::ifndef ___GFTOOLS_COMMON_H___
