#pragma once

#define BOOST_RESULT_OF_USE_DECLTYPE

#include<complex>
#include<iostream>
#include<utility>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-register" 
    #include<Eigen/Core>
    #include<Eigen/StdVector>
#pragma GCC diagnostic pop



namespace gftools {

#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#ifndef NDEBUG
#define DEBUG(MSG)            std::cout << MSG_PREFIX << MSG << std::endl;
#else
#define DEBUG(MSG)
#endif
#define INFO(MSG)             std::cout << MSG << std::endl;
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush;
#define INFO2(MSG)            std::cout << "\t"     << MSG << std::endl;
#define INFO3(MSG)            std::cout << "\t\t"   << MSG << std::endl;
#define INFO4(MSG)            std::cout << "\t\t\t" << MSG << std::endl;
#define ERROR(MSG)            std::cerr << MSG_PREFIX << MSG << std::endl;

typedef double real_type;
typedef std::complex<real_type> complex_type;

/** A short name for imaginary unit. */
static const complex_type I = complex_type(0.0,1.0);    // 'static' to prevent linking problems
/** A short name for pi. */
const real_type PI = std::atan(1.0)*4;

/*
template <typename T>
using vector_type = Eigen::Matrix<T,Eigen::Dynamic, 1>;

template <typename T>
using matrix_type = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor>;
*/
/** Dense complex matrix. */
//typedef Eigen::Matrix<complex_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
//typedef Eigen::Matrix<real_type,Eigen::Dynamic,Eigen::Dynamic,Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

/** Dense complex vector. */
//typedef Eigen::Matrix<complex_type,Eigen::Dynamic,1,Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
//typedef Eigen::Matrix<real_type,Eigen::Dynamic,1,Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
//typedef Eigen::Matrix<int,Eigen::Dynamic,1,Eigen::AutoAlign> IntVectorType;

} // end namespace gftools

