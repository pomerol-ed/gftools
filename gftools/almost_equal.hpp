#pragma once

#include <limits>
//C++11 - use boost::type_traits for C++03
#include <type_traits> 
// consider removing it - it is only used for complex_number
#include "defaults.hpp"

namespace gftools {

/// A helper structure to perform type-dependent "almost equal" comparison
/// The usage of this structure allows to make comparison for user-defined types, in which
/// case the template specialization of almost_equal_checker should be provided
template <typename, typename > struct almost_equal_checker;

/// "almost equal" comparison - performs tolerance-based comparison for provided objects
///TODO: note duplication of functionality with the more fancy is_float_equal function in tools.
template <typename T1,typename T2>
    bool almost_equal(T1 t1, T2 t2, real_type tol = 10.*std::numeric_limits<real_type>::epsilon()){ 
    return almost_equal_checker<T1,T2>::is_equal(t1,t2,tol); 
}

/// A helper structure that defines a "true_type" compile-time constant if a given type can be used as a number
/// In this realization everything that can be casted to a complex_number is a number.
template <typename T> struct is_number;

// implementations
template <typename T1, typename T2 = T1>
struct almost_equal_checker {
    static typename std::enable_if<is_number<T1>::value && is_number<T2>::value, bool>::type is_equal(T1 t1, T2 t2, real_type tol) { 
        return bool(std::abs(complex_type(t1) - complex_type(t2))<tol); };
};

template <typename T>
struct is_number : std::integral_constant<
    bool,
    std::is_arithmetic<typename std::remove_cv<T>::type>::value ||
    std::is_convertible<typename std::remove_cv<T>::type,complex_type>::value ||
    std::is_convertible<typename std::remove_cv<T>::type,double>::value
    >{};




} // end namespace gftools

