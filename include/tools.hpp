#pragma once

#include "defaults.hpp"
#include "tuple_tools.hpp"

namespace gftools {
namespace tools {

//
// type traits
//

/** Everything that can be casted to a complex_number is a number. */
template <typename T>
struct is_number : std::integral_constant<
    bool,
    std::is_arithmetic<typename std::remove_cv<T>::type>::value ||
    std::is_convertible<typename std::remove_cv<T>::type,complex_type>::value ||
    std::is_convertible<typename std::remove_cv<T>::type,double>::value
    >{};

//
// automatic argument generators
//

/** A tool to generate a type for an object T of D Args. */
template <size_t N, typename ArgType1, template <typename ...> class T, typename ...ArgTypes>  
struct ArgBackGenerator:ArgBackGenerator<N-1,ArgType1,T,ArgTypes...,ArgType1 > {};

template <typename ArgType1, template <typename ...> class T, typename ...ArgTypes> 
struct ArgBackGenerator<1, ArgType1, T, ArgTypes...> { typedef T<ArgTypes..., ArgType1> type; };

/** A tool to generate a type for a function of N Args. */
template <size_t N, typename ValueType, typename ArgType1, typename ...ArgTypes>  
struct ArgFunGenerator : ArgFunGenerator<N-1,ValueType,ArgType1, ArgTypes...,ArgType1 >{static_assert(N>1,"");};

template <typename ValueType, typename ArgType1, typename ...ArgTypes> 
struct ArgFunGenerator<1, ValueType, ArgType1, ArgTypes...> { 
    typedef std::function<ValueType(ArgTypes..., ArgType1)> type; };

//
// equality checker 
//

/** Set of tools to check float equality of numbers or tuples. */
template <typename > struct equality_checker;

template <typename T>
struct equality_checker {
    typedef T type;
    static typename std::enable_if<is_number<T>::value, bool>::type is_equal(T t1, T t2, real_type tol) { 
        return bool(std::abs(complex_type(t1) - complex_type(t2))<tol); };
};

// free function
/** Returns true, if two objects are equal up to a certain tolerance. */
template <typename T, typename = typename std::enable_if<is_number<T>::value, bool>::type> 
    bool is_float_equal(T t1, T t2, real_type tol = 10.*std::numeric_limits<real_type>::epsilon()){ return equality_checker<T>::is_equal(t1,t2,tol); };

template <typename T1, typename T2, typename = typename std::enable_if<is_number<T1>::value, bool>::type, typename = typename std::enable_if<is_number<T2>::value, bool>::type> 
    bool is_float_equal(T1 t1, T2 t2, real_type tol = 10.*std::numeric_limits<real_type>::epsilon()){ return equality_checker<T1>::is_equal(t1,t2,tol); };
} // end of namespace tools 

namespace tuple_tools { 
template <int ... S, typename ... Args>
bool is_float_equal_(extra::arg_seq<S...>, const  std::tuple<Args...>& l, const std::tuple<Args...>& r, double tol)
{
    std::array<bool,sizeof...(Args)> res = {{ (tools::is_float_equal(std::get<S>(l), std::get<S>(r), tol))... }};
    return std::all_of(res.begin(), res.end(), [](bool x){return x;});
}
}

namespace tools { 
template <typename ... Args>
bool is_float_equal (const std::tuple<Args...>& l, const std::tuple<Args...>& r, double tol)
{
    return tuple_tools::is_float_equal_(typename tuple_tools::extra::index_gen<sizeof...(Args)>::type(), l, r, tol);
}


//
/** A tool to calc an integer power function of an int. */
/*
template<int base, unsigned exponent >
struct __power {
    enum { value = base * __power<base, exponent - 1>::value };
};
template< int base >
struct __power<base,0> {
    enum { value = 1 };
};
*/

/** Function traits. */
template <typename FunctionType> struct fun_traits;
template <typename ValType, typename ... ArgTypes> 
struct fun_traits<std::function<ValType(ArgTypes...)> >
{
    static std::function<ValType(ArgTypes...)> constant(ValType c) 
    { return [c](ArgTypes...in){return c;};}
    static std::function<ValType(ArgTypes...)> getFromTupleF(std::function<ValType(std::tuple<ArgTypes...>)> f1)
    { return [f1](ArgTypes...in){return f1(std::forward_as_tuple(in...));};}
    static std::function<ValType(std::tuple<ArgTypes...>)> getTupleF(std::function<ValType(ArgTypes...)> f1)
    { return [f1](const std::tuple<ArgTypes...> &in)->ValType{ tuple_tools::unfold_tuple(f1,in);};}
};



} // end of namespace tools
} // end of namespace gftools
