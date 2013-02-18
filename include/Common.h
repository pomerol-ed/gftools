#ifndef ___FK_FK_H___
#define ___FK_FK_H___

#include<iostream>
#include<complex>
#include<string>
#include<vector>
#include<list>
#include<map>
#include<array>
#include<iomanip>
#include<fstream>

#include<memory>
#include<utility>
#include<functional>

#include<Eigen/Core>
#include<Eigen/StdVector>

#include "EigenIterator.h"

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

/** A wrapper around numbers for file output. */
template <typename T> struct __num_format;
template <typename T> std::ostream& operator<<(std::ostream& lhs, const __num_format<T> &in);
template <typename T> std::istream& operator>>(std::istream& lhs, __num_format<T> &out);
template <typename T> 
struct __num_format {
    static const int _prec = 12;
    T _v;
    __num_format(T v):_v(v){};
    operator T(){return _v;};
    void savetxt(const std::string& filename) { 
        std::cout << "Saving " << typeid(*this).name() << " to " << filename << std::endl;
        std::ofstream out; out.open(filename.c_str()); out << *this << std::endl; out.close(); 
    }; 
    friend std::ostream& operator<< <>(std::ostream& lhs, const __num_format<T> &in);
    friend std::istream& operator>> <>(std::istream& lhs, __num_format<T> &out);
};

template <typename T> 
inline std::ostream& operator<<(std::ostream& lhs, const __num_format<T> &in) {lhs << std::setprecision(in._prec) << in._v; return lhs;};
template <typename T> 
inline std::istream& operator>>(std::istream& lhs, __num_format<T> &out) {lhs >> out._v; return lhs;};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format<ComplexType> &in){lhs << std::setprecision(in._prec) << real(in._v) << " " << imag(in._v); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, __num_format<ComplexType> &out){RealType re,im; lhs >> re; lhs >> im; out._v = re+I*im; return lhs;};

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

/** A tool to calc an integer power function of an int. */
template<int base, unsigned exponent >
struct __power {
    enum { value = base * __power<base, exponent - 1>::value };
};
template< int base >
struct __power<base,0> {
    enum { value = 1 };
};

/** A tool to wrap a call a class method from a tuple. */
template<int ...> struct __seq {};
template<int N, int ...S> struct __gens : __gens<N-1, N-1, S...> {};
template<int ...S> struct __gens<0, S...>{ typedef __seq<S...> type; };

template <typename ReturnType, typename ...Args> struct __caller { 
    std::tuple<Args...> _params;
    std::function<ReturnType(Args...)> _f;
    template<int ...S> ReturnType _callf(__seq<S...>) { return _f(std::get<S>(_params)...); };
    ReturnType call(){ return _callf(typename __gens<sizeof...(Args)>::type()); };
};

template <typename Arg, size_t D, typename... Extras> struct __repeater : __repeater<Arg,D-1,Arg,Extras...>  {  
    //template <int ...S> std::array<Arg,D> _get(const std::tuple<Arg...>& in){return { std::get<S>
    template <int ...S, typename ...ArgTypes> static std::array<Arg,D> __get_arr(__seq<S...>, std::tuple<ArgTypes...> in) { return  {{ std::get<S>(in)... }}; }
    static std::array<Arg,D> get_array(const Arg& in){ return __get_arr(typename __gens<D>::type(), get_tuple(in)); };
    static typename __repeater::TupleType get_tuple(const Arg& in){ return std::tuple_cat(std::forward_as_tuple(in),__repeater<Arg,D-1>::get_tuple(in)); }; 
};

template <typename Arg, typename ... Extras> struct __repeater<Arg,1,Extras...> { 
    typedef std::tuple<Arg,Extras...> TupleType;
    static std::array<Arg,1> get_array(const Arg& in){ return {{ in }}; };
    static TupleType get_tuple(const Arg& in){ return std::forward_as_tuple(in); }
};


/** Function traits. */
template <typename FunctionType> struct __fun_traits;
template <typename ValType, typename ... ArgTypes> 
struct __fun_traits<std::function<ValType(ArgTypes...)> >
{
    static std::function<ValType(ArgTypes...)> constant(ValType c) 
    { return [c](ArgTypes...in){return c;};}
    /*static std::function<ValType(ArgTypes...)> add(const std::function<ValType(ArgTypes...)> &f1, const std::function<ValType(ArgTypes...)>& f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)+f2(in...);}; }
    static std::function<ValType(ArgTypes...)> multiply(const std::function<ValType(ArgTypes...)> &f1, const std::function<ValType(ArgTypes...)> &f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)*f2(in...);}; }
    static std::function<ValType(ArgTypes...)> subtract(const std::function<ValType(ArgTypes...)> &f1, const std::function<ValType(ArgTypes...)> &f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)-f2(in...);}; }
    static std::function<ValType(ArgTypes...)> divide(const std::function<ValType(ArgTypes...)> &f1, const std::function<ValType(ArgTypes...)> &f2)
    { return [f1,f2](ArgTypes... in){return f1(in...)*f2(in...);}; }
    */
    static std::function<ValType(ArgTypes...)> getFromTupleF(std::function<ValType(std::tuple<ArgTypes...>)> f1)
    { return [f1](ArgTypes...in){return f1(std::forward_as_tuple(in...));};}
    static std::function<ValType(std::tuple<ArgTypes...>)> getTupleF(std::function<ValType(ArgTypes...)> f1)
    { return [f1](const std::tuple<ArgTypes...> &in)->ValType{ __caller<ValType,ArgTypes...> t; t = {in, f1}; return t.call();};}
};

/** A tool to split a tuple from http://stackoverflow.com/questions/10626856/how-to-split-a-tuple. */
template <typename Target, typename Tuple, int N, bool end >
struct __split_tuple_struct
{
    template < typename ... Args >
    static Target create(Tuple const& t, Args && ... args)
    {
        return __split_tuple_struct<Target,Tuple, N+1, std::tuple_size<Tuple>::value == N+1>::create(t, std::forward<Args>(args)..., std::get<N>(t));
    }
};

template < typename Target, typename Tuple, int N >
struct __split_tuple_struct<Target,Tuple,N,true>
{
    template < typename ... Args >
    static Target create(Tuple const& t, Args && ... args) { return Target(std::forward<Args>(args)...); }
};

template < typename Head, typename ... Tail >
std::tuple<Tail...> __tuple_tail(std::tuple<Head,Tail...> const& tpl)
{
    return __split_tuple_struct<std::tuple<Tail...>, std::tuple<Head,Tail...>, 1, std::tuple_size<std::tuple<Head,Tail...>>::value == 1>::create(tpl);
}

/** Print a tuple. */
template <class T> struct __tuple_print;
template <typename...ArgTypes> struct __tuple_print<std::tuple<ArgTypes...>> 
    { static void print(std::tuple<ArgTypes...> in ){std::cout << std::get<0>(in) << " " << std::flush; auto a = __tuple_tail(in); __tuple_print<decltype(a)>::print(a);};
    };
template <typename ArgType> struct __tuple_print<std::tuple<ArgType>> 
    { static void print(std::tuple<ArgType> in ){std::cout << std::get<0>(in) << std::endl;};
    };

} // end namespace GFTools

#endif // endif::ifndef ___FK_FK_H___
