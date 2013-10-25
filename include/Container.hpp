#pragma once
#ifndef ___GFTOOLS_CONTAINER_H___
#define ___GFTOOLS_CONTAINER_H___

//#define BOOST_RESULT_OF_USE_DECLTYPE
#include "Defaults.hpp"
#include "Tools.hpp"
#include <type_traits>

#include <boost/multi_array.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace GFTools { 

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase;

template <typename ValueType, size_t N>
struct Container;

template <typename ValueType, size_t N>
struct _container_underlying_type;

// ================================================================== //

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase 
{
    constexpr static size_t _N = N;

    typedef ValueType value_t;
    typedef BoostContainerType boost_t;
    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    typedef Eigen::Map<EigenArray> EigenMap;

    mutable BoostContainerType _data;
    typedef typename _container_underlying_type<ValueType,N>::type under_type;
    typedef decltype(_data[0]) boost_under_type;
    typedef std::function<under_type(boost_under_type)> action_type;

    ContainerBase(const BoostContainerType &in):_data(in){};
//    template <typename ...Args> ContainerBase(Args... in):_data(BoostContainerType(in...)){};
    /** Copy constructor. */
    ContainerBase(const ContainerBase<ValueType,N,BoostContainerType> &rhs):_data(rhs._data){};
    ContainerBase& operator=(const ContainerBase<ValueType,N,BoostContainerType> &rhs){_data = rhs._data; return (*this);};
    template <class B2> ContainerBase& operator=(const ContainerBase<ValueType,N,B2> &rhs){_data = rhs._data; return (*this);};
    /** Move constructor. */
    ContainerBase(ContainerBase<ValueType,N,BoostContainerType> &&rhs):_data(rhs._data){};
    ContainerBase& operator=(ContainerBase<ValueType,N,BoostContainerType> &&rhs){std::swap(_data,rhs._data); return (*this);};

    /** Access operators. */
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i) { return _data[i]; };
    template<typename U = typename std::enable_if<N==1>::type> const ValueType& operator[](size_t i) const { return _data[i]; };

    //template<typename U = typename std::enable_if<N!=1>::type> 
    //    ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> operator[](size_t i);
    //template<typename U = typename std::enable_if<N!=1>::type> 
    //    const ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> operator[](size_t i) const;

    template<typename U = typename std::enable_if<N!=1>::type> 
        auto operator[](size_t i) -> ContainerBase<ValueType, N-1, decltype(_data[0])> { return  ContainerBase<ValueType, N-1, decltype(_data[0])>(_data[i]); }
    template<typename U = typename std::enable_if<N!=1>::type> 
        auto operator[](size_t i) const -> ContainerBase<ValueType, N-1, decltype(_data[0])> const { return  ContainerBase<ValueType, N-1, decltype(_data[0])>(_data[i]); }

    /** Return operators. */
    template<typename U = typename std::enable_if<N==2, bool>> MatrixType<ValueType> getAsMatrix() const;
    template<typename U = typename std::enable_if<N==2, bool>> ContainerBase<ValueType,N,BoostContainerType>& operator=(MatrixType<ValueType> rhs);
    /** Return a diagonal matrix, corresponding to the object */
    template<typename U = typename std::enable_if<N==2, bool>> MatrixType<ValueType> getAsDiagonalMatrix() const;
    /** Return a vector, corresponding to the object. */
    template<typename U = typename std::enable_if<N==2, bool>> VectorType<ValueType> getAsVector() const;
 
    /** Conjugate. */
    template <typename U = typename std::enable_if<std::is_same<ValueType, ComplexType>::value, int>::type> 
        ContainerBase conj();
    /** Recursively iterates and sums all values in the container. */
    ValueType sum(){return EigenMap(_data.origin(), _data.num_elements()).sum();};

    /** Mathematical operators. */

    template <typename T>
    using isContainer = typename std::enable_if<T::_N>=1, bool>::type;
    template <typename T>
    using isValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, bool>::type;

    ContainerBase& operator=(ValueType rhs);
    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const R &rhs);//ContainerBase<ValueType,N,R> &rhs); 
    template <typename R, isContainer<R> = 0> 
        Container<ValueType,N> operator+(const R &rhs) const;
    template <typename R2, isValue<R2> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        Container<ValueType,N> operator+(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        Container<ValueType,N> operator-(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        Container<ValueType,N> operator-(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        Container<ValueType,N> operator*(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        Container<ValueType,N> operator*(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        Container<ValueType,N> operator/(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        Container<ValueType,N> operator/(const R2& rhs) const;

    friend Container<ValueType,N> operator* (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs*lhs;};
    friend Container<ValueType,N> operator+ (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs+lhs;};
    friend Container<ValueType,N> operator- (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs*(-1.0)+lhs;};
    friend Container<ValueType,N> operator/ (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {Container<ValueType,N> out(rhs); out=lhs; return out/rhs;};
    
    /** Make the object streamable. */
    template <typename V1, size_t M, typename B>
    friend std::ostream& operator<<(std::ostream& lhs, const ContainerBase<V1, M, B> &in);

    typedef boost::transform_iterator<action_type,typename BoostContainerType::iterator> iterator; 
    typedef boost::transform_iterator<action_type,typename BoostContainerType::iterator> const_iterator; 
    /** Begin iterator. */
    template<typename U = typename std::enable_if<N!=1>::type>
    const_iterator begin() const { 
        //DEBUG("Using N!=1");
        auto f1 = [this](boost_under_type in){return ContainerBase<ValueType, N-1, boost_under_type>(in);}; //(boost_under_type in){return _act(in);}; 
        return boost::make_transform_iterator(_data.begin(),f1);
        };
    template<typename U = typename std::enable_if<N==1>::type>
    const iterator begin() const {
        //DEBUG("Using N=1");
        auto f1 = [this](boost_under_type in){return in;} ;
        return boost::make_transform_iterator(_data.begin(),f1);
    }

    /** End iterator. */
    template<typename U = typename std::enable_if<N!=1>::type>
    const_iterator end() const { 
        //action_type f1 = [&](boost_under_type in){return _act(in);}; 
        auto f1 = [this](boost_under_type in){return ContainerBase<ValueType, N-1, boost_under_type>(in);}; //(boost_under_type in){return _act(in);}; 
        //auto f1 = [&](boost_under_type in){return _act(in);}; 
        return boost::make_transform_iterator(_data.end(),f1);
    };
    template<typename U = typename std::enable_if<N==1>::type>
    const iterator end() const {
        //DEBUG("Using N=1");
        auto f1 = [this](boost_under_type in){return in;} ;
        return boost::make_transform_iterator(_data.end(),f1);
    };

    class exWrongIndex : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
};

template <typename ValueType, size_t N>
struct Container : ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>> {
    typedef boost::multi_array<ValueType, N> BoostContainerType;
    typedef ContainerBase<ValueType,N,BoostContainerType> Base;
    using Base::_data;
    Container(std::array<size_t,N> shape):ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    template <typename CT>
    Container(ContainerBase<ValueType,N,CT> in) : ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>>(in._data) {};
    template<typename U = typename std::enable_if<N==2, bool>> Container<ValueType,N> (MatrixType<ValueType> rhs);
    using Base::operator+=;
    using Base::operator-=;
    using Base::operator*=;
    using Base::operator/=;
    using Base::operator=;
};

template<typename ValueType, size_t N, size_t M = N>
using ContainerView = ContainerBase<RealType, N, typename Container<RealType, N>::boost_t::template array_view<M>::type>;

template <typename ValueType, size_t N>
struct _container_underlying_type
{
    typedef ContainerBase<ValueType,N-1,boost::detail::multi_array::sub_array<ValueType, N-1>> type;
};

template <typename ValueType>
struct _container_underlying_type<ValueType,1>
{
    typedef ValueType type;
};

}; // end of namespace GFTools
#endif

#include "Container.hxx"
