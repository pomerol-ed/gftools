#pragma once
#ifndef ___GFTOOLS_CONTAINER_H___
#define ___GFTOOLS_CONTAINER_H___

#include "Defaults.hpp"
#include "Tools.hpp"
#include <type_traits>

#include <boost/multi_array.hpp>
#include <boost/operators.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace GFTools { 

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase;

template <typename ValueType, size_t N>
struct _container_underlying_type;

// ================================================================== //

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase : boost::operators<ContainerBase<ValueType,N,BoostContainerType>>
{

    typedef ValueType VType;
//    template <size_t M> using boost_view_type = typename boost::multi_array<ValueType, N>::template array_view<M>::type;
//    typedef boost::multi_array<ValueType, N> boost_array_type;
//    typedef boost::multi_array_ref<ValueType, N> boost_array_ref_type;

    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    typedef Eigen::Map<EigenArray> EigenMap;

    mutable BoostContainerType _data;

    ContainerBase(const BoostContainerType &in):_data(in){};
//    template <typename ...Args> ContainerBase(Args... in):_data(BoostContainerType(in...)){};
    /** Copy constructor. */
    ContainerBase(const ContainerBase<ValueType,N,BoostContainerType> &rhs):_data(rhs._data){};
    ContainerBase& operator=(const ContainerBase<ValueType,N,BoostContainerType> &rhs){_data = rhs._data; return (*this);};
    /** Move constructor. */
    ContainerBase(ContainerBase<ValueType,N,BoostContainerType> &&rhs):_data(rhs._data){};
    ContainerBase& operator=(ContainerBase<ValueType,N,BoostContainerType> &&rhs){std::swap(_data,rhs._data); return (*this);};

    /** Access operators. */
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i) { return _data[i]; };
    template<typename U = typename std::enable_if<N==1>::type> const ValueType& operator[](size_t i) const { return _data[i]; };

    template<typename U = typename std::enable_if<N!=1>::type> 
        ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> operator[](size_t i);
    template<typename U = typename std::enable_if<N!=1>::type> 
        const ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> operator[](size_t i) const;

    //template<typename U = typename std::enable_if<N!=1>::type> 
    //    auto operator[](size_t i) -> ContainerBase<ValueType, N-1, decltype(_data[0])> { return  ContainerBase<ValueType, N-1, decltype(_data[0])>(_data[i]); }
    //template<typename U = typename std::enable_if<N!=1>::type> 
    //    auto operator[](size_t i) const -> ContainerBase<ValueType, N-1, decltype(_data[0])> const { return  ContainerBase<ValueType, N-1, decltype(_data[0])>(_data[i]); }

    /** Return operators. */
    template<typename U = typename std::enable_if<N==2, bool>> MatrixType<ValueType> getAsMatrix() const;
    template<typename U = typename std::enable_if<N==2, bool>> ContainerBase<ValueType,N,BoostContainerType>& operator=(MatrixType<ValueType> &&rhs);
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
    using isContainer = typename std::enable_if<std::is_convertible<typename T::VType,typename ContainerBase::VType>::value, bool>::type;
    template <typename T>
    using isValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, bool>::type;

    ContainerBase& operator=(ValueType rhs){}
    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const R &rhs);//ContainerBase<ValueType,N,R> &rhs); 
    template <typename R, isContainer<R> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator+(const R &rhs) const;
    template <typename R2, isValue<R2> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator+(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator-(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator-(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator*(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator*(const R2& rhs) const;

    template <typename R, isContainer<R> = 0>
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const R &rhs); 
    template <typename R, isContainer<R> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator/(const R &rhs) const;
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const R2& rhs);
    template <typename R2, isValue<R2> = 0> 
        ContainerBase<ValueType,N,BoostContainerType> operator/(const R2& rhs) const;
    
    //ContainerBase<ValueType,N,BoostContainerType>& operator=(const ValueType &rhs);
    /** Make the object streamable. */
    template <typename V1, size_t M, typename B>
    friend std::ostream& operator<<(std::ostream& lhs, const ContainerBase<V1, M, B> &in);

    //typedef boost::transform_iterator<double,_container_underlying_type<ValueType,N>> iterator_t; 
//    struct iterator;
//    struct const_iterator;
    /** Begin iterator. */
    //const_iterator begin() const;
    /** End iterator. */
    //const_iterator end() const;

    class exWrongIndex : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
};

//template <typename ValueType, size_t N, class BoostContainerType>
//struct ContainerBase<ValueType,N,BoostContainerType>::iterator :
//std::iterator<std::bidirectional_iterator_tag, _container_underlying_type<ValueType,N>>
//{
//};

template <typename ValueType, size_t N>
//using Container = ContainerBase<ValueType, N, typename boost::multi_array<ValueType, N>>;
struct Container : ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>> {
    typedef boost::multi_array<ValueType, N> BoostContainerType;
    typedef ContainerBase<ValueType,N,BoostContainerType> Base;
    using Base::_data;
    Container(std::array<size_t,N> shape):ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    template <typename CT>
    Container(ContainerBase<ValueType,N,CT> in) : ContainerBase<ValueType,N,typename boost::multi_array<ValueType, N>>(in._data) {};
    template<typename U = typename std::enable_if<N==2, bool>> Container<ValueType,N> (MatrixType<ValueType> &&rhs);
    using Base::operator+=;
    using Base::operator-=;
    using Base::operator*=;
    using Base::operator/=;
};

template <typename ValueType, size_t N>
struct _container_underlying_type
{
    typedef boost::multi_array_ref<ValueType, N-1> type;
};

template <typename ValueType>
struct _container_underlying_type<ValueType,1>
{
    typedef ValueType type;
};

}; // end of namespace GFTools
#endif

#include "Container.hxx"
