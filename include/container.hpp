#pragma once

//#define BOOST_RESULT_OF_USE_DECLTYPE
#include <type_traits>
#include <array>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif

#include <boost/multi_array.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace gftools { 

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base;

template <typename ValueType, size_t N>
struct container;

template<typename ValueType, size_t N, size_t M = N>
using container_view = container_base<real_type, N, typename container<real_type, N>::boost_t::template array_view<M>::type>;

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_traits
{
    typedef BoostContainerType boost_t;
    typedef decltype(std::declval<boost_t>()[0]) boost_under_type;
    typedef container_base<ValueType,N-1,boost_under_type> type;
    typedef container_base<ValueType,N-1,boost_under_type> ref_type;
};

template <typename ValueType, typename BoostContainerType>
struct container_traits<ValueType,1,BoostContainerType>
{

    typedef BoostContainerType boost_t;
    typedef decltype(std::declval<boost_t>()[0]) boost_under_type;
    typedef ValueType type;
    typedef ValueType& ref_type;
};


// ================================================================== //

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base 
{
    constexpr static size_t N_ = N;

    typedef ValueType value_t;
    typedef BoostContainerType boost_t;
    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    typedef Eigen::Map<EigenArray> EigenMap;
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, 1> VectorType;

    typedef typename container_traits<ValueType,N, BoostContainerType>::type under_type;
    typedef typename container_traits<ValueType,N, BoostContainerType>::ref_type under_ref_type;
    typedef typename container_traits<ValueType,N, BoostContainerType>::boost_under_type boost_under_type;

    typedef std::function<under_type(boost_under_type)> action_type;
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> iterator; 
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> const_iterator; 

    container_base(const boost_t &in):data_(in){};
    /** Copy constructor. */
    container_base(const container_base<ValueType,N,boost_t> &rhs):data_(rhs.data_){};
    container_base& operator=(const container_base<ValueType,N,boost_t> &rhs){data_ = rhs.data_; return (*this);};
    template <class B2> container_base& operator=(const container_base<ValueType,N,B2> &rhs){data_ = rhs.data_; return (*this);};
    /** Move constructor. */
    container_base(container_base<ValueType,N,boost_t> &&rhs):data_(rhs.data_){};
    container_base& operator=(container_base<ValueType,N,boost_t> &&rhs){std::swap(data_,rhs.data_); return (*this);};

    template <size_t N2 = N, typename DT = boost_t>
    using Is1d = typename std::enable_if<N2==1,DT>::type;
    template <size_t N2 = N, typename DT = boost_t>
    using IsNot1d = typename std::enable_if<N2!=1,DT>::type;
    template <size_t N2 = N, typename DT = boost_t>
    using Is2d = typename std::enable_if<N2==2,DT>::type;

    /** Access operators. */
    auto operator[](size_t i) -> under_ref_type { return under_ref_type(data_[i]); }
    auto operator[](size_t i) const -> under_ref_type const { return under_ref_type(data_[i]); }

    /** Return operators. */
    template<size_t N2 = N, typename U = Is2d<N2>>
    container_base<ValueType,N,boost_t>& operator=(MatrixType rhs);
    template<size_t N2 = N, typename U = Is2d<N2>>
    MatrixType as_matrix() const;
    /** Return a diagonal matrix, corresponding to the object */
    MatrixType as_diagonal_matrix() const;
    /** Return a vector, corresponding to the object. */
    VectorType as_vector() const;
 
    /** Conjugate. */
    template <typename T = ValueType, typename U = typename std::enable_if<std::is_convertible<T, complex_type>::value, int>::type> 
        container<ValueType,N> conj();
    /** Sum of all values in the container. */
    ValueType sum(){return EigenMap(data_.origin(), data_.num_elements()).sum();};

    /** Make the object streamable. */
    template <typename V1, size_t M, typename B>
    friend std::ostream& operator<<(std::ostream& lhs, const container_base<V1, M, B> &in);

    /** Begin iterator. */
    const_iterator begin() const;
    /** End iterator. */
    const_iterator end() const;
    
    /** Mathematical operators. */
    template <typename T> 
        typename std::enable_if<std::is_convertible<T,ValueType>::value,container_base&>::type operator=(T rhs);// { Base(*this) = rhs; return (*this);}

    template <typename T>
    using BaseRefIfContainer = typename std::enable_if<T::N_>=1, container_base<ValueType,N,boost_t>&>::type;
    template <typename T>
    using ContainerIfContainer = typename std::enable_if<T::N_>=1, container<ValueType,N>>::type;
    template <typename T>
    using BaseRefIfValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, container_base<ValueType,N,boost_t>&>::type;
    template <typename T>
    using ContainerIfValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, container<ValueType,N>>::type;

    template <typename R>
        BaseRefIfContainer<R> operator+=(const R &rhs); 
    template <typename R> 
        ContainerIfContainer<R> operator+(const R &rhs) const;
    template <typename R2>
        BaseRefIfValue<R2> operator+=(const R2& rhs);
    template <typename R2> 
        ContainerIfValue<R2> operator+(const R2& rhs) const;

    template <typename R>
        BaseRefIfContainer<R> operator-=(const R &rhs); 
    template <typename R> 
        ContainerIfContainer<R> operator-(const R &rhs) const;
    template <typename R2> 
        BaseRefIfValue<R2> operator-=(const R2& rhs);
    template <typename R2> 
        ContainerIfValue<R2> operator-(const R2& rhs) const;

    template <typename R>
        BaseRefIfContainer<R> operator*=(const R &rhs); 
    template <typename R> 
        ContainerIfContainer<R> operator*(const R &rhs) const;
    template <typename R2> 
        BaseRefIfValue<R2> operator*=(const R2& rhs);
    template <typename R2> 
        ContainerIfValue<R2> operator*(const R2& rhs) const;

    template <typename R>
        BaseRefIfContainer<R> operator/=(const R &rhs); 
    template <typename R> 
        ContainerIfContainer<R> operator/(const R &rhs) const;
    template <typename R2> 
        BaseRefIfValue<R2> operator/=(const R2& rhs);
    template <typename R2> 
        ContainerIfValue<R2> operator/(const R2& rhs) const;

    friend container<ValueType,N> operator* (const ValueType & lhs, const container_base<ValueType,N,boost_t> & rhs) {return rhs*lhs;};
    friend container<ValueType,N> operator+ (const ValueType & lhs, const container_base<ValueType,N,boost_t> & rhs) {return rhs+lhs;};
    friend container<ValueType,N> operator- (const ValueType & lhs, const container_base<ValueType,N,boost_t> & rhs) {return rhs*(-1.0)+lhs;};
    friend container<ValueType,N> operator/ (const ValueType & lhs, const container_base<ValueType,N,boost_t> & rhs) {container<ValueType,N> out(rhs); out=lhs; return out/rhs;};
    
    
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
    
//-----------------------------//    
    mutable boost_t data_;
};

template <typename ValueType, size_t N>
struct container : container_base<ValueType,N,typename boost::multi_array<ValueType, N>> {
    typedef boost::multi_array<ValueType, N> boost_t;
    typedef container_base<ValueType,N,boost_t> Base;
    using Base::data_;
    typedef typename Base::MatrixType MatrixType;

    container(std::array<size_t,N> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};

    template <typename CT>
        container(container_base<ValueType,N,CT> in) : container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(in.data_) {};

    template<typename ...ShapeArgs>
        container(ShapeArgs...in):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(std::array<int,N>({{in...}}))) {
            static_assert(sizeof...(in) == N,"arg mismatch");
        };

    template<size_t N2 = N, typename U = typename std::enable_if<N2==2, bool>> 
        container<ValueType,N> (MatrixType rhs);
    using Base::operator+=;
    using Base::operator-=;
    using Base::operator*=;
    using Base::operator/=;
    //using Base::operator=;
    template <typename T> typename std::enable_if<std::is_convertible<T,ValueType>::value,container&>::type operator=(T rhs) { 
        typename std::add_lvalue_reference<Base>::type(*this) = rhs; return (*this);}
};


}; // end of namespace gftools

#include "container.hxx"
