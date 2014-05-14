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

#include "tuple_tools.hpp"

namespace gftools { 

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base;

template <typename ValueType, size_t N>
struct container;

template<typename ContainerType, size_t M = ContainerType::N_>
using container_view = container_base<typename ContainerType::value_type, ContainerType::N_, typename ContainerType::boost_t::template array_view<M>::type>;

template <typename ValueType, size_t N>
using container_ref = container_base<ValueType, N, boost::multi_array_ref<ValueType,N>>;

namespace extra {
template <typename, bool View=false> struct container_traits; 

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_traits<container_base<ValueType,N,BoostContainerType>,false>
{
    typedef BoostContainerType boost_t;
    typedef typename boost_t::reference boost_under_type;
    typedef container_base<ValueType,N-1,boost_under_type> type;
    typedef container_base<ValueType,N-1,boost_under_type> ref_type;
};

template <typename ValueType, typename BoostContainerType>
struct container_traits<container_base<ValueType,1,BoostContainerType>,false>
{
    typedef BoostContainerType boost_t;
    typedef typename boost_t::reference boost_under_type;
    typedef ValueType type;
    typedef ValueType& ref_type;
};


template <typename ValueType, size_t N, typename BoostViewType>
struct container_traits<container_base<ValueType,N,BoostViewType>,true>
{
    typedef BoostViewType boost_t;
    typedef typename boost_t::reference boost_under_type;
    typedef ValueType type;
    typedef ValueType& ref_type;
};

}; // end of namespace extra

// ================================================================== //

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base 
{
    constexpr static size_t N_ = N;
    constexpr static bool is_view_ = (N != BoostContainerType::dimensionality);

    typedef ValueType value_type;
    typedef typename extra::container_traits<container_base, is_view_>::type under_type;
    typedef typename extra::container_traits<container_base, is_view_>::ref_type under_ref_type;
    typedef typename extra::container_traits<container_base, is_view_>::boost_under_type boost_under_type;

    typedef BoostContainerType boost_t;
    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    typedef Eigen::Map<EigenArray> EigenMap;
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, 1> VectorType;

    typedef std::function<under_type(boost_under_type)> action_type;
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> iterator; 
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> const_iterator; 


    container_base(const boost_t &in):storage_(in){};
    /** Copy constructor. */
    container_base(const container_base<ValueType,N,boost_t> &rhs):storage_(rhs.storage_){};
    template <typename BT2>
        container_base(const container_base<ValueType,N,BT2> &rhs):storage_(rhs.boost_container_()){};
    template <typename BT2>
        container_base& operator=(const container_base<ValueType,N,BT2> &rhs){storage_ = rhs.boost_container_(); return (*this);};
    container_base& operator=(const container_base<ValueType,N,boost_t> &rhs){storage_ = rhs.storage_; return (*this);};
    /** Move constructor. */
    container_base(container_base<ValueType,N,boost_t> &&rhs):storage_(std::forward<boost_t>(rhs.storage_)){ }
    container_base& operator=(container_base<ValueType,N,boost_t> &&rhs){std::swap(storage_,rhs.storage_); return (*this);};
    void swap(container_base& rhs) {std::swap(storage_,rhs.storage_);}

    /** Access operators. */
    auto operator[](size_t i) -> under_ref_type { return under_ref_type(storage_[i]); }
    auto operator[](size_t i) const -> under_ref_type const { return under_ref_type(storage_[i]); }

    boost_t& boost_container_() const {return storage_; }
    ValueType& operator()(std::array<size_t, N> indices){return storage_(indices);}
    const ValueType& operator() (std::array<size_t, N> indices) const {return storage_(indices);}

    container_base<ValueType,1,boost::multi_array_ref<ValueType,1>> flatten();
// TODO : put views

    template <size_t N2 = N, typename DT = boost_t>
    using Is2d = typename std::enable_if<N2==2,DT>::type;
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
    //container<ValueType,N>&& conj();
    template<typename T=ValueType> typename std::enable_if< std::is_same<T, complex_type>::value, container<ValueType,N>>::type conj() const;
    template<typename T=ValueType> typename std::enable_if<!std::is_same<T, complex_type>::value, container<ValueType,N>>::type conj() const;
    //typename std::enable_if<!std::is_convertible<ValueType, complex_type>::value, container<ValueType,N>&&>::type conj_d();
    //container<ValueType,N>&& conj();
    /** Sum of all values in the container. */
    ValueType sum() const; 
    ValueType* data() { return storage_.origin(); }
    int size() const;
    std::array<size_t,N> shape() const;
    template <typename V2, typename DT>
    double diff(const container_base<V2,N,DT>& r) const 
        { return std::sqrt(std::abs( ((*this)*(*this).conj()).sum() + (r*r.conj()).sum() - ((*this)*r.conj()).sum() - ((*this).conj()*r).sum())); }

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
    protected:
    friend struct container<ValueType,N>;
    
    mutable boost_t storage_;
};

template <typename ValueType, size_t N>
struct container : container_base<ValueType,N,typename boost::multi_array<ValueType, N>> {
    typedef boost::multi_array<ValueType, N> boost_t;
    typedef container_base<ValueType,N,boost_t> Base;
    using Base::storage_;
    typedef typename Base::MatrixType MatrixType;

    template <typename T>
    using IsNotContainer = typename std::enable_if<!(T::N_>=1)>::type;

    container(std::array<int,N> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    explicit container(std::array<size_t,N> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    container(std::initializer_list<int> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(std::array<int, N>(shape))) {};

    template <typename CT>
        container(container_base<ValueType,N,CT> in) : container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(in.storage_) {};

    template<typename ...ShapeArgs,
        typename = typename std::enable_if<sizeof...(ShapeArgs) == N 
               && (std::is_same<std::tuple<ShapeArgs...>, typename tuple_tools::repeater<int,N>::tuple_type>::value // Arguments have to be strictly ints
               || std::is_same<std::tuple<ShapeArgs...>, typename tuple_tools::repeater<size_t,N>::tuple_type>::value) // or size_t
        ,int>::type>
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
