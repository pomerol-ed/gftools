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

/// container_base is a base class, that provides an interface of multidimensional array 
/// with arithmetic operations (element-wise addition, multiplication, etc)
/// it is a wrapper over boost::multi_array and has 3 types : type of the stored numbers, rank and boost::multi_array type 
/// TODO ValueType is a duplication here
template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base;

/// container is a multidimensional array with math operations, that allocates and stores memory
template <typename ValueType, size_t N>
struct container;

/// container_view is a "view" of a container_base that allows to access elements with a different storage order
template<typename ContainerType, size_t M = ContainerType::N_>
using container_view = container_base<typename ContainerType::value_type, ContainerType::N_, typename ContainerType::boost_t::template array_view<M>::type>;

/// container_ref is a container without allocation of memory, i.e. it serves as an interface to an outside chunk of memory with all container operations
template <typename ValueType, size_t N>
using container_ref = container_base<ValueType, N, boost::multi_array_ref<ValueType,N>>;

/// container_traits is a helper ''traits'' class that provides type definitions
template <typename, bool View=false> struct container_traits; 

template <typename ValueType, size_t N, typename BoostContainerType>
struct container_base 
{
    /// total rank (number of dimensions) of the container_base
    constexpr static size_t N_ = N;
    /// a helper flag that is true, when the container_base is a view
    constexpr static bool is_view_ = (N != BoostContainerType::dimensionality);

    // typedefs
    /// typedef for stored values
    typedef ValueType value_type;
    /// typedef for a result of const operator[] operation (i.e. can be an array or a number) 
    typedef typename container_traits<container_base, is_view_>::type under_type;
    /// typedef for a result of operator[] operation (i.e. can be an array or a number) 
    typedef typename container_traits<container_base, is_view_>::ref_type under_ref_type;
    /// typedef for a resukt of operator[] of underlying boost::multi_array type
    typedef typename container_traits<container_base, is_view_>::boost_under_type boost_under_type;
    /// typedef wrapped boost::multi_array 
    typedef BoostContainerType boost_t;
    /// typedef for a flattened array (Eigen::Array type) 
    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    /// typedef for an Eigen::Map object that provides arithmetic operations on flattened arrays
    typedef Eigen::Map<EigenArray> EigenMap;
    /// typedef for a matrix with value_type values
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    /// typedef for a flattened Eigen::Vector, that provides vector*matrix operations
    typedef Eigen::Matrix<ValueType,Eigen::Dynamic, 1> VectorType;

    // helper type structues
    /// TODO : move to traits
    /// Is2d defines true_type if the objectr is 2d -> allows matrix operations
    template <size_t N2 = N, typename DT = boost_t>
    using Is2d = typename std::enable_if<N2==2,DT>::type;
    // helper typedefs for math operations (to be able to multiply container and a ref, container and container, container and a number etc)
    /// BaseRefIfContainer defines container_base type if the template parameter T is a container_base 
    template <typename T>
    using BaseRefIfContainer = typename std::enable_if<T::N_>=1, container_base<ValueType,N,boost_t>&>::type;
    /// ContaierIfContainer defines container type if the template parameter T is a container_base 
    template <typename T>
    using ContainerIfContainer = typename std::enable_if<T::N_>=1, container<ValueType,N>>::type;
    /// BaseRefIfValue defines container_base if the template parameter T is a number
    template <typename T>
    using BaseRefIfValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, container_base<ValueType,N,boost_t>&>::type;
    /// BaseRefIfValue defines container if the template parameter T is a number
    template <typename T>
    using ContainerIfValue = typename std::enable_if<std::is_convertible<T,ValueType>::value, container<ValueType,N>>::type;

    // iterators
    /// typedef for an operation that converts a operator[] of wrapped multi_array into operator[] of container
    typedef std::function<under_type(boost_under_type)> action_type;
    /// iterator (operates as a wrapper over multi_array)
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> iterator; 
    /// const_iterator
    typedef boost::transform_iterator<action_type,typename boost_t::iterator> const_iterator; 

    // constructors
    /// Constructor from the boost::multi_array (view, ref) type
    container_base(const boost_t &in):storage_(in){};
    /// Copy constructor
    container_base(const container_base<ValueType,N,boost_t> &rhs):storage_(rhs.storage_){};
    /// copy constructor from a different container_base
    template <typename BT2>
        container_base(const container_base<ValueType,N,BT2> &rhs):storage_(rhs.boost_container_()){};
    /// move constructor (C++11)
    container_base(container_base<ValueType,N,boost_t> &&rhs):storage_(std::forward<boost_t>(rhs.storage_)){ }
    /// swap operation
    void swap(container_base& rhs) {std::swap(storage_,rhs.storage_);}

    // assigments
    /// assignment operator from a differenet container base
    template <typename BT2>
        container_base& operator=(const container_base<ValueType,N,BT2> &rhs){storage_ = rhs.boost_container_(); return (*this);};
    /// assignment operator
    container_base& operator=(const container_base<ValueType,N,boost_t> &rhs){storage_ = rhs.storage_; return (*this);};
    /// move assignment (C++11)
    container_base& operator=(container_base<ValueType,N,boost_t> &&rhs){std::swap(storage_,rhs.storage_); return (*this);};
    /// assign from a number
    template <typename T> 
        typename std::enable_if<std::is_convertible<T,ValueType>::value,container_base&>::type operator=(T rhs);// { Base(*this) = rhs; return (*this);}
    /// assign from matrix (2d containers only)
    template<size_t N2 = N, typename U = Is2d<N2>>
    container_base<ValueType,N,boost_t>& operator=(MatrixType rhs);

    // access operations
    /// access rank N-1 object
    auto operator[](size_t i) -> under_ref_type { return under_ref_type(storage_[i]); }
    /// const access rank N-1 object
    auto operator[](size_t i) const -> under_ref_type const { return under_ref_type(storage_[i]); }
    /// return value ref from an array of indices
    ValueType& operator()(std::array<size_t, N> indices){return storage_(indices);}
    /// return value const-ref from an array of indices
    const ValueType& operator() (std::array<size_t, N> indices) const {return storage_(indices);}

    // view modifiers
    /// copy into a matrix 
    /// TODO consider viewing using Eigen::Map instead of copying
    template<size_t N2 = N, typename U = Is2d<N2>>
    MatrixType as_matrix() const;
    /// copy into a diagonal matrix
    MatrixType as_diagonal_matrix() const;
    /// copy flattened container into a vector
    VectorType as_vector() const;

    // global operations
    /// return a reference to this container as a flattened array 
    container_base<ValueType,1,boost::multi_array_ref<ValueType,1>> flatten();
    /// conjugate copy for complex number objects
    template<typename T=ValueType> typename std::enable_if< std::is_same<T, complex_type>::value, container<ValueType,N>>::type conj() const;
    /// conjugate copy non-complex number objects (does nothing)
    template<typename T=ValueType> typename std::enable_if<!std::is_same<T, complex_type>::value, container<ValueType,N>>::type conj() const;
    /// return transposed copy (2d only)
    template<size_t N2 = N, typename U = Is2d<N2>>
    container<ValueType,N> transpose() const { return container<ValueType, N>( this->as_matrix().transpose()); }
    /// return hermite conjugate copy (2d only)
    template<size_t N2 = N, typename U = Is2d<N2>>
    container<ValueType,N> hermite_conj() const { return this->transpose().conj(); }

    /// Sum of all values in the container.
    ValueType sum() const; 
    /// return the total size (sum of all dimensions) of the container
    int size() const;
    /// returns the shape of the container
    std::array<size_t,N> shape() const;
    /// returns the squared norm difference between 2 containers 
    template <typename V2, typename DT> double diff(const container_base<V2,N,DT>& r, bool norm = true) const; 

    /// returns a pointer to the first element in the container
    ValueType* data() { return storage_.origin(); }
    /// returns a const pointer to the first element in the container
    const ValueType* data() const { return storage_.origin(); }
    /** Make the object streamable. */
    template <typename V1, size_t M, typename B>
    friend std::ostream& operator<<(std::ostream& lhs, const container_base<V1, M, B> &in);

    // iterator accessors
    /// Begin const iterator./
    const_iterator begin() const;
    /// Begin iterator./
    iterator begin();
    /// End const iterator.
    const_iterator end() const;
    /// End iterator.
    iterator end();
    
    // Mathematical operations
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
    friend container<ValueType,N> operator/ (const ValueType & lhs, const container_base<ValueType,N,boost_t> & rhs) {
        container<ValueType,N> out(rhs); out=lhs; return out/rhs;};
    
    /// An exception provided for incorrect indices 
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
    
    /// return underlying boost container
    /// underscore is added to discourage from using the method
    boost_t& boost_container_() const {return storage_; }
protected:
    /// allow other container t access internal storage
    friend struct container<ValueType,N>;
    /// wrapped boost multi_array
    // it is mutable, because Eigen::Map operates with * pointers
    // FIXME - remove mutable
    mutable boost_t storage_;
};

template <typename ValueType, size_t N>
struct container : container_base<ValueType,N,typename boost::multi_array<ValueType, N>> {
    typedef boost::multi_array<ValueType, N> boost_t;
    typedef container_base<ValueType,N,boost_t> Base;
    using Base::storage_;
    typedef typename Base::MatrixType MatrixType;

    /// IsNotContainer defines true_type if the template parameter T is not a container_base type
    template <typename T>
    using IsNotContainer = typename std::enable_if<!(T::N_>=1)>::type;

    /// construct container from a given shape of ints (initialize with zeros)
    container(std::array<int,N> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    /// construct container from a given shape of size_t (initialize with zeros) 
    explicit container(std::array<size_t,N> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(shape)) {};
    /// construct container from an initializer list of ints, aka container<double, 2> a({{1,2}})
    container(std::initializer_list<int> shape):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(std::array<int, N>(shape))) {};

    /// construct from a different container_base 
    // using value here is safe with a move constructor
    template <typename CT>
        container(container_base<ValueType,N,CT> in) : container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(in.storage_) {};

    /// construct container from given variable amount of numbers, aka container<double, 4> a(1,2,4,2)
    template<typename ...ShapeArgs,
        typename = typename std::enable_if<sizeof...(ShapeArgs) == N 
               && (std::is_same<std::tuple<ShapeArgs...>, typename tuple_tools::repeater<int,N>::tuple_type>::value // Arguments have to be strictly ints
               || std::is_same<std::tuple<ShapeArgs...>, typename tuple_tools::repeater<size_t,N>::tuple_type>::value) // or size_t
        ,int>::type>
        container(ShapeArgs...in):container_base<ValueType,N,typename boost::multi_array<ValueType, N>>(boost::multi_array<ValueType, N>(std::array<int,N>({{in...}}))) {
            static_assert(sizeof...(in) == N,"arg mismatch");
        };
    /// construct 2d container from matrix
    template<size_t N2 = N, typename U = typename std::enable_if<N2==2, bool>::type> 
        container<ValueType,N> (MatrixType rhs);

    // inherit math from base
    using Base::operator+=;
    using Base::operator-=;
    using Base::operator*=;
    using Base::operator/=;
    /// assign from a number
    template <typename T> typename std::enable_if<std::is_convertible<T,ValueType>::value,container&>::type operator=(T rhs) { 
        typename std::add_lvalue_reference<Base>::type(*this) = rhs; return (*this);}
};

/// type traits for container/container_base of rank>1
template <typename ValueType, size_t N, typename BoostContainerType>
struct container_traits<container_base<ValueType,N,BoostContainerType>,false>
{
    /// typedef boost multi_array type
    typedef BoostContainerType boost_t;
    /// typedef for a result of const operator[] operation (i.e. can be an array or a number) 
    typedef typename boost_t::reference boost_under_type;
    /// typedef for a result of operator[] operation (i.e. can be an array or a number) 
    typedef container_base<ValueType,N-1,boost_under_type> type;
    /// typedef for a resukt of operator[] of underlying boost::multi_array type
    typedef container_base<ValueType,N-1,boost_under_type> ref_type;
};

/// type traits for container/container_base of rank==1
template <typename ValueType, typename BoostContainerType>
struct container_traits<container_base<ValueType,1,BoostContainerType>,false>
{
    typedef BoostContainerType boost_t;
    typedef typename boost_t::reference boost_under_type;
    typedef ValueType type;
    typedef ValueType& ref_type;
};

/// type traits for container views
template <typename ValueType, size_t N, typename BoostViewType>
struct container_traits<container_base<ValueType,N,BoostViewType>,true>
{
    typedef BoostViewType boost_t;
    typedef typename boost_t::reference boost_under_type;
    typedef ValueType type;
    typedef ValueType& ref_type;
};



}; // end of namespace gftools

#include "container.hxx"
