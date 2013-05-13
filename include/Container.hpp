#pragma once
#ifndef ___GFTOOLS_CONTAINER_H___
#define ___GFTOOLS_CONTAINER_H___

#include "Defaults.hpp"
#include "Tools.hpp"
#include <type_traits>

#include <boost/multi_array.hpp>
#include <boost/operators.hpp>
#include <boost/iterator/permutation_iterator.hpp>

namespace GFTools { 

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase;

/** This class is a N-dimensional view of container to store the data of
 * the ValueType type. */
//template <typename ValueType, size_t N, size_t M=N>
//struct ContainerView;

/** This class is a N-dimensional container to store the data of
 * the ValueType type. */
//template <typename ValueType, size_t N>
//struct Container;

// ================================================================== //

template <typename ValueType, size_t N, typename BoostContainerType>
struct ContainerBase : boost::operators<ContainerBase<ValueType,N,BoostContainerType>>
{
    template <size_t M> using boost_view_type = typename boost::multi_array<ValueType, N>::template array_view<M>::type;
    typedef boost::multi_array<ValueType, N> boost_array_type;
    typedef boost::multi_array_ref<ValueType, N> boost_array_ref_type;

    typedef Eigen::Array<ValueType, Eigen::Dynamic, 1> EigenArray;
    typedef Eigen::Map<EigenArray> EigenMap;

    mutable BoostContainerType _data;

    ContainerBase(BoostContainerType in):_data(in){};
    template <typename ...Args> ContainerBase(Args... in):_data(BoostContainerType(in...)){};
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

    template<typename N2 = typename std::enable_if<N==2, bool>> MatrixType<ValueType> getAsMatrix() const;
    template<typename N2 = typename std::enable_if<N==2, bool>> ContainerBase<ValueType,N,BoostContainerType> (MatrixType<ValueType> &&rhs);
    template<typename N2 = typename std::enable_if<N==2, bool>> ContainerBase<ValueType,N,BoostContainerType>& operator=(MatrixType<ValueType> &&rhs);
 
    /** Conjugate. */
    //template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type=0> 
    //    Container<N, ValueType> conj();
    /** Recursively iterates and sums all values in the container. */
    //ValueType sum();

    /** Mathematical operators. */
    ContainerBase& operator=(ValueType rhs){}
    template <typename OtherContainerType>
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const ContainerBase<ValueType,N,OtherContainerType> &rhs); 
    template <typename OtherContainerType> 
        ContainerBase<ValueType,N,BoostContainerType> operator+(const ContainerBase<ValueType,N,OtherContainerType> &rhs) const;
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType>& operator+=(const RhsArg& rhs);
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType> operator+(const RhsArg& rhs) const;

    template <typename OtherContainerType>
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const ContainerBase<ValueType,N,OtherContainerType> &rhs); 
    template <typename OtherContainerType> 
        ContainerBase<ValueType,N,BoostContainerType> operator-(const ContainerBase<ValueType,N,OtherContainerType> &rhs) const;
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType>& operator-=(const RhsArg& rhs);
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType> operator-(const RhsArg& rhs) const;

    template <typename OtherContainerType>
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const ContainerBase<ValueType,N,OtherContainerType> &rhs); 
    template <typename OtherContainerType> 
        ContainerBase<ValueType,N,BoostContainerType> operator*(const ContainerBase<ValueType,N,OtherContainerType> &rhs) const;
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType>& operator*=(const RhsArg& rhs);
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType> operator*(const RhsArg& rhs) const;

    template <typename OtherContainerType>
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const ContainerBase<ValueType,N,OtherContainerType> &rhs); 
    template <typename OtherContainerType> 
        ContainerBase<ValueType,N,BoostContainerType> operator/(const ContainerBase<ValueType,N,OtherContainerType> &rhs) const;
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType>& operator/=(const RhsArg& rhs);
    template <typename RhsArg> 
        ContainerBase<ValueType,N,BoostContainerType> operator/(const RhsArg& rhs) const;
    
    //ContainerBase<ValueType,N,BoostContainerType>& operator=(const ValueType &rhs);
    /** Make the object streamable. */
    template <typename V1, size_t M, typename B>
    friend std::ostream& operator<<(std::ostream& lhs, const ContainerBase<V1, M, B> &in);

    /** Begin iterator. */
    //const_iterator begin() const;
    /** End iterator. */
    //const_iterator end() const;

    class exWrongIndex : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
};


template <typename ValueType, size_t N>
using Container = ContainerBase<ValueType, N, typename boost::multi_array<ValueType, N>>;

/*template <typename ValueType, size_t N>
struct Container : public ContainerBase<ValueType, N, typename boost::multi_array<ValueType, N>>
{
    typedef boost::multi_array_ref<ValueType, N> DataType;
    typedef ContainerBase<ValueType, N, boost::multi_array_ref<ValueType, N>> BaseType;
    DataType _data;
}
*/

/*
template <typename ValueType, size_t N, size_t M>
struct ContainerView : public ContainerBase<ValueType, N, typename boost::multi_array<ValueType, N>::template array_view<M>::type > {
    typedef typename boost::multi_array<ValueType, N>::template array_view<M>::type boost_view_type;
    //using ContainerBase<ValueType, N, boost_view_type>::_data;
    boost_view_type _data;
    template<typename U = typename std::enable_if<N!=1>::type> ContainerView<ValueType, N-1>& operator[](size_t i);
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i);
    ContainerView(boost_view_type in);
    using ContainerBase<ValueType, N, typename boost::multi_array<ValueType, N>::template array_view<M>::type>::exWrongIndex;
};

template <typename ValueType, size_t N>
struct ContainerRef : public ContainerBase<ValueType, N, boost::multi_array_ref<ValueType, N>>
{
    typedef boost::multi_array_ref<ValueType, N> DataType;
    typedef ContainerBase<ValueType, N, boost::multi_array_ref<ValueType, N>> BaseType;
    DataType _data;
    ContainerRef(typename boost::multi_array<ValueType,N+1>::reference in);
    template<typename U = typename std::enable_if<N!=1>::type> ContainerRef<ValueType, N-1> operator[](size_t i);
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i);
    using typename BaseType::exWrongIndex;
};

*/
//template <typename ValueType, size_t N>
//struct SubContainer : public ContainerBase<ValueType, N, typename boost::detail::multi_array::template sub_array<ValueType,N>>
//{
    /*typedef boost::multi_array_ref<ValueType, N> DataType;
    typedef ContainerBase<ValueType, N, boost::multi_array_ref<ValueType, N>> BaseType;
    DataType _data;
    ContainerRef(typename boost::multi_array<ValueType,N+1>::reference in);
    template<typename U = typename std::enable_if<N!=1>::type> ContainerRef<ValueType, N-1> operator[](size_t i);
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i);
    using typename BaseType::exWrongIndex;
*/
//};

/**
template <typename ValueType, size_t N>
struct Container : public ContainerBase<ValueType, N, boost::multi_array<ValueType, N>> {
//static_assert(N>1, "N=1 has a specialized variant");
    friend struct Container<ValueType,N+1>;
    typedef boost::multi_array<ValueType, N> DataType;
    typedef ContainerBase<ValueType, N, DataType> BaseType;

    DataType _data;
    
    template<typename U = typename std::enable_if<N!=1>::type> SubContainer<ValueType, N-1> operator[](size_t i);
    template<typename U = typename std::enable_if<N==1>::type> ValueType& operator[](size_t i);
    //std::array<size_t, N> _shape;

    // Descriptors
    
public:
    //typedef typename std::vector<Container<N-1, ValueType>>::iterator iterator;
    //typedef typename std::vector<Container<N-1, ValueType>>::const_iterator const_iterator;
    *//** Constructor from the std::array of size_t. */
    //Container ( const std::array<size_t, N> &in);
    /** Copy constructor. */
    //Container(const ContainerBase<ValueType,N,BoostContainerType> &rhs):_vals(rhs._vals){};
    //template<std::enable_if<N==1, bool>::type = 0> Container<ValueType, N-1> operator[](size_t i);
        /** Returns the size of the _vals. */
    //size_t getSize() const { return _vals.size(); };
    /** Algebraic operators. */
    //ContainerBase<ValueType,N,BoostContainerType>& operator=(ContainerBase<ValueType,N,BoostContainerType> &&rhs);
    //ContainerBase<ValueType,N,BoostContainerType>& operator=(const ContainerBase<ValueType,N,BoostContainerType> &rhs);
    //ContainerBase<ValueType,N,BoostContainerType>& operator=(const ValueType &rhs);
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType>& operator+=(const RhsArg &rhs); 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType> operator+(const RhsArg &rhs) const; 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType>& operator*=(const RhsArg &rhs); 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType> operator*(const RhsArg &rhs) const; 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType>& operator/=(const RhsArg &rhs); 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType> operator/(const RhsArg &rhs) const; 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType>& operator-=(const RhsArg &rhs); 
    //template <typename RhsArg> ContainerBase<ValueType,N,BoostContainerType> operator-(const RhsArg &rhs) const; 
    //ContainerBase<ValueType,N,BoostContainerType>& operator+=(const ContainerBase<ValueType,N,BoostContainerType> &rhs); 
    //ContainerBase<ValueType,N,BoostContainerType>& operator*=(const ContainerBase<ValueType,N,BoostContainerType> &rhs); 
    //ContainerBase<ValueType,N,BoostContainerType>& operator/=(const ContainerBase<ValueType,N,BoostContainerType> &rhs); 
    //friend inline ContainerBase<ValueType,N,BoostContainerType> operator* (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs*lhs;};
    //friend inline ContainerBase<ValueType,N,BoostContainerType> operator+ (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs+lhs;};
    //friend inline ContainerBase<ValueType,N,BoostContainerType> operator- (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {return rhs*(-1.0)+lhs;};
    //friend inline ContainerBase<ValueType,N,BoostContainerType> operator/ (const ValueType & lhs, const ContainerBase<ValueType,N,BoostContainerType> & rhs) {ContainerBase<ValueType,N,BoostContainerType> out(rhs); out=lhs; return out/rhs;};


    
//    using typename ContainerBase<ValueType, N, boost::multi_array<ValueType, N>>::exWrongIndex;
//};

}; // end of namespace GFTools
#endif

#include "Container.hxx"
