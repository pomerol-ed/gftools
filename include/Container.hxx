#ifndef ___GFTOOLS_CONTAINER_HPP___
#define ___GFTOOLS_CONTAINER_HPP___

#include "Container.hpp"
#include <fstream>

namespace GFTools {

/*
template <typename ValueType, size_t N, typename BoostContainerType> 
template<typename U>
inline ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> ContainerBase<ValueType,N,BoostContainerType>::operator[](size_t i)
{
    std::array<size_t, N-1> shape_tail;
    std::copy(_data.shape()+1, _data.shape()+N, shape_tail.begin());
    boost::multi_array_ref<ValueType, N-1> out(_data[i].origin(),shape_tail);
    return ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>>(out);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template<typename U>
inline const ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>> ContainerBase<ValueType,N,BoostContainerType>::operator[](size_t i) const
{
    std::array<size_t, N-1> shape_tail;
    std::copy(_data.shape()+1, _data.shape()+N, shape_tail.begin());
    boost::multi_array_ref<ValueType, N-1> out(_data[i].origin(),shape_tail);
    return ContainerBase<ValueType, N-1, boost::multi_array_ref<ValueType, N-1>>(out);
}
*/

template <typename ValueType, size_t N, typename BoostContainerType> 
std::ostream& operator<<(std::ostream& lhs, const ContainerBase<ValueType,N, BoostContainerType> &in)
{
    lhs << "[";
    std::ostream_iterator<decltype(in[0])> out_it (lhs,", ");
    for (size_t i=0; i<in._data.size(); ++i) {*out_it = in[i]; out_it++;};
    //std::copy(in._data.begin(),in._data.end(),out_it);
    //std::copy(in.begin(),in.end(),out_it);
    lhs << "]";
    return lhs;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template<typename N2>
MatrixType<ValueType> ContainerBase<ValueType,N, BoostContainerType>::getAsMatrix() const
{
    Eigen::Map<MatrixType<ValueType>> Map1 (_data.origin(), _data.shape()[0],_data.shape()[1]);
    return Map1;
}

template <typename ValueType, size_t N> 
template<typename U>
Container<ValueType,N>::Container(MatrixType<ValueType> rhs):
    Container<ValueType,N>(std::array<size_t,2>({{static_cast<size_t>(rhs.rows()), static_cast<size_t>(rhs.cols()) }}))
{
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), _data.origin());
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template<typename N2>
ContainerBase<ValueType,N, BoostContainerType>& ContainerBase<ValueType,N, BoostContainerType>::operator=(MatrixType<ValueType> rhs)
{
    assert(rhs.rows() == _data.shape()[0] && rhs.cols() == _data.shape()[1]);
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), _data.origin());
    return *this;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename U>
ContainerBase<ValueType,N, BoostContainerType> ContainerBase<ValueType,N, BoostContainerType>::conj()
{
    ContainerBase <ValueType,N,BoostContainerType> out(*this);
    //for (iterator it1 = this->begin(); it1!=this->end(); it1++) {  
    //    *it1=it1->conj();
    //}
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
ContainerBase<ValueType,N, BoostContainerType>& ContainerBase<ValueType,N, BoostContainerType>::operator=(ValueType rhs)
{
    std::fill(_data.origin(), _data.origin()+_data.num_elements(), rhs);
    return *this;
}


template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>>
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator+=(const R &rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    EigenMap map2(rhs._data.origin(),rhs._data.num_elements());
    map1+=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>>
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator+=(const R2& rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    map1+=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator+(const R &rhs) const
{
    Container<ValueType,N> out(*this); 
    out+=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator+(const R2& rhs) const
{
    Container<ValueType,N> out(*this); 
    out+=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>>
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator-=(const R &rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    EigenMap map2(rhs._data.origin(),rhs._data.num_elements());
    map1-=map2;
    return (*this);
}
 
template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator-=(const R2& rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    map1-=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator-(const R &rhs) const
{
    Container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator-(const R2& rhs) const
{
    Container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>>
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator*=(const R &rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    EigenMap map2(rhs._data.origin(),rhs._data.num_elements());
    map1*=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator*=(const R2& rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    map1*=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator*(const R &rhs) const
{
    Container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator*(const R2& rhs) const
{
    Container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}




template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>>
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator/=(const R &rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    EigenMap map2(rhs._data.origin(),rhs._data.num_elements());
    map1/=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
ContainerBase<ValueType,N,BoostContainerType>& ContainerBase<ValueType,N,BoostContainerType>::operator/=(const R2& rhs)
{
    EigenMap map1(_data.origin(),_data.num_elements());
    map1/=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R, ContainerBase<ValueType,N,BoostContainerType>::isContainer<R>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator/(const R &rhs) const
{
    Container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename R2, ContainerBase<ValueType,N,BoostContainerType>::isValue<R2>> 
Container<ValueType,N> ContainerBase<ValueType,N,BoostContainerType>::operator/(const R2& rhs) const
{
    Container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename U>
MatrixType<ValueType> ContainerBase<ValueType,N,BoostContainerType>::getAsDiagonalMatrix() const
{
    size_t size1 = _data.shape()[0];
    Eigen::DiagonalMatrix<ValueType, Eigen::Dynamic> out(size1);
    Eigen::Map<const VectorType<ValueType>> v(_data.origin(), size1);
    out.diagonal() = v;
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename U>
VectorType<ValueType> ContainerBase<ValueType,N,BoostContainerType>::getAsVector() const
{
    size_t size1 = _data.shape()[0];
    VectorType<ValueType> out(size1);
    Eigen::Map<const VectorType<ValueType>> v(_data.origin(), size1);
    out = v;
    return out;
}

} // end of namespace GFTools
#endif
