#ifndef ___GFTOOLS_CONTAINER_HPP___
#define ___GFTOOLS_CONTAINER_HPP___

#include "container.hpp"
#include <fstream>

namespace gftools {

template <typename ValueType, size_t N, typename BoostCType> 
std::ostream& operator<<(std::ostream& lhs, const container_base<ValueType,N, BoostCType> &in)
{
    lhs << "[";
    std::ostream_iterator<decltype(in[0])> out_it (lhs,", ");
    for (size_t i=0; i<in.data_.size(); ++i) {*out_it = in[i]; out_it++;};
    lhs << "]";
    return lhs;
}

template <typename ValueType, size_t N, typename BoostCType> 
template<size_t N2, typename U>
    typename container_base<ValueType,N, BoostCType>::MatrixType 
    container_base<ValueType,N, BoostCType>::as_matrix() const
{
    Eigen::Map<MatrixType> Map1 (data_.origin(), data_.shape()[0],data_.shape()[1]);
    return Map1;
}

template <typename ValueType, size_t N> 
template <size_t, typename>
container<ValueType,N>::container(MatrixType rhs):
    container<ValueType,N>(std::array<size_t,2>({{static_cast<size_t>(rhs.rows()), static_cast<size_t>(rhs.cols()) }}))
{
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), data_.origin());
}

template <typename ValueType, size_t N, typename BoostCType> 
template<size_t N2, typename U>
container_base<ValueType,N, BoostCType>&
    container_base<ValueType,N, BoostCType>::operator=(MatrixType rhs)
{
    assert(rhs.rows() == data_.shape()[0] && rhs.cols() == data_.shape()[1]);
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), data_.origin());
    return *this;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType,N, BoostCType>::MatrixType 
    container_base<ValueType,N,BoostCType>::as_diagonal_matrix() const
{
    size_t size1 = data_.size();
    Eigen::DiagonalMatrix<ValueType, Eigen::Dynamic> out(size1);
    Eigen::Map<const VectorType> v(data_.origin(), size1);
    out.diagonal() = v;
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType,N, BoostCType>::VectorType container_base<ValueType,N,BoostCType>::as_vector() const
{
    size_t size1 = data_.size();
    VectorType out(size1);
    Eigen::Map<const VectorType> v(data_.origin(), size1);
    out = v;
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename T, typename U>
container<ValueType,N> container_base<ValueType,N, BoostCType>::conj()
{
    container <ValueType,N> out(*this);
    std::transform(out.data_.data(), out.data_.data() + out.data_.size(), out.data_.data(), [](ValueType in){return std::conj(in);});
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename T> 
    typename std::enable_if<std::is_convertible<T,ValueType>::value,container_base<ValueType,N,BoostCType>&>::type 
    container_base<ValueType,N,BoostCType>::operator=(T rhs)
{
    std::fill(data_.origin(), data_.origin()+data_.num_elements(), rhs);
    return *this;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::const_iterator container_base<ValueType, N, BoostCType>::begin() const 
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(data_.begin(),f1);
};

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::const_iterator container_base<ValueType, N, BoostCType>::end() const 
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(data_.end(),f1);
}


template <typename ValueType, size_t N, typename BoostCType> 
template <typename R> 
typename container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator+=(const R &rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    EigenMap map2(rhs.data_.origin(),rhs.data_.num_elements());
    map1+=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2> 
container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator+=(const R2& rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    map1+=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator+(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out+=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator+(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out+=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator-=(const R &rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    EigenMap map2(rhs.data_.origin(),rhs.data_.num_elements());
    map1-=map2;
    return (*this);
}
 
template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator-=(const R2& rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    map1-=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator-(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator-(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator*=(const R &rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    EigenMap map2(rhs.data_.origin(),rhs.data_.num_elements());
    map1*=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator*=(const R2& rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    map1*=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator*(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2> 
container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator*(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator/=(const R &rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    EigenMap map2(rhs.data_.origin(),rhs.data_.num_elements());
    map1/=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator/=(const R2& rhs)
{
    EigenMap map1(data_.origin(),data_.num_elements());
    map1/=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator/(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator/(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

} // end of namespace GFTools
#endif
