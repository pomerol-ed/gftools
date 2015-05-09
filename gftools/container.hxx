#pragma once

#include "container.hpp"
#include <fstream>

namespace gftools {

template <typename ValueType, size_t N, typename BoostCType> 
std::ostream& operator<<(std::ostream& lhs, const container_base<ValueType,N, BoostCType> &in)
{
    lhs << "[";
    std::ostream_iterator<decltype(in[0])> out_it (lhs,", ");
    for (size_t i=0; i<in.storage_.size(); ++i) {*out_it = in[i]; out_it++;};
    lhs << "]";
    return lhs;
}


template <typename ValueType, size_t N, typename BoostCType> 
ValueType container_base<ValueType,N, BoostCType>::sum() const 
{
    return EigenMap(storage_.origin(), storage_.num_elements()).sum();
};

template <typename ValueType, size_t N, typename BoostCType> 
int container_base<ValueType,N, BoostCType>::size() const 
{
    return storage_.num_elements();
};

template <typename ValueType, size_t N, typename BoostCType> 
template <typename V2, typename DT>
double container_base<ValueType,N, BoostCType>::diff(const container_base<V2,N,DT>& r, bool norm) const
{   
    assert(this->size() == r.size());
    double d = 0;
    for (int i=0; i<this->size(); i++) {
        d+=std::abs(*(this->data() + i) - *(r.data()+i));
        }
    return d/double(norm?this->size():1);
}

template <typename ValueType, size_t N, typename BoostCType> 
std::array<size_t,N> container_base<ValueType,N, BoostCType>::shape() const 
{
    std::array<size_t, N> out; for (int i=0; i<N; i++) out[i] = storage_.shape()[i]; return out; 
}

template <typename ValueType, size_t N, typename BoostCType> 
container_base<ValueType,1,boost::multi_array_ref<ValueType,1>> container_base<ValueType,N, BoostCType>::flatten()
{
    size_t m = this->size();
    typedef boost::multi_array_ref<ValueType,1> arr_ref;
    return container_base<ValueType,1,arr_ref>(arr_ref(storage_.origin(),std::array<size_t,1>({{m}})));
};

template <typename ValueType, size_t N, typename BoostCType> 
template<size_t N2, typename U>
    typename container_base<ValueType,N, BoostCType>::MatrixType 
    container_base<ValueType,N, BoostCType>::as_matrix() const
{
    Eigen::Map<MatrixType> Map1 (storage_.origin(), storage_.shape()[0],storage_.shape()[1]);
    return Map1;
}

template <typename ValueType, size_t N> 
template <size_t N2, typename>
container<ValueType,N>::container(MatrixType rhs):
    container<ValueType,N>(std::array<size_t,2>({{static_cast<size_t>(rhs.rows()), static_cast<size_t>(rhs.cols()) }}))
{
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), storage_.origin());
}

template <typename ValueType, size_t N, typename BoostCType> 
template<size_t N2, typename U>
container_base<ValueType,N, BoostCType>&
    container_base<ValueType,N, BoostCType>::operator=(MatrixType rhs)
{
    assert(rhs.rows() == storage_.shape()[0] && rhs.cols() == storage_.shape()[1]);
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), storage_.origin());
    return *this;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType,N, BoostCType>::MatrixType 
    container_base<ValueType,N,BoostCType>::as_diagonal_matrix() const
{
    size_t size1 = storage_.size();
    Eigen::DiagonalMatrix<ValueType, Eigen::Dynamic> out(size1);
    Eigen::Map<const VectorType> v(storage_.origin(), size1);
    out.diagonal() = v;
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType,N, BoostCType>::VectorType container_base<ValueType,N,BoostCType>::as_vector() const
{
    size_t size1 = storage_.size();
    VectorType out(size1);
    Eigen::Map<const VectorType> v(storage_.origin(), size1);
    out = v;
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template<typename T> typename std::enable_if< std::is_same<T, complex_type>::value, container<ValueType,N>>::type 
    container_base<ValueType,N,BoostCType>::conj() const
{
    container <ValueType,N> out(*this);
    std::transform(out.storage_.origin(), out.storage_.origin() + out.size(), out.storage_.origin(), [](ValueType in){return std::conj(in);});
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template<typename T> typename std::enable_if<!std::is_same<T, complex_type>::value, container<ValueType,N>>::type 
    container_base<ValueType,N,BoostCType>::conj() const
{
    container <ValueType,N> out(*this);
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename T> 
    typename std::enable_if<std::is_convertible<T,ValueType>::value,container_base<ValueType,N,BoostCType>&>::type 
    container_base<ValueType,N,BoostCType>::operator=(T rhs)
{
    std::fill(storage_.origin(), storage_.origin()+storage_.num_elements(), rhs);
    return *this;
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::const_iterator container_base<ValueType, N, BoostCType>::begin() const 
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(storage_.begin(),f1);
};

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::iterator container_base<ValueType, N, BoostCType>::begin()  
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(storage_.begin(),f1);
};

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::const_iterator container_base<ValueType, N, BoostCType>::end() const 
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(storage_.end(),f1);
}

template <typename ValueType, size_t N, typename BoostCType> 
typename container_base<ValueType, N, BoostCType>::iterator container_base<ValueType, N, BoostCType>::end()  
{
    auto f1 = [this](boost_under_type in){return under_type(in);}; 
    return boost::make_transform_iterator(storage_.end(),f1);
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename R> 
typename container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator+=(const R &rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    typename R::EigenMap map2(rhs.boost_container_().origin(),rhs.boost_container_().num_elements());
    map1+=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2> 
typename container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator+=(const R2& rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
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
typename container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator+(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out+=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator-=(const R &rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    typename R::EigenMap map2(rhs.boost_container_().origin(),rhs.boost_container_().num_elements());
    map1-=map2;
    return (*this);
}
 
template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator-=(const R2& rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    map1-=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator-(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
typename container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator-(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out-=rhs; 
    return out;
}



template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator*=(const R &rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    typename R::EigenMap map2(rhs.boost_container_().origin(),rhs.boost_container_().num_elements());

    map1*=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator*=(const R2& rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    map1*=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator*(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2> 
typename container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator*(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out*=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfContainer<R> container_base<ValueType,N,BoostCType>::operator/=(const R &rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    typename R::EigenMap map2(rhs.storage_.origin(),rhs.storage_.num_elements());
    map1/=map2;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
typename container_base<ValueType,N,BoostCType>::template BaseRefIfValue<R2> container_base<ValueType,N,BoostCType>::operator/=(const R2& rhs)
{
    EigenMap map1(storage_.origin(),storage_.num_elements());
    map1/=rhs;
    return (*this);
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R>
typename container_base<ValueType,N,BoostCType>::template ContainerIfContainer<R> container_base<ValueType,N,BoostCType>::operator/(const R &rhs) const
{
    container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

template <typename ValueType, size_t N, typename BoostCType> 
template <typename R2>
typename container_base<ValueType,N,BoostCType>::template ContainerIfValue<R2> container_base<ValueType,N,BoostCType>::operator/(const R2& rhs) const
{
    container<ValueType,N> out(*this); 
    out/=rhs; 
    return out;
}

} // end of namespace GFTools
