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
inline MatrixType<ValueType> ContainerBase<ValueType,N, BoostContainerType>::getAsMatrix() const
{
    Eigen::Map<MatrixType<ValueType>> Map1 (_data.origin(), _data.shape()[0],_data.shape()[1]);
    return Map1;
}

template <typename ValueType, size_t N> 
template<typename U>
inline Container<ValueType,N>::Container(MatrixType<ValueType> rhs):
    Container<ValueType,N>(std::array<size_t,2>({{static_cast<size_t>(rhs.rows()), static_cast<size_t>(rhs.cols()) }}))
{
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), _data.origin());
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template<typename N2>
inline ContainerBase<ValueType,N, BoostContainerType>& ContainerBase<ValueType,N, BoostContainerType>::operator=(MatrixType<ValueType> rhs)
{
    assert(rhs.rows() == _data.shape()[0] && rhs.cols() == _data.shape()[1]);
    std::copy(rhs.data(), rhs.data()+rhs.rows()*rhs.cols(), _data.origin());
    return *this;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename U>
inline ContainerBase<ValueType,N, BoostContainerType> ContainerBase<ValueType,N, BoostContainerType>::conj()
{
    ContainerBase <ValueType,N,BoostContainerType> out(*this);
    //for (iterator it1 = this->begin(); it1!=this->end(); it1++) {  
    //    *it1=it1->conj();
    //}
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
inline ContainerBase<ValueType,N, BoostContainerType>& ContainerBase<ValueType,N, BoostContainerType>::operator=(ValueType rhs)
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
inline MatrixType<ValueType> ContainerBase<ValueType,N,BoostContainerType>::getAsDiagonalMatrix() const
{
    size_t size1 = _data.shape()[0];
    Eigen::DiagonalMatrix<ValueType, Eigen::Dynamic> out(size1);
    Eigen::Map<const VectorType<ValueType>> v(_data.origin(), size1);
    out.diagonal() = v;
    return out;
}

template <typename ValueType, size_t N, typename BoostContainerType> 
template <typename U>
inline VectorType<ValueType> ContainerBase<ValueType,N,BoostContainerType>::getAsVector() const
{
    size_t size1 = _data.shape()[0];
    VectorType<ValueType> out(size1);
    Eigen::Map<const VectorType<ValueType>> v(_data.origin(), size1);
    out = v;
    return out;
}





/*
template <typename ValueType, size_t N, typename BoostContainerType> 
typename ContainerBase<ValueType,N, BoostContainerType>::iterator ContainerBase<ValueType,N, BoostContainerType>::begin()
{
    return _data.begin();
}
*/


//
// Container, N!=1
//

/*
template <typename ValueType, size_t N>
inline Container<ValueType,N>::Container ( const std::array<size_t, N> &in):_data(in)
{
}

template <typename ValueType, size_t N>
template<typename U>
inline ContainerRef<ValueType,N-1> Container<ValueType,N>::operator[](size_t i)
{
    if (i>=_data.shape()[0]) throw exWrongIndex();
    ContainerRef<ValueType, N-1> t1(_data[i]);
}

template <typename ValueType, size_t N>
template <typename U>
inline ValueType& Container<ValueType,N>::operator[](size_t i)
{
    if (i>=_data.shape()[0]) throw exWrongIndex();
    return _data[i];
}

template <typename ValueType, size_t N>
template <typename U>
inline ValueType& ContainerRef<ValueType,N>::operator[](size_t i)
{
    if (i>=_data.shape()[0]) throw exWrongIndex();
    return _data[i];
}


template <typename ValueType, size_t N>
ContainerRef<ValueType,N>::ContainerRef(typename boost::multi_array<ValueType,N+1>::reference in):_data(in)
{
}

template <typename ValueType, size_t N, size_t M>
ContainerView<ValueType,N,M>::ContainerView(boost_view_type in):_data(in)
{
}
*/
/*
template <typename ValueType, size_t N> 
typename Container<ValueType,N>::iterator Container<ValueType,N>::begin()
{
    return _vals.begin();
}

template <typename ValueType, size_t N> 
typename Container<ValueType,N>::iterator Container<ValueType,N>::end()
{
    return _vals.end();
}

template <typename ValueType, size_t N> 
inline ValueType Container<ValueType,N>::sum()
{
    ValueType out=0.0;
    out = std::accumulate(_vals.begin(), _vals.end(), out, [](ValueType x, Container<N-1,ValueType> &in){return x+in.sum();});
    return out;
}


//
// Container, N=1
//

template <typename ValueType>
template <size_t M>
inline Container<1,ValueType>::Container ( const std::array<size_t, M> &in):_vals(std::vector<ValueType>(std::get<M-1>(in)))
{
//    DEBUG("Constructing from size_t.");
}

template <typename ValueType>
inline Container<1,ValueType>::Container ( size_t size ):_vals(std::vector<ValueType>(size))
{
//    DEBUG("Constructing from size_t.");
}

template <typename ValueType>
inline ValueType& Container<1,ValueType>::operator[] (size_t i)
{
    return _vals[i];
}

template <typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<1,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<ValueType> out_it (lhs,", ");
    std::copy(in._vals.begin(),in._vals.end(), out_it);
    lhs << "]";
    return lhs;
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::begin()
{
    return _vals.begin();
}

template <typename ValueType> 
typename Container<1,ValueType>::iterator Container<1,ValueType>::end()
{
    return _vals.end();
}


template <typename ValueType> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
inline Container <1,ValueType> Container<1,ValueType>::conj()
{
    Container <1,ValueType> out(*this);
    for_each(out._vals.begin(), out._vals.end(), [&](ComplexType &x){x=std::conj(x);});
    return out;
}

template <typename ValueType> 
inline ValueType Container<1,ValueType>::sum()
{
    ValueType out=0.0;
    out = std::accumulate(_vals.begin(), _vals.end(), out, std::plus<ValueType>());
    return out;
}

template <typename ValueType> 
void Container<1,ValueType>::savetxt(const std::string& fname)
{
    std::ofstream out;
    out.open(fname.c_str());
    for (auto x : _vals)
        {
            out << std::scientific << __num_format<decltype(x)>(x) << std::endl;
        }

    out.close();
}

//
// Algebraic operators
//

// Operator=

template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator=(const Container<ValueType,N> &rhs)
{
    _vals = rhs._vals; 
    return (*this);
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator=(const Container<1,ValueType> &rhs)
{
    _vals = rhs._vals; 
    return (*this);
}

template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator=(Container<ValueType,N> &&rhs)
{
    _vals.swap(rhs._vals); 
    return (*this);
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator=(Container<1,ValueType> &&rhs)
{
    _vals.swap(rhs._vals); 
    return (*this);
}



template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator=(const ValueType &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x=rhs;});
    return (*this);
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator=(const ValueType &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x=rhs;});
    return (*this);
}



// Operator+=
template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N>& Container<ValueType,N>::operator+=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x+=rhs;});
    return *this;
}

template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator+=(const Container<ValueType,N> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    //std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), [](Container<N-1,ValueType>&x, const Container<N-1,ValueType>&y){x+=y;} );
  //  std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), [](Container<N-1,ValueType>&x, const Container<N-1,ValueType>&y){x+=y;} );
    //std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::plus<Container<N-1,ValueType> >());
    for (int i=0; i<_vals.size(); ++i) _vals[i]+=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x+=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator+=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::plus<ValueType>());
    return *this;
}

//Operator +

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N> Container<ValueType,N>::operator+(const RhsArg &rhs) const
{
    Container<ValueType,N> out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out+=rhs;
    return out;
}

//
// Operator*=
//
template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N>& Container<ValueType,N>::operator*=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x*=rhs;});
    return *this;
}

template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator*=(const Container<ValueType,N> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    for (int i=0; i<_vals.size(); ++i) _vals[i]*=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x*=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator*=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::multiplies<ValueType>());
    return *this;
}


//
// Operator*
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N> Container<ValueType,N>::operator*(const RhsArg &rhs) const
{
    Container<ValueType,N> out(*this);
    out*=rhs;
    return out;
}

//
// Operator-=
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*ValueType(-1);
    return *this;
}

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N>& Container<ValueType,N>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*ValueType(-1);
    return *this;
}

//
// Operator-
//

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N> Container<ValueType,N>::operator-(const RhsArg &rhs) const
{
    Container<ValueType,N> out(*this);
    out-=rhs;
    return out;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out-=rhs;
    return out;
}

//
// Operator/=
//

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N>& Container<ValueType,N>::operator/=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x/=rhs;});
    return *this;
}

template <typename ValueType, size_t N> 
inline Container<ValueType,N>& Container<ValueType,N>::operator/=(const Container<ValueType,N> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    for (int i=0; i<_vals.size(); ++i) _vals[i]/=rhs._vals[i];
    return *this;
}

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType>& Container<1,ValueType>::operator/=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](ValueType &x){x/=rhs;});
    return *this;
}

template <typename ValueType> 
inline Container<1,ValueType>& Container<1,ValueType>::operator/=(const Container<1,ValueType> &rhs)
{
    assert(this->_vals.size() == rhs._vals.size());
    std::transform(_vals.begin(), _vals.end(), rhs._vals.begin(), _vals.begin(), std::divides<ValueType>());
    return *this;
}

//
// Operator /
//

template <typename ValueType> 
template <typename RhsArg>
inline Container<1,ValueType> Container<1,ValueType>::operator/(const RhsArg &rhs) const
{
    Container<1,ValueType> out(*this);
    out/=rhs;
    return out;
}

template <typename ValueType, size_t N> 
template <typename RhsArg>
inline Container<ValueType,N> Container<ValueType,N>::operator/(const RhsArg &rhs) const
{
    Container<ValueType,N> out(*this);
    out/=rhs;
    return out;
}

*/


} // end of namespace GFTools
#endif
