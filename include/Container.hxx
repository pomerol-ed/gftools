#ifndef ___GFTOOLS_CONTAINER_HPP___
#define ___GFTOOLS_CONTAINER_HPP___

#include "Container.h"
#include <fstream>

namespace GFTools {


//
// Container, N!=1
//

template <size_t N, typename ValueType>
template <size_t M>
inline Container<N,ValueType>::Container ( const std::array<size_t, M> &in):
    _vals(std::vector<Container<N-1,ValueType> >(std::get<M-N>(in), Container<N-1,ValueType>(in)))
{
}

template <size_t N, typename ValueType>
inline auto Container<N,ValueType>::operator[](size_t i)->decltype(_vals[0])
{
    if (i>=_vals.size()) throw exWrongIndex();
    return _vals[i];
}

template <size_t N, typename ValueType> 
std::ostream& operator<<(std::ostream& lhs, const Container<N,ValueType> &in)
{
    lhs << "[";
    std::ostream_iterator<Container<N-1,ValueType> > out_it (lhs,", ");
    std::copy(&in._vals.data()[0],&in._vals.data()[in._vals.size()], out_it);
    lhs << "]";
    return lhs;
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::begin()
{
    return _vals.begin();
}

template <size_t N, typename ValueType> 
typename Container<N,ValueType>::iterator Container<N,ValueType>::end()
{
    return _vals.end();
}

template <size_t N, typename ValueType> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
inline Container <N,ValueType> Container<N,ValueType>::conj()
{
    Container <N,ValueType> out(*this);
    for (iterator it1 = this->begin(); it1!=this->end(); it1++) {  
        *it1=it1->conj();
    }
    return out;
}

template <size_t N, typename ValueType> 
inline ValueType Container<N,ValueType>::sum()
{
    ValueType out=0.0;
    out = std::accumulate(_vals.begin(), _vals.end(), out, [](ValueType x, Container<N-1,ValueType> &in){return x+in.sum();});
    return out;
}

template <size_t N, typename ValueType> 
template<typename N2>
inline MatrixType<ValueType> Container<N,ValueType>::getAsMatrix() const
{
    size_t rows = _vals.size();
    size_t cols = _vals[0].getSize();
    MatrixType<ValueType> out(rows, cols);
    for (size_t i=0; i<rows; ++i) { 
        const ValueType *d = _vals[i]._vals.data();
        Eigen::Map<const VectorType<ValueType>> v(d, cols);
        out.row(i) = v;
        }
    return out;
}

template <size_t N, typename ValueType> 
template<typename N2>
inline Container<N,ValueType>::Container(const MatrixType<ValueType> &rhs)
{
    std::array<size_t, 2> shape = {{ static_cast<size_t>(rhs.rows()), static_cast<size_t>(rhs.cols()) }};
    *this = Container<2,ValueType>(shape);
    for (size_t i=0; i<shape[0]; ++i) { 
            std::copy(rhs.row(i).data(), rhs.row(i).data()+shape[1],_vals[i]._vals.data());
        }
}

template <size_t N, typename ValueType> 
template<typename N2>
inline Container<N,ValueType>& Container<N,ValueType>::operator=(MatrixType<ValueType> &&rhs)
{
    assert(rhs.rows() == _vals.size() && rhs.cols() == _vals[0].getSize());
    for (size_t i=0; i<rhs.rows(); ++i) { 
            std::copy(rhs.row(i).data(), rhs.row(i).data()+rhs.cols(),_vals[i]._vals.data());
        }
    return *this;
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

template <typename ValueType> 
inline MatrixType<ValueType> Container<1,ValueType>::getAsDiagonalMatrix() const
{
    size_t size1 = _vals.size();
    Eigen::DiagonalMatrix<ValueType, Eigen::Dynamic> out(size1);
    const ValueType *d = _vals.data();
    Eigen::Map<const VectorType<ValueType>> v(d, size1);
    out.diagonal() = v;
    return out;
}

template <typename ValueType> 
inline VectorType<ValueType> Container<1,ValueType>::getAsVector() const
{
    size_t size1 = _vals.size();
    VectorType<ValueType> out(size1);
    const ValueType *d = _vals.data();
    Eigen::Map<const VectorType<ValueType>> v(d, size1);
    out = v;
    return out;
}


//
// Algebraic operators
//

// Operator=

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator=(const Container<N,ValueType> &rhs)
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

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator=(Container<N,ValueType> &&rhs)
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



template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator=(const ValueType &rhs)
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
template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x+=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator+=(const Container<N,ValueType> &rhs)
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

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator+(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
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
template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x*=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator*=(const Container<N,ValueType> &rhs)
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

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator*(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
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

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator-=(const RhsArg &rhs)
{
    (*this)+=rhs*ValueType(-1);
    return *this;
}

//
// Operator-
//

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator-(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
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

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType>& Container<N,ValueType>::operator/=(const RhsArg &rhs)
{
    std::for_each(_vals.begin(), _vals.end(), [&](Container<N-1,ValueType> &x){x/=rhs;});
    return *this;
}

template <size_t N, typename ValueType> 
inline Container<N,ValueType>& Container<N,ValueType>::operator/=(const Container<N,ValueType> &rhs)
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

template <size_t N, typename ValueType> 
template <typename RhsArg>
inline Container<N,ValueType> Container<N,ValueType>::operator/(const RhsArg &rhs) const
{
    Container<N,ValueType> out(*this);
    out/=rhs;
    return out;
}




} // end of namespace GFTools
#endif
