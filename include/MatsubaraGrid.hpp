#ifndef ___GFTOOLS_MATSUBARA_GRID_HPP___
#define ___GFTOOLS_MATSUBARA_GRID_HPP___

#include "Grid.hpp"
#include <numeric>

namespace GFTools { 

/** A Grid of fermionic Matsubara frequencies. */
template <bool Fermion>
class MatsubaraGrid : public Grid<ComplexType, MatsubaraGrid<Fermion>>
{
public:
    using Grid<ComplexType, MatsubaraGrid<Fermion>>::_vals;
    using typename Grid<ComplexType, MatsubaraGrid<Fermion>>::exWrongIndex;
    //typedef typename Grid<ComplexType, MatsubaraGrid<Fermion>>::point point;
    using typename Grid<ComplexType, MatsubaraGrid<Fermion>>::point;
    /** Inverse temperature. */
    const RealType _beta;
    /** Spacing between values. */
    const RealType _spacing;
    /** Min and max numbers of freq. - useful for searching. */
    const int _w_min, _w_max;
    MatsubaraGrid(int min, int max, RealType beta);
    MatsubaraGrid(const MatsubaraGrid &rhs);
    MatsubaraGrid(MatsubaraGrid&& rhs);
    int getNumber(ComplexType in) const;
    std::tuple <bool, size_t, RealType> find (ComplexType in) const ;
    template <class Obj> auto integrate(const Obj &in) const -> decltype(in(_vals[0]));
    template <class Obj> auto prod(const Obj &in) const -> decltype(in(_vals[0]));
    //template <class Obj> auto gridIntegrate(const std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, ComplexType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
};

typedef MatsubaraGrid<1> FMatsubaraGrid;
typedef MatsubaraGrid<0> BMatsubaraGrid;

template <>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format< typename FMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v._val); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, __num_format<typename FMatsubaraGrid::point> &out){RealType im; lhs >> im; out._v._val = I*im; return lhs;};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format< typename BMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v._val); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, __num_format<typename BMatsubaraGrid::point> &out){RealType im; lhs >> im; out._v._val = I*im; return lhs;};

//
// MatsubaraGrid implementations
//

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(int min, int max, RealType beta):
    //Grid<ComplexType, MatsubaraGrid<F>>(min,max,std::bind(Matsubara<F>, std::placeholders::_1, beta)),
    Grid<ComplexType, MatsubaraGrid<F>>(min,max,[beta](int n){return Matsubara<F>(n,beta);}),
    _beta(beta), 
    _spacing(PI/beta), 
    _w_min(min),
    _w_max(max)
{
}

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(const MatsubaraGrid<F> &rhs) : 
    Grid<ComplexType, MatsubaraGrid<F>>(rhs._vals),
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(MatsubaraGrid<F>&& rhs):
    Grid<ComplexType, MatsubaraGrid>(rhs._vals), 
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
template <class Obj> 
auto MatsubaraGrid<F>::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(this->_vals[0])) R = in(this->_vals[0]);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {return y+in(x);}); 
    return R/_beta;
}

template <bool F>
template <class Obj> 
auto MatsubaraGrid<F>::prod(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R = in(_vals[0]);
    R=std::accumulate(_vals.begin()+1, _vals.end(), R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {return y*in(x);}); 
    //decltype(in(_vals[0])) R = in(_vals[_vals.size()/2]);
    //R=std::accumulate(_vals.begin()+1+_vals.size()/2, _vals.end(), R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    //R=std::accumulate(_vals.begin(), _vals.begin()+_vals.size()/2, R,[&](decltype(in(_vals[0]))& y,decltype(_vals[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    return R;
}


/*
template <class Obj> 
auto MatsubaraGrid::gridIntegrate(const std::vector<Obj> &in) const -> Obj
{
    decltype(in[0]) R = in[0];
    R=std::accumulate(_vals.begin()+1, _vals.end(),R,[&](decltype(in[0])& y, decltype(in[0]) &x) {return y+x;}); 
    return R/_beta;
}
*/

template <bool F>
inline std::tuple <bool, size_t, RealType> MatsubaraGrid<F>::find (ComplexType in) const
{
    int n=getNumber(in);
    #ifndef NDEBUG
    DEBUG("Invoking matsubara find");
    #endif
    if (n<_w_min) { 
        #ifndef NDEBUG
        ERROR("Frequency to find is out of bounds, " << in << "<" << FMatsubara(_w_min,_beta)); 
        #endif
        return std::make_tuple(0,0,0); 
        };
    if (n>=_w_max) { 
        #ifndef NDEBUG
        ERROR("Frequency to find is out of bounds, " << in << ">" << FMatsubara(_w_max,_beta)); 
        #endif
        return std::make_tuple(0,_vals.size(),0); 
        };
    return std::make_tuple (1,n-_w_min,1);
}

template <bool F>
inline int MatsubaraGrid<F>::getNumber(ComplexType in) const
{
    assert (std::abs(real(in))<std::numeric_limits<RealType>::epsilon());
    return std::lround(imag(in)/_spacing-F)/2;
};

template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, ComplexType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (exWrongIndex()); } 
    return in[std::get<1>(find_result)];
}


template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, MatsubaraGrid::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, ComplexType(x)); 
         };
}

} // end of namespace GFTools
#endif // endif # ifndef __GFTOOLS_MATSUBARAGRID_HPP_
