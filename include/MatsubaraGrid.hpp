#ifndef ___GFTOOLS_MATSUBARA_GRID_HPP___
#define ___GFTOOLS_MATSUBARA_GRID_HPP___

#include "Grid.hpp"
#include <numeric>

namespace gftools { 

template <bool Fermion> inline complex_type Matsubara(int n, real_type beta){return PI*I/beta*complex_type(2*n+Fermion);};
template <bool Fermion> inline int MatsubaraIndex(complex_type in, real_type beta){return std::round((beta*imag(in)/PI-Fermion)/2.0);};

inline complex_type FMatsubara(int n, real_type beta){return Matsubara<1>(n,beta);};
inline complex_type BMatsubara(int n, real_type beta){return Matsubara<0>(n,beta);};
inline int FMatsubaraIndex(complex_type in, real_type beta){return MatsubaraIndex<1>(in,beta);};
inline int BMatsubaraIndex(complex_type in, real_type beta){return MatsubaraIndex<0>(in,beta);};

/** A Grid of fermionic Matsubara frequencies. */
template <bool Fermion>
class MatsubaraGrid : public Grid<complex_type, MatsubaraGrid<Fermion>>
{
public:
    using Grid<complex_type, MatsubaraGrid<Fermion>>::vals_;
    using typename Grid<complex_type, MatsubaraGrid<Fermion>>::ex_wrong_index;
    //typedef typename Grid<complex_type, MatsubaraGrid<Fermion>>::point point;
    using typename Grid<complex_type, MatsubaraGrid<Fermion>>::point;
    /** Inverse temperature. */
    const real_type _beta;
    /** Spacing between values. */
    const real_type _spacing;
    /** Min and max numbers of freq. - useful for searching. */
    const int _w_min, _w_max;
    MatsubaraGrid(int min, int max, real_type beta);
    MatsubaraGrid(const MatsubaraGrid &rhs);
    MatsubaraGrid(MatsubaraGrid&& rhs);
    int getNumber(complex_type in) const;
    std::tuple <bool, size_t, real_type> find (complex_type in) const ;
    template <class Obj> auto integrate(const Obj &in) const -> decltype(in(vals_[0]));
    template <class Obj> auto prod(const Obj &in) const -> decltype(in(vals_[0]));
    //template <class Obj> auto gridIntegrate(const std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, complex_type x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
};

typedef MatsubaraGrid<1> FMatsubaraGrid;
typedef MatsubaraGrid<0> BMatsubaraGrid;

template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io< typename FMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v.val_); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<typename FMatsubaraGrid::point> &out){real_type im; lhs >> im; out._v.val_ = I*im; return lhs;};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io< typename BMatsubaraGrid::point> &in){lhs << std::setprecision(in._prec) << imag(in._v.val_); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<typename BMatsubaraGrid::point> &out){real_type im; lhs >> im; out._v.val_ = I*im; return lhs;};

//
// MatsubaraGrid implementations
//

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(int min, int max, real_type beta):
    //Grid<complex_type, MatsubaraGrid<F>>(min,max,std::bind(Matsubara<F>, std::placeholders::_1, beta)),
    Grid<complex_type, MatsubaraGrid<F>>(min,max,[beta](int n){return Matsubara<F>(n,beta);}),
    _beta(beta), 
    _spacing(PI/beta), 
    _w_min(min),
    _w_max(max)
{
}

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(const MatsubaraGrid<F> &rhs) : 
    Grid<complex_type, MatsubaraGrid<F>>(rhs.vals_),
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
MatsubaraGrid<F>::MatsubaraGrid(MatsubaraGrid<F>&& rhs):
    Grid<complex_type, MatsubaraGrid>(rhs.vals_), 
    _beta(rhs._beta), 
    _spacing(rhs._spacing), 
    _w_min(rhs._w_min), 
    _w_max(rhs._w_max)
{
}

template <bool F>
template <class Obj> 
auto MatsubaraGrid<F>::integrate(const Obj &in) const -> decltype(in(vals_[0]))
{
    decltype(in(this->vals_[0])) R = in(this->vals_[0]);
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {return y+in(x);}); 
    return R/_beta;
}

template <bool F>
template <class Obj> 
auto MatsubaraGrid<F>::prod(const Obj &in) const -> decltype(in(vals_[0]))
{
    decltype(in(vals_[0])) R = in(vals_[0]);
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {return y*in(x);}); 
    //decltype(in(vals_[0])) R = in(vals_[vals_.size()/2]);
    //R=std::accumulate(vals_.begin()+1+vals_.size()/2, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    //R=std::accumulate(vals_.begin(), vals_.begin()+vals_.size()/2, R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    return R;
}


/*
template <class Obj> 
auto MatsubaraGrid::gridIntegrate(const std::vector<Obj> &in) const -> Obj
{
    decltype(in[0]) R = in[0];
    R=std::accumulate(vals_.begin()+1, vals_.end(),R,[&](decltype(in[0])& y, decltype(in[0]) &x) {return y+x;}); 
    return R/_beta;
}
*/

template <bool F>
inline std::tuple <bool, size_t, real_type> MatsubaraGrid<F>::find (complex_type in) const
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
        return std::make_tuple(0,vals_.size(),0); 
        };
    return std::make_tuple (1,n-_w_min,1);
}

template <bool F>
inline int MatsubaraGrid<F>::getNumber(complex_type in) const
{
    assert (std::abs(real(in))<std::numeric_limits<real_type>::epsilon());
    return std::lround(imag(in)/_spacing-F)/2;
};

template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, complex_type x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (ex_wrong_index()); } 
    return in[std::get<1>(find_result)];
}


template <bool F>
template <class Obj>
inline auto MatsubaraGrid<F>::getValue(Obj &in, MatsubaraGrid::point x) const ->decltype(in[0]) 
{
    if (x.index_ < vals_.size() && x == vals_[x.index_])
    return in[x.index_];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, complex_type(x)); 
         };
}

} // end of namespace gftools
#endif // endif # ifndef __GFTOOLS_MATSUBARAGRID_HPP_
