#pragma once 

#include "grid_base.hpp"
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
class matsubara_grid : public grid_base<complex_type, matsubara_grid<Fermion>>
{
public:
    typedef grid_base<complex_type, matsubara_grid<Fermion>> base;
    typedef typename base::point point;
    using base::vals_;
    using typename base::ex_wrong_index;
    //typedef typename grid_base<complex_type, matsubara_grid<Fermion>>::point point;
    //using typename grid_base<complex_type, matsubara_grid<Fermion>>::point;
    /** Inverse temperature. */
    const real_type beta_;
    /** Spacing between values. */
    const real_type spacing_;
    /** Min and max numbers of freq. - useful for searching. */
    const int w_min_, w_max_;
    matsubara_grid(int min, int max, real_type beta);
    matsubara_grid(const matsubara_grid &rhs);
    matsubara_grid(matsubara_grid&& rhs);
    matsubara_grid(std::vector<complex_type> const& in);
    int getNumber(complex_type in) const;
    double beta() const { return beta_; }
    int min_n() const { return w_min_; }
    int max_n() const { return w_max_; }

    point find_nearest(complex_type in) const;
    template <class Obj> auto integrate(const Obj &in) const -> decltype(in(vals_[0]));
    template <class Obj> auto prod(const Obj &in) const -> decltype(in(vals_[0]));
    //template <class Obj> auto gridIntegrate(const std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto eval(Obj &in, complex_type x) const ->decltype(in[0]);
    //template <class Obj> auto eval(Obj &in, point x) const ->decltype(in[0]);
    //using base::eval;
};

typedef matsubara_grid<1> fmatsubara_grid;
typedef matsubara_grid<0> bmatsubara_grid;

//
// matsubara_grid implementations
//

template <bool F>
matsubara_grid<F>::matsubara_grid(int min, int max, real_type beta):
    //grid_base<complex_type, matsubara_grid<F>>(min,max,std::bind(Matsubara<F>, std::placeholders::_1, beta)),
    grid_base<complex_type, matsubara_grid<F>>(min,max,[beta](int n){return Matsubara<F>(n,beta);}),
    beta_(beta), 
    spacing_(PI/beta), 
    w_min_(min),
    w_max_(max)
{
}

template <bool F>
matsubara_grid<F>::matsubara_grid(const matsubara_grid<F> &rhs) : 
    grid_base<complex_type, matsubara_grid<F>>(rhs.vals_),
    beta_(rhs.beta_), 
    spacing_(rhs.spacing_), 
    w_min_(rhs.w_min_), 
    w_max_(rhs.w_max_)
{
}

template <bool F>
matsubara_grid<F>::matsubara_grid(matsubara_grid<F>&& rhs):
    grid_base<complex_type, matsubara_grid>(rhs.vals_), 
    beta_(rhs.beta_), 
    spacing_(rhs.spacing_), 
    w_min_(rhs.w_min_), 
    w_max_(rhs.w_max_)
{
}

template <bool F>
matsubara_grid<F>::matsubara_grid(std::vector<complex_type> const& in):
    base(in),
    beta_( PI / std::abs((in[1] - in[0])/2.) ),
    spacing_(std::abs((in[1] - in[0])/2.)),
    w_min_(MatsubaraIndex<F>(in[0], beta_)),
    w_max_(MatsubaraIndex<F>(in[in.size() - 1], beta_))
{
}

template <bool F>
template <class Obj> 
auto matsubara_grid<F>::integrate(const Obj &in) const -> decltype(in(vals_[0]))
{
    decltype(in(this->vals_[0])) R = in(this->vals_[0]);
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {return y+in(x);}); 
    return R/beta_;
}

template <bool F>
template <class Obj> 
auto matsubara_grid<F>::prod(const Obj &in) const -> decltype(in(vals_[0]))
{
    // fix prod for more numerical stability
    decltype(in(vals_[0])) R = in(vals_[0]);
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {return y*in(x);}); 
    //decltype(in(vals_[0])) R = in(vals_[vals_.size()/2]);
    //R=std::accumulate(vals_.begin()+1+vals_.size()/2, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    //R=std::accumulate(vals_.begin(), vals_.begin()+vals_.size()/2, R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) &x) {DEBUG(x << "|" << in(x) << "|" << y << "->" << y*in(x)); return y*in(x);}); 
    return R;
}

template <bool F>
inline typename matsubara_grid<F>::point matsubara_grid<F>::find_nearest (complex_type in) const
{
    int n=getNumber(in);
    #ifndef NDEBUG
    //DEBUG("Invoking matsubara find");
    #endif
    if (n>=w_min_ && n<w_max_) { return vals_[n-w_min_]; }
    else { 
        #ifndef NDEBUG
        //ERROR("Couldn't find the point");
        #endif
        if (n<w_min_) return vals_[0];
        else return vals_[vals_.size()-1];
        };
}

template <bool F>
inline int matsubara_grid<F>::getNumber(complex_type in) const
{
    assert (std::abs(real(in))<std::numeric_limits<real_type>::epsilon());
    return std::lround(imag(in)/spacing_-F)/2;
};

template <bool F>
template <class Obj>
inline auto matsubara_grid<F>::eval(Obj &in, complex_type x) const ->decltype(in[0]) 
{
    const auto find_result=this->find_nearest(x);
    if (!tools::is_float_equal<complex_type>(x,find_result.val_)) { throw (ex_wrong_index()); } 
    return in[find_result.index_];
}

} // end of namespace gftools
