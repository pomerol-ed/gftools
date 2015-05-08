#pragma once

#include <map>
#include <numeric>

#include "grid_base.hpp"
#include "almost_equal.hpp"

namespace gftools { 

class kmesh : public grid_base<real_type, kmesh>
{
public:
    typedef grid_base<real_type, kmesh> base;
    using grid_base<real_type, kmesh>::vals_;
    kmesh(size_t n_points, real_type len = 2.0*M_PI);
    kmesh(const kmesh& rhs):grid_base<real_type, kmesh>(rhs),domain_len_(rhs.domain_len_),npoints_(rhs.npoints_){}
    kmesh(kmesh &&rhs):base(std::forward<base>(rhs)),domain_len_(rhs.domain_len_),npoints_(rhs.npoints_){}
    kmesh():npoints_(0){};
    kmesh(std::vector<real_type> const& in);
    kmesh& operator=(kmesh &&rhs) {npoints_ = rhs.npoints_; domain_len_ = rhs.domain_len_; vals_.swap(rhs.vals_); return (*this);};
    kmesh& operator=(const kmesh &rhs) {npoints_ = rhs.npoints_; domain_len_ = rhs.domain_len_;vals_ = rhs.vals_; return (*this);};
    std::tuple <bool, size_t, real_type> find (real_type in) const ;
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(vals_[0]));
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto eval(Obj &in, real_type x) const ->decltype(in[0]);
    template <class Obj> auto eval(Obj &in, point x) const ->decltype(in[0]) { return base::eval(in,x); }
    real_type shift(real_type in,real_type shift_arg) const;
    point shift(point in, real_type shift_arg) const { return static_cast<const base*>(this)->shift(in,shift_arg); }
    point shift(point in, point shift_arg) const { return static_cast<const base*>(this)->shift(in,shift_arg); }
protected:
    ///number of equidistantly spaced points in k-space
    int npoints_;
    ///extent of domain in k-space, usually 2*PI
    real_type domain_len_ = 2.0*M_PI;
};

//
// kmesh
//

inline kmesh::kmesh(size_t n_points, real_type len):
grid_base<real_type,kmesh>(0,n_points,[n_points,len](size_t in){return len/n_points*in;}),
domain_len_(len),
npoints_(n_points)
{
}

inline kmesh::kmesh(std::vector<real_type> const& in):
    base(in)
{
    if (vals_.size() <= 1) throw std::logic_error("Can't construct a kmesh with <= 1 element");
    double diff = vals_[1] - vals_[0];
    bool ok = true;
    for (int i=1; i<vals_.size() && ok; i++) { 
         ok = almost_equal(vals_[i] - vals_[i-1], diff, diff * 1e-6); 
        }
    if (!ok) throw std::logic_error("kmesh should be uniform");
    domain_len_ = vals_.size() * diff;
    npoints_ = vals_.size();
}


inline std::tuple <bool, size_t, real_type> kmesh::find (real_type in) const
{
    assert(in>=0 && in < domain_len_);
    int n = std::lround(in/domain_len_*npoints_);
    if (n<0) { ERROR("kmesh point is out of bounds, " << in << "<" << 0); return std::make_tuple(0,0,0); };
    if (n==npoints_) n=0; 
    if (n>npoints_) { ERROR("kmesh point is out of bounds, " << in << ">" << domain_len_); return std::make_tuple(0,npoints_,0); };
    real_type weight=in/domain_len_*npoints_-real_type(n);
    return std::make_tuple (1,n,weight);
}


template <class Obj>
inline auto kmesh::eval(Obj &in, real_type x) const ->decltype(in[0]) 
{
	auto p = find_nearest(x);
    return eval(in,p);
}


template <class Obj> 
auto kmesh::integrate(const Obj &in) const -> decltype(in(vals_[0]))
{
    decltype(in(vals_[0])) R = in(real_type(vals_[0]));
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) & x) {return y+in(x);}); 
    return R/npoints_;
}

inline real_type kmesh::shift(real_type in, real_type shift_arg) const
{
    assert (in>=0 && in < domain_len_);
    real_type out;
    out = in + real_type(shift_arg); 
    out-= domain_len_*(almost_equal(out, domain_len_, num_io<double>::tolerance())?1.0:std::floor(out/domain_len_));
    return out;
}

} // end of namespace gftools
