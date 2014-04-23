#pragma once

#include <numeric>

#include "grid_base.hpp"
//#include <unsupported/Eigen/Splines>

namespace gftools { 

/** A grid of real values. */
class real_grid : public grid_base<real_type, real_grid>
{
    real_type min_;
    real_type max_;
    bool is_uniform_ = false;
public:
    /** Generates a uniform grid.
     * \param[in] min Minimal point
     * \param[in] max Maximal point
     * \param[in] npoints Number of points
     * \param[in] include_last True, if the max point needs to be included
     */
    real_grid(real_type min, real_type max, size_t npoints, bool include_last = true);
    real_grid(int min, int max, const std::function<real_type(int)> &f, bool include_last = true);
    real_grid(std::vector<real_type>&& in);
    real_grid(const std::vector<real_type>& in);

    real_type max() const { return max_; }
    real_type min() const { return min_; }

    std::tuple <bool, size_t, real_type> find (real_type in) const ;
    template <class Obj> auto evaluate(Obj &in, real_type x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0);

    template <class Obj, typename ...OtherArgTypes> 
        auto integrate(Obj &&in, OtherArgTypes... Args) const -> 
            typename std::remove_reference<typename std::result_of<Obj(value_type,OtherArgTypes...)>::type>::type;

    template <class Obj>// decltype (std::declval<Obj>()[0])>
        auto integrate(Obj &&in) const -> 
            typename std::remove_reference<decltype (std::declval<Obj>()[0])>::type;
        
    bool check_uniform_();


    using grid_base<real_type, real_grid>::vals_;
};

template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io< typename real_grid::point> &in){lhs << std::setprecision(in.prec_) << in.value_.val_; return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<typename real_grid::point> &out){real_type im; lhs >> im; out.value_.val_ = im; return lhs;};

//
// real_grid implementation
//

inline real_grid::real_grid(real_type min, real_type max, size_t n_points, bool include_last):
    grid_base<real_type,real_grid>(0,n_points,[n_points,max,min,include_last](size_t in){return (max-min)/(n_points-include_last)*in+min;}),
    min_(min),
    max_((include_last?max:vals_[n_points-1])),
    is_uniform_(true)
{
}

inline real_grid::real_grid(int min, int max, const std::function<real_type (int)> &f, bool include_last):
    grid_base(min,max+include_last,f),
    min_(f(min)),
    max_(f(max-include_last))
{
    check_uniform_();
}

inline real_grid::real_grid(std::vector<real_type>&& in):
grid_base<real_type, real_grid>(in)
{
    auto in2(in);
    std::sort(in2.begin(), in2.end());
    size_t npts = in2.size();
    for (int i=0; i<npts; ++i) vals_[i]=point(in2[i],i);
    min_ = in2[0]; max_ = in2[npts-1];
    check_uniform_();
}


inline real_grid::real_grid(const std::vector<real_type>& in)
{
    vals_.reserve(in.size());
    auto in2(in);
    std::sort(in2.begin(), in2.end());
    size_t npts = in2.size();
    for (int i=0; i<npts; ++i) vals_.emplace_back(in2[i],i);
    min_ = in2[0]; max_ = in2[npts-1];
    check_uniform_();
}

template <class Obj, typename ...OtherArgTypes> 
auto real_grid::integrate(Obj &&f, OtherArgTypes... Args) const -> 
    typename std::remove_reference<typename std::result_of<Obj(value_type,OtherArgTypes...)>::type>::type
{
    typedef typename std::result_of<Obj(value_type,OtherArgTypes...)>::type R;
    std::function<R(value_type)> tmp = [&](value_type x){return f(x,Args...);};
    return this->integrate(gftools::extra::function_proxy<Obj,real_grid>(f,*this));
}

bool real_grid::check_uniform_()  
{ 
    bool is_uniform = true; 
    for (int i=0; i<vals_.size()-2 && is_uniform; i++) { 
        is_uniform = is_uniform && tools::is_float_equal(vals_[i+2].val_-vals_[i+1].val_, vals_[i+1].val_ - vals_[i].val_);
        }
    is_uniform_ = is_uniform;
    return is_uniform;
}

template <class Obj>// decltype (std::declval<Obj>()[0])>
    auto real_grid::integrate(Obj &&v) const -> 
        typename std::remove_reference<decltype (std::declval<Obj>()[0])>::type 
{
    typedef typename std::remove_reference<decltype (std::declval<Obj>()[0])>::type R;
    static_assert(std::is_convertible<R,std::complex<double>>::value,"can't integrate");
    int n = this->size();
    R s = 0.0;
    if (!is_uniform_ || n<8) {
        DEBUG("Using trapezoidal integration");
        for (int i=0; i<n-1; i++) s+=(v[i+1] + v[i])*(vals_[i+1].val_ - vals_[i].val_);
        return 0.5*s;
        }
    DEBUG("Using simpson");
    for (int i=4;i<n-4;i++) { s += 48.*v[i]; }
    auto dx = vals_[1].val_ - vals_[0].val_;
    return dx/48.*(17.*v[0] + 59.*v[1] + 43.*v[2]+49.*v[3] + s + 49. *v[n-4] + 43.*v[n-3] + 59.*v[n-2] + 17.*v[n-1]);
}


inline std::tuple <bool, size_t, real_type> real_grid::find (real_type in) const
{
    #ifndef NDEBUG
    DEBUG("Invoking find");
    #endif
    if (in<min_) { ERROR("Point to find is out of bounds, " << in << "<" << min_ ); return std::make_tuple(0,0,0); };
    if (in>max_) { ERROR("Point to find is out of bounds, " << in << ">" << max_ ); return std::make_tuple(0,vals_.size(),0); };
    auto out = std::lower_bound (vals_.begin(), vals_.end(), in);
    size_t i = size_t(out-vals_.begin());
    i--;
    if (i==vals_.size()-1) return std::make_tuple(1,i,1.0);
    real_type val_i = vals_[i];
    real_type weight=(in-val_i)/(vals_[i+1] - val_i);
    return std::make_tuple (1,i,weight);
}


template <class Obj>
inline auto real_grid::evaluate(Obj &in, real_type x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0)
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (ex_wrong_index()); 
// linear spline
    auto prev_index = std::get<1>(find_result);
    auto prevval_ue = in[prev_index];
    auto weight = std::get<2>(find_result);
    auto nextval_ue = in[prev_index+1];
    auto out = prevval_ue + (nextval_ue - prevval_ue)*weight;
    return out;
/*
    VectorType<typename Eigen::Spline2d::PointType> a1(5);
    DEBUG(a1.cols());
    DEBUG(a1[0]);
    exit(0);
    //auto SplineF = Eigen::SplineFitting<Eigen::Spline<real_type,2>>::Interpolate(a1,Eigen::Dynamic);
    //return SplineF(x)[0];
*/
}

} // end of namespace gftools
