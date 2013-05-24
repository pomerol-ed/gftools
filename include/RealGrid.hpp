#ifndef ___GFTOOLS_REALGRID_HPP___
#define ___GFTOOLS_REALGRID_HPP___

#include "Grid.hpp"
#include <numeric>
#include <unsupported/Eigen/Splines>

namespace GFTools { 

/** A grid of real values. */
class RealGrid : public Grid<RealType, RealGrid>
{
    RealType _min;
    RealType _max;
public:
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    /** Generates a uniform grid.
     * \param[in] min Minimal point
     * \param[in] max Maximal point
     * \param[in] npoints Number of points
     * \param[in] include_last True, if the max point needs to be included
     */
    RealGrid(RealType min, RealType max, size_t npoints, bool include_last = true);
    RealGrid(int min, int max, const std::function<RealType(int)> &f, bool include_last = true);
    RealGrid(std::vector<RealType>&& in);
    RealGrid(const std::vector<RealType>& in);
    std::tuple <bool, size_t, RealType> find (RealType in) const ;
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto getValue(Obj &in, RealGrid::point x) const -> decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, RealType x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0);
    //template <class Obj> auto getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]);
};


//
// RealGrid implementation
//

inline RealGrid::RealGrid(RealType min, RealType max, size_t n_points, bool include_last):
    Grid<RealType,RealGrid>(0,n_points,[n_points,max,min,include_last](size_t in){return (max-min)/(n_points-include_last)*in+min;}),
    _min(min),
    _max((include_last?max:_vals[n_points-1]))
{
}

inline RealGrid::RealGrid(int min, int max, const std::function<RealType (int)> &f, bool include_last):
    Grid(min,max+include_last,f),
    _min(f(min)),
    _max(f(max-include_last))
{
}

inline RealGrid::RealGrid(std::vector<RealType>&& in)
{
    auto in2(in);
    std::sort(in2.begin(), in2.end());
    size_t npts = in2.size();
    _vals.resize(npts);
    for (int i=0; i<npts; ++i) _vals[i]=point(in2[i],i);
    _min = in2[0]; _max = in2[npts-1];
}


inline RealGrid::RealGrid(const std::vector<RealType>& in)
{
    auto in2(in);
    std::sort(in2.begin(), in2.end());
    size_t npts = in2.size();
    _vals.resize(npts);
    for (int i=0; i<npts; ++i) _vals[i]=point(in2[i],i);
    _min = in2[0]; _max = in2[npts-1];
}

template <class Obj> 
inline auto RealGrid::integrate(const Obj &in) const -> decltype(in(_vals[0]))
{
    decltype(in(_vals[0])) R=0.0;
    for (int i=0; i<_vals.size()-1; ++i) {
        R+=0.5*(in(_vals[i])+in(_vals[i+1]))*(_vals[i+1]-_vals[i]);
        }
    return R;
}

template <class Obj, typename ...OtherArgTypes> 
inline auto RealGrid::integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...))
{
    decltype(in(_vals[0],Args...)) R=0.0;

    for (int i=0; i<_vals.size()-1; ++i) {
        R+=0.5*(in(_vals[i],Args...)+in(_vals[i+1],Args...))*(_vals[i+1]-_vals[i]);
        }
    return R;
}

inline std::tuple <bool, size_t, RealType> RealGrid::find (RealType in) const
{
    #ifndef NDEBUG
    DEBUG("Invoking find");
    #endif
    if (in<_min) { ERROR("Point to find is out of bounds, " << in << "<" << _min ); return std::make_tuple(0,0,0); };
    if (in>_max) { ERROR("Point to find is out of bounds, " << in << ">" << _max ); return std::make_tuple(0,_vals.size(),0); };
    auto out = std::lower_bound (_vals.begin(), _vals.end(), in);
    size_t i = size_t(out-_vals.begin());
    i--;
    if (i==_vals.size()-1) return std::make_tuple(1,i,1.0);
    RealType val_i = _vals[i];
    RealType weight=(in-val_i)/(_vals[i+1] - val_i);
    return std::make_tuple (1,i,weight);
}


template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealType x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0)
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (exWrongIndex()); 
// linear spline
    auto prev_index = std::get<1>(find_result);
    auto prev_value = in[prev_index];
    auto weight = std::get<2>(find_result);
    auto next_value = in[prev_index+1];
    auto out = prev_value + (next_value - prev_value)*weight;
    return out;
/*
    VectorType<typename Eigen::Spline2d::PointType> a1(5);
    DEBUG(a1.cols());
    DEBUG(a1[0]);
    exit(0);
    //auto SplineF = Eigen::SplineFitting<Eigen::Spline<RealType,2>>::Interpolate(a1,Eigen::Dynamic);
    //return SplineF(x)[0];
*/
}

template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]) 
{
    if (checkPoint(x)) return in[x._index];
    else { ERROR ("Point not found"); throw exWrongIndex(); };
}

} // end of namespace GFTools
#endif // endif :: #ifndef ___GFTOOLS_REALGRID_HPP___
