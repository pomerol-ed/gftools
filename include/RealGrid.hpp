#ifndef ___FK_GRID_HPP___
#define ___FK_GRID_HPP___

#include "Grid.hpp"
#include <numeric>

namespace GFTools { 

/** A grid of real values. */
class RealGrid : public Grid<RealType, RealGrid>
{
    RealType _min;
    RealType _max;
public:
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    RealGrid(RealType min, RealType max, size_t npoints);
    RealGrid(int min, int max, const std::function<RealType(int)> &f);
    std::tuple <bool, size_t, RealType> find (RealType in) const ;
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto getValue(Obj &in, RealType x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]);
};


//
// RealGrid implementation
//

inline RealGrid::RealGrid(RealType min, RealType max, size_t n_points):
    Grid<RealType,RealGrid>(0,n_points,[n_points,max,min](size_t in){return (max-min)/n_points*in+min;}),
    _min(min),
    _max(max)
{
}

inline RealGrid::RealGrid(int min, int max, const std::function<RealType (int)> &f):
    Grid(min,max,f),
    _min(f(min)),
    _max(f(max))
{
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
    if (in>=_max) { ERROR("Point to find is out of bounds, " << in << ">" << _max ); return std::make_tuple(0,_vals.size(),0); };
    auto out = std::lower_bound (_vals.begin(), _vals.end(), in);
    size_t i = size_t(out-_vals.begin());
    RealType val_i = _vals[i];
    RealType weight=(in-val_i)/(_vals[i+1]/val_i);
    return std::make_tuple (1,i,weight);
}


template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealType x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (exWrongIndex()); 
    return in[std::get<1>(find_result)];
}

template <class Obj>
inline auto RealGrid::getValue(Obj &in, RealGrid::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { ERROR ("Point not found"); return this->getValue(in, RealType(x)); };
}

} // end of namespace GFTools
#endif
