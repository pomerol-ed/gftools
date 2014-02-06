#ifndef ___GFTOOLS_ENUMERATEGRID_HPP___
#define ___GFTOOLS_ENUMERATEGRID_HPP___

#include "Grid.hpp"

namespace GFTools { 

/** A grid of real values. */
class EnumerateGrid : public Grid<int,EnumerateGrid>
{
public:
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(_vals[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(_vals[0],Args...));
    /** Generates a uniform grid.
     * \param[in] min Minimal point
     * \param[in] max Maximal point
     * \param[in] npoints Number of points
     * \param[in] include_last True, if the max point needs to be included
     */
    EnumerateGrid(int min, int max, bool include_last = false);
    EnumerateGrid(const std::vector<int>& in);
    std::tuple <bool, size_t, int> find (int in) const ;
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto getValue(Obj &in, EnumerateGrid::point x) const -> decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, int x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0);
    //template <class Obj> auto getValue(Obj &in, EnumerateGrid::point x) const ->decltype(in[0]);
};

/*
template <>
inline std::ostream& operator<<(std::ostream& lhs, const __num_format< typename EnumerateGrid::point> &in){lhs << int(in._v._val); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, __num_format<typename EnumerateGrid::point> &out){int im; lhs >> im; out._v._val = im; return lhs;};
*/

//
// EnumerateGrid implementation
//

inline EnumerateGrid::EnumerateGrid(int min, int max, bool include_last):
    Grid<int, EnumerateGrid>(min,max+include_last,[](int n){return n;})
{}

inline std::tuple<bool, size_t, int> EnumerateGrid::find (int in) const
{
    if (in<_vals[0]._val) { ERROR("out of bounds"); return std::make_tuple(0,0,0);};
    if (in > _vals[_vals.size()-1]._val) { ERROR("out of bounds"); return std::make_tuple(0,_vals.size(),0);}; 
    return std::make_tuple (1,in-_vals[0]._val,1);
}

template <class Obj>
inline auto EnumerateGrid::getValue(Obj &in, int x) const -> 
    decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (exWrongIndex()); } 
    return in[std::get<1>(find_result)];
}


template <class Obj>
inline auto EnumerateGrid::getValue(Obj &in, EnumerateGrid::point x) const ->decltype(in[0]) 
{
    if (x._index < _vals.size() && x == _vals[x._index])
    return in[x._index];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, int(x)); 
         };
}



} // end of namespace GFTools
#endif // endif :: #ifndef ___GFTOOLS_REALGRID_HPP___
