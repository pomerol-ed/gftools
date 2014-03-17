#ifndef ___GFTOOLS_ENUMERATEGRID_HPP___
#define ___GFTOOLS_ENUMERATEGRID_HPP___

#include "Grid.hpp"

namespace gftools { 

// A wrapper around int to avoid weird gcc buf
struct int_wrap_enumerate_grid
{int v_; operator int() const{return v_;}; int_wrap_enumerate_grid(int i=0):v_(i) {}; };

/** A grid of real values. */
class EnumerateGrid : public Grid<int_wrap_enumerate_grid, EnumerateGrid>
{
public:
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(vals_[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(vals_[0],Args...));
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
    template <class Obj> auto evaluate(Obj &in, EnumerateGrid::point x) const -> decltype(in[0]);
    template <class Obj> auto evaluate(Obj &in, int x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0);
    //template <class Obj> auto evaluate(Obj &in, EnumerateGrid::point x) const ->decltype(in[0]);
};



template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io< typename EnumerateGrid::point> &in){lhs << int(in.value_.val_); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io<typename EnumerateGrid::point> &out){int im; lhs >> im; out.value_.val_ = im; return lhs;};

//
// EnumerateGrid implementation
//

inline EnumerateGrid::EnumerateGrid(int min, int max, bool include_last):
Grid<int_wrap_enumerate_grid, EnumerateGrid>(min,max+include_last,[](int n){return n;})
{
}

inline std::tuple<bool, size_t, int> EnumerateGrid::find (int in) const
{
    if (in<vals_[0].val_) { ERROR("out of bounds"); return std::make_tuple(0,0,0);};
    if (in > vals_[vals_.size()-1].val_) { ERROR("out of bounds"); return std::make_tuple(0,vals_.size(),0);}; 
    return std::make_tuple (1,in-vals_[0].val_,1);
}

template <class Obj>
inline auto EnumerateGrid::evaluate(Obj &in, int x) const -> 
    decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (std::logic_error("Wrong index")); } 
    return in[std::get<1>(find_result)];
}


template <class Obj>
inline auto EnumerateGrid::evaluate(Obj &in, EnumerateGrid::point x) const ->decltype(in[0]) 
{
    if (x.index_ < vals_.size() && x == vals_[x.index_])
    return in[x.index_];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->evaluate(in, int(x)); 
         };
}



} // end of namespace gftools
#endif // endif :: #ifndef ___GFTOOLS_REALGRID_HPP___
