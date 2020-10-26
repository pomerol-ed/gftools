#ifndef ___GFTOOLS_ENUMERATEGRID_HPP___
#define ___GFTOOLS_ENUMERATEGRID_HPP___

#include <boost/functional/hash.hpp>
#include "grid_base.hpp"

namespace gftools { 

/// A wrapper around int to avoid weird gcc bug
struct int_wrap_enumerate_grid
{
    int v_; 
    operator int() const{return v_;}; 
    explicit operator double() { return double(v_);};
    operator int&(){return v_;}; 
    int_wrap_enumerate_grid(int i=0):v_(i) {}; 
};
/// Make the wrapper hashable
inline std::size_t hash_value(int_wrap_enumerate_grid const& b) 
{ boost::hash<int> h; return h(static_cast<int>(b)); }

/** A grid of real values. */
class enum_grid : public grid_base<int_wrap_enumerate_grid, enum_grid>
{
public:
    typedef grid_base<int_wrap_enumerate_grid, enum_grid> base; 
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(vals_[0]));
    template <class Obj, typename ...OtherArgTypes> auto integrate(const Obj &in, OtherArgTypes... Args) const -> decltype(in(vals_[0],Args...));
    /** Generates a uniform grid.
     * \param[in] min Minimal point
     * \param[in] max Maximal point
     * \param[in] npoints Number of points
     * \param[in] include_last True, if the max point needs to be included
     */
    enum_grid(int min, int max, bool include_last = false);
    enum_grid(const std::vector<int>& in):base( ([&in](){std::vector<typename point::value_type> out(in.size()); for (size_t i=0; i<out.size(); ++i) out[i]=in[i]; return out; })()) {}
    enum_grid(const std::vector<int_wrap_enumerate_grid>& in):base(in){}
    enum_grid(enum_grid&& rhs);
    enum_grid(const enum_grid& rhs):base(rhs){}
    enum_grid& operator=(enum_grid &&rhs);
    std::tuple <bool, size_t, int> find (int in) const ;
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) -> Obj;
    template <class Obj> auto eval(Obj &in, enum_grid::point x) const -> decltype(in[0]);
    template <class Obj> auto eval(Obj &in, int x) const -> decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0);
    //template <class Obj> auto eval(Obj &in, EnumerateGrid::point x) const ->decltype(in[0]);
};





//
// enum_grid implementation
//

inline enum_grid::enum_grid(int min, int max, bool include_last):
grid_base<int_wrap_enumerate_grid, enum_grid>(min,max+include_last,[](int n){return n;})
{
}

inline enum_grid::enum_grid(enum_grid&& rhs):
    base(std::forward<base>(rhs))
{
};

inline enum_grid& enum_grid::operator=(enum_grid &&rhs) 
{
    base(*this) = base(rhs); 
    return (*this);
}; 

inline std::tuple<bool, size_t, int> enum_grid::find (int in) const
{
    if (in<vals_[0].value()) { ERROR("out of bounds"); return std::make_tuple(0,0,0);};
    if (in > vals_[vals_.size()-1].value()) { ERROR("out of bounds"); return std::make_tuple(0,vals_.size(),0);}; 
    return std::make_tuple (1,in-vals_[0].value(),1);
}

template <class Obj>
inline auto enum_grid::eval(Obj &in, int x) const -> 
    decltype(std::declval<typename std::remove_reference<decltype(in[0])>::type>()*1.0) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) { throw (ex_not_found(x,*this)); } 
    return in[std::get<1>(find_result)];
}


template <class Obj>
inline auto enum_grid::eval(Obj &in, enum_grid::point x) const ->decltype(in[0]) 
{
    if (x.index() < vals_.size() && x == vals_[x.index()])
    return in[x.index()];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->eval(in, int(x)); 
         };
}

} // end of namespace gftools
#endif // endif :: #ifndef ___GFTOOLS_REALGRID_HPP___
