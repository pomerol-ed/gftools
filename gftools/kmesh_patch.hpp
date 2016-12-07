#pragma once

#include <map>
#include <numeric>

#include "kmesh.hpp"

namespace gftools { 

struct kmesh_patch : public kmesh 
{
    std::map<size_t,size_t> mapvals_;
public:
    const kmesh& _parent;
    size_t _npoints;
    using kmesh::vals_;
    using typename kmesh::ex_not_found;
    kmesh_patch(const kmesh& parent, std::vector<size_t> indices);
    kmesh_patch(const kmesh& parent);
    template <class Obj> auto eval(Obj &in, real_type x) const ->decltype(in[0]);
    template <class Obj> auto eval(Obj &in, kmesh::point x) const ->decltype(in[0]);
    size_t get_index(kmesh::point x) const;
};

//
// kmesh_patch
//


inline kmesh_patch::kmesh_patch(const kmesh& parent, std::vector<size_t> indices):
    kmesh(indices.size()),
    _parent(parent),
    _npoints(indices.size())
{
    for (size_t i=0; i<_npoints; ++i) {
        vals_[i]=_parent[indices[i]]; 
        mapvals_[size_t(vals_[i])] = i;
        }
}

inline kmesh_patch::kmesh_patch(const kmesh& parent):
    _parent(parent),
    _npoints(parent.size())
{
    vals_ = parent.points();
     for (size_t i=0; i<_npoints; ++i) {
        mapvals_[size_t(vals_[i])] = i;
        }
}

template <class Obj> 
inline auto kmesh_patch::eval(Obj &in, real_type x) const ->decltype(in[0])
{
    const auto find_result=_parent.find_nearest(x);
    if (!almost_equal(find_result.value(), x)) throw ex_not_found(x, *this); 
    return _parent.eval(in, find_result);
}

template <class Obj> 
inline auto kmesh_patch::eval(Obj &in, kmesh::point x) const ->decltype(in[0])
{
    return in[this->get_index(x)];
}

inline size_t kmesh_patch::get_index(kmesh::point x) const
{
    auto f1 = mapvals_.find(size_t(x));
    if (f1!=mapvals_.end()) { return f1->second; }
    else throw ex_not_found(x,*this); 
}

} // end of namespace gftools
