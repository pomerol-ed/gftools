#ifndef ___GFTOOLS_KMESH_HPP___
#define ___GFTOOLS_KMESH_HPP___

#include "Grid.hpp"
#include <numeric>

namespace gftools { 

class KMesh : public Grid<real_type, KMesh>
{
    mutable real_type _domain_len = 2.0*PI;
public:
    int _points;
    using Grid<real_type, KMesh>::vals_;
    KMesh(size_t n_points, real_type len = 2.0*PI);
    KMesh(const KMesh& rhs):Grid<real_type, KMesh>(rhs),_domain_len(rhs._domain_len),_points(rhs._points){}
    KMesh(KMesh &&rhs):Grid<real_type, KMesh>(rhs),_domain_len(rhs._domain_len),_points(rhs._points){}
    KMesh():_points(0){};
    KMesh& operator=(KMesh &&rhs) {_points = rhs._points; _domain_len = rhs._domain_len; vals_.swap(rhs.vals_); return (*this);};
    KMesh& operator=(const KMesh &rhs) {_points = rhs._points; _domain_len = rhs._domain_len;vals_ = rhs.vals_; return (*this);};
    std::tuple <bool, size_t, real_type> find (real_type in) const ;
    template <class Obj> auto integrate(const Obj &in) const ->decltype(in(vals_[0]));
    //template <class Obj> auto gridIntegrate(std::vector<Obj> &in) const -> Obj;
    template <class Obj> auto getValue(Obj &in, real_type x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, point x) const ->decltype(in[0]);
    point shift(point in, point shift_arg) const;
    template <class ArgType> point shift(point in, ArgType shift_arg) const;
    template <class ArgType> real_type shift(real_type in, ArgType shift_arg) const;
};

struct KMeshPatch : public KMesh 
{
    std::map<size_t,size_t> mapvals_;
public:
    const KMesh& _parent;
    size_t _npoints;
    using KMesh::vals_;
    KMeshPatch(const KMesh& parent, std::vector<size_t> indices);
    KMeshPatch(const KMesh& parent);
    template <class Obj> auto getValue(Obj &in, real_type x) const ->decltype(in[0]);
    template <class Obj> auto getValue(Obj &in, KMesh::point x) const ->decltype(in[0]);
    size_t getIndex(KMesh::point x) const;
};
template <>
inline std::ostream& operator<<(std::ostream& lhs, const num_io< typename KMesh::point> &in){lhs << std::setprecision(in._prec) << real_type(in._v.val_); return lhs;};
template <>
inline std::istream& operator>>(std::istream& lhs, num_io< typename KMesh::point> &in){real_type v; lhs >> v; in._v.val_ = v; return lhs;};

//
// KMesh
//

inline KMesh::KMesh(size_t n_points, real_type len):
Grid<real_type,KMesh>(0,n_points,[n_points,len](size_t in){return len/n_points*in;}),
_domain_len(len),
_points(n_points)
{
}
/*
KMesh::KMesh(const KMesh& rhs):Grid(rhs.vals_),_points(rhs._points)
{
}

KMesh::KMesh(KMesh &&rhs):Grid(rhs.vals_),_points(rhs._points)
{
}
*/
inline std::tuple <bool, size_t, real_type> KMesh::find (real_type in) const
{
    assert(in>=0 && in < _domain_len);
    int n = std::lround(in/_domain_len*_points);
    if (n<0) { ERROR("KMesh point is out of bounds, " << in << "<" << 0); return std::make_tuple(0,0,0); };
    if (n==_points) n=0; 
    if (n>_points) { ERROR("KMesh point is out of bounds, " << in << ">" << _domain_len); return std::make_tuple(0,_points,0); };
    real_type weight=in/_domain_len*_points-real_type(n);
    return std::make_tuple (1,n,weight);
}

template <class Obj>
inline auto KMesh::getValue(Obj &in, real_type x) const ->decltype(in[0]) 
{
    const auto find_result=this->find(x);
    if (!std::get<0>(find_result)) throw (ex_wrong_index()); 
    return in[std::get<1>(find_result)];
}


template <class Obj>
inline auto KMesh::getValue(Obj &in, KMesh::point x) const ->decltype(in[0]) 
{
    if (x.index_ < vals_.size() && x == vals_[x.index_])
    return in[x.index_];
    else { 
        #ifndef NDEBUG
        ERROR ("Point not found"); 
        #endif
        return this->getValue(in, real_type(x)); 
         };
}

template <class Obj> 
auto KMesh::integrate(const Obj &in) const -> decltype(in(vals_[0]))
{
    decltype(in(vals_[0])) R = in(real_type(vals_[0]));
    R=std::accumulate(vals_.begin()+1, vals_.end(), R,[&](decltype(in(vals_[0]))& y,decltype(vals_[0]) & x) {return y+in(x);}); 
    return R/_points;
}

template <class ArgType>
inline real_type KMesh::shift(real_type in, ArgType shift_arg) const
{
    assert (in>=0 && in < _domain_len);
    //if (std::abs(real_type(shift_arg))<std::numeric_limits<real_type>::epsilon()) return in;
    real_type out;
    out = in + real_type(shift_arg); 
    out-= std::floor(out/_domain_len)*_domain_len;
    return out;
}


template <class ArgType>
inline typename KMesh::point KMesh::shift(point in, ArgType shift_arg) const
{
    if (std::abs(real_type(shift_arg))<std::numeric_limits<real_type>::epsilon()) return in;
    
    double val = this->shift(in.val_, shift_arg);
    auto find_result = this->find(val);
    if (!std::get<0>(find_result)) throw (ex_wrong_index());
    size_t index = std::get<1>(find_result);
    return (*this)[index];
}

inline typename KMesh::point KMesh::shift(point in, point shift_arg) const
{
    size_t index = (in.index_ + shift_arg.index_)%_points;
    #ifndef NDEBUG
    real_type val = this->shift(in.val_, shift_arg.val_);
    if (std::abs(val - vals_[index])>1e-3) throw (ex_wrong_index()); 
    #endif
    //out.val_ = vals_[out.index_];
    return vals_[index];
}



//
// KMeshPatch
//


inline KMeshPatch::KMeshPatch(const KMesh& parent, std::vector<size_t> indices):
    KMesh(indices.size()),
    _parent(parent),
    _npoints(indices.size())
{
    for (size_t i=0; i<_npoints; ++i) {
        vals_[i]=_parent[indices[i]]; 
        mapvals_[size_t(vals_[i])] = i;
        }
}

inline KMeshPatch::KMeshPatch(const KMesh& parent):
    _parent(parent),
    _npoints(parent.getSize())
{
    vals_ = parent.getPoints();
     for (size_t i=0; i<_npoints; ++i) {
        mapvals_[size_t(vals_[i])] = i;
        }
}

template <class Obj> 
inline auto KMeshPatch::getValue(Obj &in, real_type x) const ->decltype(in[0])
{
    const auto find_result=_parent.find(x);
    if (!std::get<0>(find_result)) throw (ex_wrong_index()); 
    return getValue(in, KMesh::point(std::get<1>(find_result), x));
}

template <class Obj> 
inline auto KMeshPatch::getValue(Obj &in, KMesh::point x) const ->decltype(in[0])
{
    return in[getIndex(x)];
}

inline size_t KMeshPatch::getIndex(KMesh::point x) const
{
    auto f1 = mapvals_.find(size_t(x));
    if (f1!=mapvals_.end()) { return f1->second; }
    else throw (ex_wrong_index());
}

} // end of namespace gftools
#endif // endif :: ifndef ___GFTOOLS_KMESH_HPP___
