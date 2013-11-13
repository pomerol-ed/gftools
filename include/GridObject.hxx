#ifndef ___GFTOOLS_GRIDOBJECT_HXX___
#define ___GFTOOLS_GRIDOBJECT_HXX___

#include "GridObject.hpp"
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>

namespace GFTools {
//
// Specification of integration with GridObject
//


//
// GridObject::ContainerExtractor
//

// ContainerExtactor::get_ref

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes>
inline ValueType& GridObject<ValueType,GridTypes...>::ContainerExtractor<Nc,CT,ArgType1,ArgTypes...>::get_ref(
    CT &data, const std::tuple<GridTypes...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    static_assert(std::is_same<ArgType1,typename std::tuple_element<N-Nc,std::tuple<GridTypes...>>::type::point>::value,"Argument to a reference operator is not a point");
    if (!grid.checkPoint(arg1)) throw (exPointMismatch());
    auto &tmp = data[arg1._index];
    return ContainerExtractor<Nc-1, decltype(tmp), ArgTypes...>::get_ref(tmp,grids,args...);
};

template< typename ValueType, typename ...GridTypes>
template <typename CT, typename ArgType1> 
inline ValueType& GridObject<ValueType,GridTypes...>::ContainerExtractor<1,CT,ArgType1>::get_ref(
    CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    static_assert(std::is_same<ArgType1,typename std::tuple_element<N-1,std::tuple<GridTypes...>>::type::point>::value,"Argument to a reference operator is not a point");
    if (!grid.checkPoint(arg1)) throw (exPointMismatch());
    auto &tmp = data[arg1._index];
    return tmp;
}

// ContainerExtactor::get

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes>
inline ValueType GridObject<ValueType,GridTypes...>::ContainerExtractor<Nc,CT,ArgType1,ArgTypes...>::get(
    CT &data, const std::tuple<GridTypes...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    auto tmp = grid.getValue(data, arg1);
    return ContainerExtractor<Nc-1, decltype(tmp), ArgTypes...>::get(tmp,grids,args...);
};

template< typename ValueType, typename ...GridTypes>
template <typename CT, typename ArgType1> 
inline ValueType GridObject<ValueType,GridTypes...>::ContainerExtractor<1,CT,ArgType1>::get(
    CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    auto tmp = grid.getValue(data, arg1);
    return tmp;
}

//
// GridObject
//
    
template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( const std::tuple<GridTypes...> &grids, const Container<ValueType, sizeof...(GridTypes)>& data):
    _grids(grids),
    _dims(GetGridSizes<N>::TupleSizeToArray(_grids)),
    _data(data),
    _f(__fun_traits<FunctionType>::constant(0.0)) 
{
};

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( const std::tuple<GridTypes...> &in):
    _grids(in),
    _dims(GetGridSizes<N>::TupleSizeToArray(_grids)),
    _data(_dims),
    _f(__fun_traits<FunctionType>::constant(0.0))
{
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( const GridObject<ValueType, GridTypes...>& rhs):
    _grids(rhs._grids), 
    _dims(rhs._dims),
    _data(rhs._data),
    _f(rhs._f)
{
}; 

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>::GridObject( GridObject<ValueType,GridTypes...> && rhs):
    _grids(rhs._grids),
    _data(rhs._data)
{
    _dims.swap(rhs._dims);
    _f.swap(rhs._f);
}

template <typename ValueType, typename ...GridTypes> 
const std::tuple<GridTypes...> GridObject<ValueType,GridTypes...>::getGrids() const 
{ 
    return _grids; 
};

template <typename ValueType, typename ...GridTypes> 
inline auto GridObject<ValueType,GridTypes...>::operator[](size_t i)->decltype(_data[i])
{
    return _data[i];
}

/*
template <typename ValueType, typename ...GridTypes> 
template <int M> 
ValueType& GridObject<ValueType,GridTypes...>::operator[](const std::array<size_t,M>& in)
{
}*/

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridTypes...>::get(const ArgTypes&... in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    return ContainerExtractor<sizeof...(GridTypes), Container<ValueType, N>, ArgTypes...>::get_ref(_data,_grids,in...);
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& GridObject<ValueType,GridTypes...>::get(const std::tuple<ArgTypes...> &in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    std::function<ValueType&(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType&{return this->get<ArgTypes...>(in1...); };// ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(_data,_grids,in...);};
    __caller<ValueType&,ArgTypes...> t = {in,f1};
    return t.call();
//ContainerExtractor<sizeof...(GridTypes), Container<ValueType, N>, ArgTypes...>::get_ref(_data,_grids,in);
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType& GridObject<ValueType,GridTypes...>::get(const PointTupleType& in)
{
        auto indices = getIndicesFromPoints(in);
        return _data._data(indices);
}



template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const ArgTypes&... in) const
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    try { return ContainerExtractor<sizeof...(GridTypes), Container<ValueType, N>, ArgTypes...>::get(_data,_grids,in...); }
    //catch (std::exception &e) { 
    catch (...) { 
        #ifndef NDEBUG
        DEBUG("Using analytical expression");
        #endif
        return _f(in...);
     };

}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const std::tuple<ArgTypes...>& in) const
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "GridObject call, number of input parameters mismatch."); 
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return this->template operator()<ArgTypes...>(in1...); };// ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::get(_data,_grids,in...);};
    __caller<ValueType,ArgTypes...> t = {in,f1};
    return t.call();
}

/*
template <typename ValueType, typename ...GridTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const ArgTupleType& in) const
{
    try {
    return ContainerExtractor<sizeof...(GridTypes), Container<ValueType, N>, ArgTupleType>::get(_data,_grids,in);
    }
    catch (...) { 
        #ifndef NDEBUG
        DEBUG("Using analytical expression");
        #endif
        __caller_tuple<ValueType,ArgTupleType> t = {in,_f};
        return t.call();
     };
}*/

template <typename ValueType, typename ...GridTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::operator()(const PointTupleType& in) const
{
    try {
        //__tuple_print<PointTupleType>::print(in);
        auto indices = getIndicesFromPoints(in);
        //for (auto x:indices) INFO_NONEWLINE(x); INFO("");
        return _data._data(indices);
    }
    catch (exWrongIndex)
    {
        return this->operator()(ArgTupleType(in));
    };
}



template <typename ValueType, typename ...GridTypes> 
std::ostream& operator<<(std::ostream& lhs, const GridObject<ValueType,GridTypes...> &in)
{
    lhs << (in._data);
    return lhs;
}

template <typename ValueType, typename ...GridTypes> 
template <size_t M>
auto GridObject<ValueType,GridTypes...>::getGrid() const -> const decltype(std::get<M>(_grids))
{
    return std::get<M>(_grids);
}



template <typename ValueType, typename ...GridTypes> 
auto GridObject<ValueType,GridTypes...>::getGrid() const -> const decltype(std::get<0>(_grids))
{
    return std::get<0>(_grids);
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
void GridObject<ValueType,GridTypes...>::fill(const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(FunctionType(in));
}

/*
template <typename ValueType, typename ...GridTypes> 
inline void GridObject<ValueType,GridTypes...>::fill(const typename GridObject<ValueType,GridTypes...>::PointFunctionType& in)
{

    __genContainerExtractor<sizeof...(GridTypes), ContainerExtractor, std::tuple<GridTypes...>>::type::set(_data,_grids,in);

    //ContainerExtractor<sizeof...(GridTypes), ArgTypes...>::set(_data,_grids,in);
    //_f = in;
}
*/

template <typename ValueType, typename ...GridTypes> 
inline const size_t GridObject<ValueType,GridTypes...>::getTotalContainerSize() const
{
    size_t out = 1;
    for (auto i : _dims) out*=i;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename GridObject<ValueType,GridTypes...>::PointIndices GridObject<ValueType,GridTypes...>::_getPointsIndices(const size_t index) const
{
    PointIndices indices;
    size_t t = index;
    for (int i=N-1; i>=0; i--) { 
        indices[i]=t%_dims[i];
        t-=indices[i];
        t/=_dims[i];
        }
    return indices;
}


template <typename ValueType, typename ...GridTypes>
inline typename GridObject<ValueType,GridTypes...>::ArgTupleType GridObject<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    return this->getArgsFromIndices<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename GridObject<ValueType,GridTypes...>::ArgTupleType GridObject<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    ArgTupleType out;
    auto t1 = std::get<N-1>(_grids)[in[N-1]]._val;
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename GridObject<ValueType,GridTypes...>::ArgTupleType GridObject<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    auto out = getArgsFromIndices<M-1>(in);
    auto t1 = std::get<N-1-M>(_grids)[in[N-1-M]];
    std::get<N-1-M>(out) = t1._val;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename GridObject<ValueType,GridTypes...>::PointTupleType GridObject<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    return this->getPointsFromIndices<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename GridObject<ValueType,GridTypes...>::PointTupleType GridObject<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    PointTupleType out;
    auto t1 = std::get<N-1>(_grids)[in[N-1]];
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename GridObject<ValueType,GridTypes...>::PointTupleType GridObject<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    auto out = getPointsFromIndices<M-1>(in);
    auto t1 = std::get<N-1-M>(_grids)[in[N-1-M]];
    std::get<N-1-M>(out) = t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename GridObject<ValueType,GridTypes...>::PointIndices GridObject<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    return this->getIndicesFromPoints<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename GridObject<ValueType,GridTypes...>::PointIndices GridObject<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    PointIndices out;
    auto t1 = std::get<N-1>(in);
    if (std::get<N-1>(_grids).getSize()<=t1._index) throw exWrongIndex();
    if (std::get<N-1>(_grids)[t1._index]._index != t1._index) { throw exWrongIndex(); };
    out[N-1]=t1._index;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename GridObject<ValueType,GridTypes...>::PointIndices GridObject<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    auto out = getIndicesFromPoints<M-1>(in);
    auto t1 = std::get<N-1-M>(in);
    if (std::get<N-1-M>(_grids).getSize()<=t1._index) throw exWrongIndex();
    if (std::get<N-1-M>(_grids)[t1._index]._index != t1._index) throw exWrongIndex();
    out[N-1-M]=t1._index;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::fill_tuple(const std::function<ValueType(ArgTupleType)>& in)
{
    size_t total_size = this->getTotalContainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        ArgTupleType args = this->getArgsFromIndices(pts_index);
        auto val = in(args);
        _data._data(pts_index) = val;
        };
    _f = __fun_traits<FunctionType>::getFromTupleF(in); 

}

template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::fill(const typename GridObject<ValueType,GridTypes...>::FunctionType& in)
{
    size_t total_size = this->getTotalContainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        ArgTupleType args = this->getArgsFromIndices(pts_index);
        typename GridArgTypeExtractor<ValueType, std::tuple<GridTypes...> >::arg_function_wrapper t = {args, in};
        auto val = t.call();
        _data._data(pts_index) = val;
        };
    _f = in;
}

template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::fill(const typename GridObject<ValueType,GridTypes...>::PointFunctionType& in)
{
    size_t total_size = this->getTotalContainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        PointTupleType args = this->getPointsFromIndices(pts_index);
        typename GridPointExtractor<ValueType, std::tuple<GridTypes...> >::point_function_wrapper t = {args, in};
        auto val = t.call();
        _data._data(pts_index) = val;
        };
}

template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::fill_tuple(const std::function<ValueType(PointTupleType)>& in)
{
    size_t total_size = this->getTotalContainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        PointTupleType args = this->getPointsFromIndices(pts_index);
        auto val = in(args);
        _data._data(pts_index) = val;
        };
    //_f = __fun_traits<FunctionType>::getFromTupleF(std::function<ValueType(ArgTupleType)>(in)); 
}




template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
RealType GridObject<ValueType,GridTypes...>::diff(const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject outObj(_grids);
    auto f1 = [&](PointTupleType in){return std::abs((*this)(in) - rhs(in));};
    PointFunctionType f = __fun_traits<PointFunctionType>::getFromTupleF(f1); 
    outObj.fill(f);
    RealType norm = 1.0;
    for (auto v : _dims) { norm*=v; };
    return std::real(outObj.sum())/norm;
}

template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, RealType>::value, int>::type>
RealType GridObject<ValueType,GridTypes...>::diff(const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject outObj(_grids);
    auto f1 = [&](PointTupleType in){return std::abs((*this)(in) - rhs(in));};
    PointFunctionType f = __fun_traits<PointFunctionType>::getFromTupleF(f1); 
    outObj.fill(f);
    RealType norm = 1.0;
    for (auto v : _dims) { norm*=v; };
    return std::real(outObj.sum())/norm;
}


template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type>
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::conj()
{
    GridObject<ValueType,GridTypes...> out(*this);
    (out._data) = out._data.conj();
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::sum()
{
    return _data.sum();
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType GridObject<ValueType,GridTypes...>::__get_f(const std::tuple<ArgTypes...>& in) const
{
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return this->_f(in1...); };
    __caller<ValueType,ArgTypes...> t = {in,f1};
    return t.call();
}



template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::shift(ArgTypes... args) const
{
    return this->shift(std::forward_as_tuple(args...));
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::shift(const std::tuple<ArgTypes...>& shift_args) const
{
    GridObject<ValueType,GridTypes...> out(_grids);
    std::function<ValueType(PointTupleType)> ShiftFunction = [&](PointTupleType args1)->ValueType { 
        PointTupleType out_args = this->_shiftArgs(args1, shift_args);
    //    __tuple_print<PointTupleType>::print(args1); 
    //    INFO_NONEWLINE("+");  __tuple_print<std::tuple<ArgTypes...>>::print(shift_args); 
    //    INFO_NONEWLINE("-->");__tuple_print<PointTupleType>::print(out_args);
        return (*this)(out_args);
        };
    PointFunctionType fillF = __fun_traits<PointFunctionType>::getFromTupleF(ShiftFunction);
    //out.fill(fillF);
    DEBUG("!");
    out.fill_tuple(ShiftFunction);
    DEBUG("!!");
    
    static std::function<ValueType(ArgTupleType)> ShiftAnalyticF;
    ShiftAnalyticF = [this, shift_args](const ArgTupleType& in)->ValueType {
        ArgTupleType out_args = _shiftArgs(in,shift_args); 
        return __get_f(out_args);
    };
    
    FunctionType tailF = __fun_traits<FunctionType>::getFromTupleF(ShiftAnalyticF);
    out._f = tailF;
    
    return out;
}

/*
template <typename ValueType, typename ...GridTypes> 
    template <typename OrigArg1, typename ...OrigArgs, typename ArgType1, typename ...ArgTypes, 
        typename std::enable_if<sizeof...(OrigArgs)==sizeof...(ArgTypes)>,  
        typename std::enable_if<sizeof...(OrigArgs)!=0> > 
std::tuple<OrigArg1, OrigArgs...> GridObject<ValueType,GridTypes...>::_shiftArgs(const std::tuple<OrigArg1, OrigArgs...>&in, const std::tuple<ArgType1, ArgTypes...>& shift_args) const
{
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-sizeof...(ArgTypes)-1>(_grids).shift(arg1,shift_arg1); 


    auto f_o = [&in](OrigArg1 arg1, OrigArgs... others)->std::tuple<OrigArgs...>{ return std::forward_as_tuple(others...);};
    auto f_s = [&shift_args](ArgType1 arg1, ArgTypes... others)->std::tuple<ArgTypes...>{ return std::forward_as_tuple(others...);};

    __caller<std::tuple<OrigArgs...>,OrigArg1,OrigArgs...> t_o = {in,f_o};
    std::tuple<OrigArgs...> other_orig_args(t_o.call());

    __caller<std::tuple<ArgTypes...>,ArgType1,ArgTypes...> t_s = {shift_args,f_s};
    std::tuple<ArgTypes...> other_shift_args(t_s.call());

    return std::tuple_cat(std::forward_as_tuple(out1),this->_shiftArgs(other_orig_args,other_shift_args));
}
*/
template <typename ValueType, typename ...GridTypes> 
template <typename OrigArg1, typename ArgType1> 
inline std::tuple<OrigArg1> GridObject<ValueType,GridTypes...>::_shiftArgs(const std::tuple<OrigArg1>&in, const std::tuple<ArgType1>& shift_args) const
{
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-1>(_grids).shift(arg1,shift_arg1); 
    return std::forward_as_tuple(out1);
}

template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::savetxt(const std::string& fname) const
{
    INFO("Saving " << typeid(*this).name() << " to " << fname);
    std::ofstream out;
    out.open(fname.c_str());
    size_t total_size = this->getTotalContainerSize();
    size_t last_grid_size = std::get<N-1>(_grids).getSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        //ArgTupleType args = this->getArgsFromIndices(pts_index);
        PointTupleType pts = this->getPointsFromIndices(pts_index);
        auto val = (*this)(pts);
        out << std::scientific << __tuple_print<PointTupleType>::serialize(pts) << "    " << __num_format<ValueType>(val) << std::endl;
        if (N > 1 && i && (i+1)%last_grid_size==0) out << std::endl;
        };
    out.close();
}

template <typename ValueType, typename ...GridTypes> 
void GridObject<ValueType,GridTypes...>::loadtxt(const std::string& fname, RealType tol)
{
    INFO("Loading " << typeid(*this).name() << " from " << fname);
    std::ifstream in;
    in.open(fname.c_str());
    if (in.fail()) { ERROR("Couldn't open file " << fname); throw exIOProblem(); };
    size_t total_size = this->getTotalContainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);

        PointTupleType pts = this->getPointsFromIndices(pts_index);
        //ArgTupleType args = this->getArgsFromIndices(pts_index);
        ArgTupleType pts2 = __tuple_print<PointTupleType>::read(in);
        if (!__is_equal<ArgTupleType>(pts,pts2,tol)) throw (exIOProblem());

        __num_format<ValueType> tmp2(this->get(pts));
        in >> tmp2;
        this->get(pts) = tmp2._v;
        };


    in.close();
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::copyInterpolate (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //_data=rhs._data;
    _f = rhs._f;
    const std::function<ValueType(ArgTupleType)> bindf = [&](ArgTupleType in){return rhs(in);};
    this->fill_tuple(bindf);

    return *this;
}


//
// Operators
//

/*
template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(in);
    return *this;
}
*/

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data=rhs._data;
    _f = rhs._f;
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator= (
    const ValueType& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data=rhs;
    _f = __fun_traits<FunctionType>::constant(rhs);
    return *this;
}



template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data+=rhs._data;
    //_f=__fun_traits<FunctionType>::add(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator+= (
    const ValueType & rhs)
{
    _data+=rhs;
    //_f=__fun_traits<FunctionType>::add(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data*=rhs._data;
    //_f=__fun_traits<FunctionType>::multiply(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator*= (
    const ValueType & rhs)
{
    _data*=rhs;
    //_f=__fun_traits<FunctionType>::multiply(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data/=rhs._data;
    //_f=__fun_traits<FunctionType>::divide(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator/= (
    const ValueType & rhs)
{
    _data/=rhs;
    //_f=__fun_traits<FunctionType>::divide(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const GridObject<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs._grids == _grids, "Grid mismatch");
    _data-=rhs._data;
    //_f=__fun_traits<FunctionType>::subtract(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...>& GridObject<ValueType,GridTypes...>::operator-= (
    const ValueType & rhs)
{
    _data-=rhs;
    //_f=__fun_traits<FunctionType>::subtract(_f, __fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator+ (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator+ (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out+=rhs;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator* (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator* (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator/ (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out/=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator/ (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out/=rhs;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator- (
    const GridObject<ValueType,GridTypes...>& rhs) const
{
    GridObject out(*this);
    out-=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
GridObject<ValueType,GridTypes...> GridObject<ValueType,GridTypes...>::operator- (
    const ValueType & rhs) const
{
    GridObject out(*this);
    out-=rhs;
    return out;
}


} // end of namespace GFTools
#endif
