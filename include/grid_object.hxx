#ifndef ___GFTOOLS_GRIDOBJECT_HXX___
#define ___GFTOOLS_GRIDOBJECT_HXX___

#include "grid_object.hpp"
#include <iomanip>
//#include <boost/foreach.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>

namespace gftools {

//
// grid_object::containerExtractor
//

// containerExtactor::get_ref

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes>
inline ValueType& grid_object<ValueType,GridTypes...>::containerExtractor<Nc,CT,ArgType1,ArgTypes...>::get_ref(
    CT &data, const std::tuple<GridTypes...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    static_assert(std::is_same<ArgType1,typename std::tuple_element<N-Nc,std::tuple<GridTypes...>>::type::point>::value,"Argument to a reference operator is not a point");
    if (!grid.checkPoint(arg1)) throw (exPointMismatch());
    auto &tmp = data[arg1.index_];
    return containerExtractor<Nc-1, decltype(tmp), ArgTypes...>::get_ref(tmp,grids,args...);
};

template< typename ValueType, typename ...GridTypes>
template <typename CT, typename ArgType1> 
inline ValueType& grid_object<ValueType,GridTypes...>::containerExtractor<1,CT,ArgType1>::get_ref(
    CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    static_assert(std::is_same<ArgType1,typename std::tuple_element<N-1,std::tuple<GridTypes...>>::type::point>::value,"Argument to a reference operator is not a point");
    if (!grid.checkPoint(arg1)) throw (exPointMismatch());
    auto &tmp = data[arg1.index_];
    return tmp;
}

// containerExtactor::get

template< typename ValueType, typename ...GridTypes> 
template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes>
inline ValueType grid_object<ValueType,GridTypes...>::containerExtractor<Nc,CT,ArgType1,ArgTypes...>::get(
    CT &data, const std::tuple<GridTypes...> &grids, 
    const ArgType1& arg1, const ArgTypes&... args) 
{
    const auto & grid=std::get<N-Nc>(grids);
    auto tmp = grid.evaluate(data, arg1);
    return containerExtractor<Nc-1, decltype(tmp), ArgTypes...>::get(tmp,grids,args...);
};

template< typename ValueType, typename ...GridTypes>
template <typename CT, typename ArgType1> 
inline ValueType grid_object<ValueType,GridTypes...>::containerExtractor<1,CT,ArgType1>::get(
    CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1)
{
    const auto & grid=std::get<N-1>(grids);
    auto tmp = grid.evaluate(data, arg1);
    return tmp;
}

//
// grid_object
//
    
template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>::grid_object( const std::tuple<GridTypes...> &grids, const container<ValueType, sizeof...(GridTypes)>& data):
    grids_(grids),
    dims_(GetGridSizes<N>::TupleSizeToArray(grids_)),
    data_(data),
    _f(tools::fun_traits<FunctionType>::constant(0.0)) 
{
};

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>::grid_object( const std::tuple<GridTypes...> &in):
    grids_(in),
    dims_(GetGridSizes<N>::TupleSizeToArray(grids_)),
    data_(dims_),
    _f(tools::fun_traits<FunctionType>::constant(0.0))
{
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>::grid_object( const grid_object<ValueType, GridTypes...>& rhs):
    grids_(rhs.grids_), 
    dims_(rhs.dims_),
    data_(rhs.data_),
    _f(rhs._f)
{
}; 

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>::grid_object( grid_object<ValueType,GridTypes...> && rhs):
    grids_(rhs.grids_),
    data_(rhs.data_)
{
    dims_.swap(rhs.dims_);
    _f.swap(rhs._f);
}

template <typename ValueType, typename ...GridTypes> 
const std::tuple<GridTypes...> grid_object<ValueType,GridTypes...>::getGrids() const 
{ 
    return grids_; 
};

template <typename ValueType, typename ...GridTypes> 
inline auto grid_object<ValueType,GridTypes...>::operator[](size_t i)->decltype(data_[i])
{
    return data_[i];
}

/*
template <typename ValueType, typename ...GridTypes> 
template <int M> 
ValueType& grid_object<ValueType,GridTypes...>::operator[](const std::array<size_t,M>& in)
{
}*/

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& grid_object<ValueType,GridTypes...>::get(const ArgTypes&... in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "grid_object call, number of input parameters mismatch."); 
    return containerExtractor<sizeof...(GridTypes), container<ValueType, N>, ArgTypes...>::get_ref(data_,grids_,in...);
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType& grid_object<ValueType,GridTypes...>::get(const std::tuple<ArgTypes...> &in)
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "grid_object call, number of input parameters mismatch."); 
    std::function<ValueType&(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType&{return this->get<ArgTypes...>(in1...); };// containerExtractor<sizeof...(GridTypes), ArgTypes...>::get(data_,grids_,in...);};
    tuple_tools::extra::tuple_caller<ValueType&,ArgTypes...> t = {in,f1};
    return t.call();
//containerExtractor<sizeof...(GridTypes), container<ValueType, N>, ArgTypes...>::get_ref(data_,grids_,in);
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType& grid_object<ValueType,GridTypes...>::get(const PointTupleType& in)
{
        auto indices = getIndicesFromPoints(in);
        return data_(indices);
}



template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType grid_object<ValueType,GridTypes...>::operator()(const ArgTypes&... in) const
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "grid_object call, number of input parameters mismatch."); 
    try { return containerExtractor<sizeof...(GridTypes), container<ValueType, N>, ArgTypes...>::get(data_,grids_,in...); }
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
inline ValueType grid_object<ValueType,GridTypes...>::operator()(const std::tuple<ArgTypes...>& in) const
{
    static_assert(sizeof...(ArgTypes) == sizeof...(GridTypes), "grid_object call, number of input parameters mismatch."); 
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return this->template operator()<ArgTypes...>(in1...); };
    return tuple_tools::unfold_tuple(f1,in);
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType grid_object<ValueType,GridTypes...>::operator()(const PointTupleType& in) const
{
    try {
        auto indices = getIndicesFromPoints(in);
        return data_(indices);
        }
    catch (ex_wrong_index) {
        return this->operator()(ArgTupleType(in));
        };
}



template <typename ValueType, typename ...GridTypes> 
std::ostream& operator<<(std::ostream& lhs, const grid_object<ValueType,GridTypes...> &in)
{
    lhs << (in.data_);
    return lhs;
}

template <typename ValueType, typename ...GridTypes> 
template <size_t M>
auto grid_object<ValueType,GridTypes...>::getGrid() const -> const typename std::tuple_element<M, std::tuple<GridTypes...>>::type& 
{
    return std::get<M>(grids_);
}



template <typename ValueType, typename ...GridTypes> 
auto grid_object<ValueType,GridTypes...>::getGrid() const -> const typename std::tuple_element<0, std::tuple<GridTypes...>>::type&
{
    return std::get<0>(grids_);
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
void grid_object<ValueType,GridTypes...>::fill(const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(FunctionType(in));
}

/*
template <typename ValueType, typename ...GridTypes> 
inline void grid_object<ValueType,GridTypes...>::fill(const typename grid_object<ValueType,GridTypes...>::PointFunctionType& in)
{

    __gencontainerExtractor<sizeof...(GridTypes), containerExtractor, std::tuple<GridTypes...>>::type::set(data_,grids_,in);

    //containerExtractor<sizeof...(GridTypes), ArgTypes...>::set(data_,grids_,in);
    //_f = in;
}
*/

template <typename ValueType, typename ...GridTypes> 
inline size_t grid_object<ValueType,GridTypes...>::getTotalcontainerSize() const
{
    size_t out = 1;
    for (auto i : dims_) out*=i;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename grid_object<ValueType,GridTypes...>::PointIndices grid_object<ValueType,GridTypes...>::_getPointsIndices(const size_t index) const
{
    PointIndices indices;
    size_t t = index;
    for (int i=N-1; i>=0; i--) { 
        indices[i]=t%dims_[i];
        t-=indices[i];
        t/=dims_[i];
        }
    return indices;
}


template <typename ValueType, typename ...GridTypes>
inline typename grid_object<ValueType,GridTypes...>::ArgTupleType grid_object<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    return this->getArgsFromIndices<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_object<ValueType,GridTypes...>::ArgTupleType grid_object<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    ArgTupleType out;
    auto t1 = std::get<N-1>(grids_)[in[N-1]].val_;
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_object<ValueType,GridTypes...>::ArgTupleType grid_object<ValueType,GridTypes...>::getArgsFromIndices(PointIndices in) const
{
    auto out = getArgsFromIndices<M-1>(in);
    auto t1 = std::get<N-1-M>(grids_)[in[N-1-M]];
    std::get<N-1-M>(out) = t1.val_;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename grid_object<ValueType,GridTypes...>::PointTupleType grid_object<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    return this->getPointsFromIndices<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_object<ValueType,GridTypes...>::PointTupleType grid_object<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    PointTupleType out;
    auto t1 = std::get<N-1>(grids_)[in[N-1]];
    std::get<N-1>(out)=t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_object<ValueType,GridTypes...>::PointTupleType grid_object<ValueType,GridTypes...>::getPointsFromIndices(PointIndices in) const
{
    auto out = getPointsFromIndices<M-1>(in);
    auto t1 = std::get<N-1-M>(grids_)[in[N-1-M]];
    std::get<N-1-M>(out) = t1;
    return out;
}

template <typename ValueType, typename ...GridTypes>
inline typename grid_object<ValueType,GridTypes...>::PointIndices grid_object<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    return this->getIndicesFromPoints<N-1>(in);
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M ==0, bool>::type>
inline typename grid_object<ValueType,GridTypes...>::PointIndices grid_object<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    PointIndices out;
    auto t1 = std::get<N-1>(in);
    if (std::get<N-1>(grids_).size()<=t1.index_) throw ex_wrong_index();
    if (std::get<N-1>(grids_)[t1.index_].index_ != t1.index_) { throw ex_wrong_index(); };
    out[N-1]=t1.index_;
    return out;
}

template <typename ValueType, typename ...GridTypes>
template <int M, typename std::enable_if<M >= 1, bool>::type >
inline typename grid_object<ValueType,GridTypes...>::PointIndices grid_object<ValueType,GridTypes...>::getIndicesFromPoints(PointTupleType in) const
{
    auto out = getIndicesFromPoints<M-1>(in);
    auto t1 = std::get<N-1-M>(in);
    if (std::get<N-1-M>(grids_).size()<=t1.index_) throw ex_wrong_index();
    if (std::get<N-1-M>(grids_)[t1.index_].index_ != t1.index_) throw ex_wrong_index();
    out[N-1-M]=t1.index_;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::fill_tuple(const std::function<ValueType(ArgTupleType)>& in)
{
    size_t total_size = this->getTotalcontainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        ArgTupleType args = this->getArgsFromIndices(pts_index);
        auto val = in(args);
        data_(pts_index) = val;
        };
    _f = tools::fun_traits<FunctionType>::getFromTupleF(in); 

}

template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::fill(const typename grid_object<ValueType,GridTypes...>::FunctionType& in)
{
    size_t total_size = this->getTotalcontainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        ArgTupleType args = this->getArgsFromIndices(pts_index);
        auto val = tuple_tools::unfold_tuple(in, args);
        data_(pts_index) = val;
        };
    _f = in;
}

template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::fill(const typename grid_object<ValueType,GridTypes...>::PointFunctionType& in)
{
    size_t total_size = this->getTotalcontainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        PointTupleType args = this->getPointsFromIndices(pts_index);
        auto val = tuple_tools::unfold_tuple(in, args);
        data_(pts_index) = val;
        };
}

template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::fill_tuple(const std::function<ValueType(PointTupleType)>& in)
{
    size_t total_size = this->getTotalcontainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        PointTupleType args = this->getPointsFromIndices(pts_index);
        auto val = in(args);
        data_(pts_index) = val;
        };
    //_f = tools::fun_traits<FunctionType>::getFromTupleF(std::function<ValueType(ArgTupleType)>(in)); 
}

template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, complex_type>::value, int>::type>
real_type grid_object<ValueType,GridTypes...>::diff(const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object outObj(grids_);
    auto f1 = [&](PointTupleType in){return std::abs((*this)(in) - rhs(in));};
    PointFunctionType f = tools::fun_traits<PointFunctionType>::getFromTupleF(f1); 
    outObj.fill(f);
    real_type norm = 1.0;
    for (auto v : dims_) { norm*=v; };
    return std::real(outObj.sum())/norm;
}

template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_same<U, real_type>::value, int>::type>
real_type grid_object<ValueType,GridTypes...>::diff(const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object outObj(grids_);
    auto f1 = [&](PointTupleType in){return std::abs((*this)(in) - rhs(in));};
    PointFunctionType f = tools::fun_traits<PointFunctionType>::getFromTupleF(f1); 
    outObj.fill(f);
    real_type norm = 1.0;
    for (auto v : dims_) { norm*=v; };
    return std::real(outObj.sum())/norm;
}


template <typename ValueType, typename ...GridTypes> 
template <typename U, typename std::enable_if<std::is_convertible<U, complex_type>::value, int>::type>
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::conj()
{
    grid_object<ValueType,GridTypes...> out(*this);
    (out.data_) = out.data_.conj();
    return out;
}

template <typename ValueType, typename ...GridTypes> 
inline ValueType grid_object<ValueType,GridTypes...>::sum()
{
    return data_.sum();
}

template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline ValueType grid_object<ValueType,GridTypes...>::__get_f(const std::tuple<ArgTypes...>& in) const
{
    std::function<ValueType(ArgTypes...)> f1 = [&](ArgTypes... in1)->ValueType{return this->_f(in1...); };
    tuple_tools::extra::tuple_caller<ValueType,ArgTypes...> t = {in,f1};
    return t.call();
}



template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::shift(ArgTypes... args) const
{
    return this->shift(std::forward_as_tuple(args...));
}


template <typename ValueType, typename ...GridTypes> 
template <typename ...ArgTypes> 
inline grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::shift(const std::tuple<ArgTypes...>& shift_args) const
{
    grid_object<ValueType,GridTypes...> out(grids_);
    std::function<ValueType(PointTupleType)> ShiftFunction = [&](PointTupleType args1)->ValueType { 
        PointTupleType out_args = this->_shiftArgs(args1, shift_args);
    //    __tuple_print<PointTupleType>::print(args1); 
    //    INFO_NONEWLINE("+");  __tuple_print<std::tuple<ArgTypes...>>::print(shift_args); 
    //    INFO_NONEWLINE("-->");__tuple_print<PointTupleType>::print(out_args);
        return (*this)(out_args);
        };
    PointFunctionType fillF = tools::fun_traits<PointFunctionType>::getFromTupleF(ShiftFunction);
    //out.fill(fillF);
    out.fill_tuple(ShiftFunction);
    
    static std::function<ValueType(ArgTupleType)> ShiftAnalyticF;
    ShiftAnalyticF = [this, shift_args](const ArgTupleType& in)->ValueType {
        ArgTupleType out_args = _shiftArgs(in,shift_args); 
        return __get_f(out_args);
    };
    
    FunctionType tailF = tools::fun_traits<FunctionType>::getFromTupleF(ShiftAnalyticF);
    out._f = tailF;
    
    return out;
}

/*
template <typename ValueType, typename ...GridTypes> 
    template <typename OrigArg1, typename ...OrigArgs, typename ArgType1, typename ...ArgTypes, 
        typename std::enable_if<sizeof...(OrigArgs)==sizeof...(ArgTypes)>,  
        typename std::enable_if<sizeof...(OrigArgs)!=0> > 
std::tuple<OrigArg1, OrigArgs...> grid_object<ValueType,GridTypes...>::_shiftArgs(const std::tuple<OrigArg1, OrigArgs...>&in, const std::tuple<ArgType1, ArgTypes...>& shift_args) const
{
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-sizeof...(ArgTypes)-1>(grids_).shift(arg1,shift_arg1); 


    auto f_o = [&in](OrigArg1 arg1, OrigArgs... others)->std::tuple<OrigArgs...>{ return std::forward_as_tuple(others...);};
    auto f_s = [&shift_args](ArgType1 arg1, ArgTypes... others)->std::tuple<ArgTypes...>{ return std::forward_as_tuple(others...);};

    tuple_tools::extra::tuple_caller<std::tuple<OrigArgs...>,OrigArg1,OrigArgs...> t_o = {in,f_o};
    std::tuple<OrigArgs...> other_orig_args(t_o.call());

    tuple_tools::extra::tuple_caller<std::tuple<ArgTypes...>,ArgType1,ArgTypes...> t_s = {shift_args,f_s};
    std::tuple<ArgTypes...> other_shift_args(t_s.call());

    return std::tuple_cat(std::forward_as_tuple(out1),this->_shiftArgs(other_orig_args,other_shift_args));
}
*/
template <typename ValueType, typename ...GridTypes> 
template <typename OrigArg1, typename ArgType1> 
inline std::tuple<OrigArg1> grid_object<ValueType,GridTypes...>::_shiftArgs(const std::tuple<OrigArg1>&in, const std::tuple<ArgType1>& shift_args) const
{
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-1>(grids_).shift(arg1,shift_arg1); 
    return std::forward_as_tuple(out1);
}

template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::savetxt(const std::string& fname) const
{
    INFO("Saving " << typeid(*this).name() << " to " << fname);
    std::ofstream out;
    out.open(fname.c_str());
    size_t total_size = this->getTotalcontainerSize();
    size_t last_grid_size = std::get<N-1>(grids_).size();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);
        //ArgTupleType args = this->getArgsFromIndices(pts_index);
        PointTupleType pts = this->getPointsFromIndices(pts_index);
        auto val = (*this)(pts);
        out << std::scientific << tuple_tools::serialize_tuple(pts) << "    " << num_io<ValueType>(val) << std::endl;
        if (N > 1 && i && (i+1)%last_grid_size==0) out << std::endl;
        };
    out.close();
}

template <typename ValueType, typename ...GridTypes> 
void grid_object<ValueType,GridTypes...>::loadtxt(const std::string& fname, real_type tol)
{
    INFO("Loading " << typeid(*this).name() << " from " << fname);
    std::ifstream in;
    in.open(fname.c_str());
    if (in.fail()) { ERROR("Couldn't open file " << fname); throw exIOProblem(); };
    size_t total_size = this->getTotalcontainerSize();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = _getPointsIndices(i);

        PointTupleType pts = this->getPointsFromIndices(pts_index);
        //ArgTupleType args = this->getArgsFromIndices(pts_index);
        ArgTupleType pts2 = tuple_tools::read_tuple<ArgTupleType>(in);
        if (!tools::is_float_equal<ArgTupleType>(pts,pts2,tol)) throw (exIOProblem());

        num_io<ValueType> tmp2(this->get(pts));
        in >> tmp2;
        this->get(pts) = tmp2.value_;
        };


    in.close();
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::copyInterpolate (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //data_=rhs.data_;
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
inline grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator= (
    const std::function<ValueType(ArgTypes...)> & in)
{
    this->fill(in);
    return *this;
}
*/

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator= (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_=rhs.data_;
    _f = rhs._f;
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator= (
    const ValueType& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_=rhs;
    _f = tools::fun_traits<FunctionType>::constant(rhs);
    return *this;
}



template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator+= (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_+=rhs.data_;
    //_f=tools::fun_traits<FunctionType>::add(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator+= (
    const ValueType & rhs)
{
    data_+=rhs;
    //_f=tools::fun_traits<FunctionType>::add(_f, tools::fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator*= (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_*=rhs.data_;
    //_f=tools::fun_traits<FunctionType>::multiply(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator*= (
    const ValueType & rhs)
{
    data_*=rhs;
    //_f=tools::fun_traits<FunctionType>::multiply(_f, tools::fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator/= (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_/=rhs.data_;
    //_f=tools::fun_traits<FunctionType>::divide(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator/= (
    const ValueType & rhs)
{
    data_/=rhs;
    //_f=tools::fun_traits<FunctionType>::divide(_f, tools::fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator-= (
    const grid_object<ValueType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_-=rhs.data_;
    //_f=tools::fun_traits<FunctionType>::subtract(_f, rhs._f);
    return *this;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...>& grid_object<ValueType,GridTypes...>::operator-= (
    const ValueType & rhs)
{
    data_-=rhs;
    //_f=tools::fun_traits<FunctionType>::subtract(_f, tools::fun_traits<FunctionType>::constant(rhs));
    return *this;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator+ (
    const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object out(*this);
    out+=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator+ (
    const ValueType & rhs) const
{
    grid_object out(*this);
    out+=rhs;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator* (
    const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator* (
    const ValueType & rhs) const
{
    grid_object out(*this);
    out*=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator/ (
    const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object out(*this);
    out/=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator/ (
    const ValueType & rhs) const
{
    grid_object out(*this);
    out/=rhs;
    return out;
}


template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator- (
    const grid_object<ValueType,GridTypes...>& rhs) const
{
    grid_object out(*this);
    out-=rhs;
    return out;
}

template <typename ValueType, typename ...GridTypes> 
grid_object<ValueType,GridTypes...> grid_object<ValueType,GridTypes...>::operator- (
    const ValueType & rhs) const
{
    grid_object out(*this);
    out-=rhs;
    return out;
}


} // end of namespace GFTools
#endif
