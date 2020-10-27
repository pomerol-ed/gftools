#pragma once

#include <boost/functional/hash.hpp>
#include <unordered_set>
#include "grid_object.hpp"

namespace gftools {

//
// grid_object_base
//

//
// constructors
//

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( const std::tuple<GridTypes...> &grids):
    grids_(grids),
    dims_(trs::get_dimensions(grids)),
    data_(dims_),
    tail_(tools::fun_traits<function_type>::constant(0.0))
{
}

template <typename ContainerType, typename ...GridTypes>
template <typename CType>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( const std::tuple<GridTypes...> &grids, CType& data):
    grids_(grids),
    dims_(trs::get_dimensions(grids)),
    data_(data),
    tail_(tools::fun_traits<function_type>::constant(0.0))
{
    if (dims_ != data.shape()) throw gftools::ex_generic("Dimensions mismatch when creating grid_object from existing data");
};

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( const std::tuple<GridTypes...> &grids, ContainerType&& data):
    grids_(grids),
    dims_(trs::get_dimensions(grids)),
    data_(std::forward<ContainerType>(data)),
    tail_(tools::fun_traits<function_type>::constant(0.0))
{
    if (dims_ != data.shape()) throw gftools::ex_generic("Dimensions mismatch when creating grid_object from existing data");
};

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( grid_object_base<ContainerType,GridTypes...> && rhs):
    grids_(rhs.grids_),
    dims_(rhs.dims_),
    data_(std::forward<ContainerType>(rhs.data_))
{
    tail_.swap(rhs.tail_);
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( const grid_object_base<ContainerType, GridTypes...>& rhs):
    grids_(rhs.grids_),
    dims_(rhs.dims_),
    data_(rhs.data_),
    tail_(rhs.tail_)
{
};


template <typename ContainerType, typename ...GridTypes>
template <typename CType>
grid_object_base<ContainerType,GridTypes...>::grid_object_base( const grid_object_base<CType, GridTypes...>& rhs):
    grids_(rhs.grids_),
    dims_(rhs.dims_),
    data_(rhs.data_),
    tail_(rhs.tail_)
{
}
//
// Assignment
//

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator= (
    const grid_object_base<ContainerType,GridTypes...>& rhs)
{
    if (rhs.size() != this->size()) throw gftools::ex_generic("Assigning objects of different sizes");
    assert(dims_ == rhs.dims_);
    data_=rhs.data_;
    tail_ = rhs.tail_;
    return *this;
}


template <typename ContainerType, typename ...GridTypes>
    template <typename CType>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator= (
    const grid_object_base<CType,GridTypes...>& rhs)
{
    if (rhs.size() != this->size()) throw gftools::ex_generic("Assigning objects of different sizes");
    data_=rhs.data_;
    tail_ = rhs.tail_;
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator= (
    const value_type& rhs)
{
    data_=rhs;
    tail_ = tools::fun_traits<function_type>::constant(rhs);
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator= (
    grid_object_base<ContainerType,GridTypes...>&& rhs)
{
    if (rhs.size() != this->size()) throw gftools::ex_generic("Assigning objects of different sizes");
    data_.swap(rhs.data_);
    tail_.swap(rhs.tail_);
    return *this;
}


//
// shift
//

template <typename ContainerType, typename ...GridTypes>
template <typename ...ArgTypes>
typename std::enable_if<
    (std::is_convertible<std::tuple<ArgTypes...>, typename grid_object_base<ContainerType,GridTypes...>::arg_tuple>::value ||
     std::is_same<std::tuple<ArgTypes...>, typename grid_object_base<ContainerType,GridTypes...>::point_tuple>::value),
    grid_object<typename grid_object_base<ContainerType,GridTypes...>::value_type, GridTypes...>>::type
        grid_object_base<ContainerType,GridTypes...>::shift (const std::tuple<ArgTypes...>& shift_args) const
{
    grid_object_base<ContainerType,GridTypes...> out(grids_);
    std::function<value_type(point_tuple)> ShiftFunction = [&](point_tuple args1)->value_type {
        try { point_tuple out_args = trs::shift(args1, shift_args,grids_); return (*this)(out_args); }
        catch (...)
            { arg_tuple out_args = trs::shift(arg_tuple(args1), arg_tuple(shift_args) ,grids_); return (*this)(out_args); }
    //    __tuple_print<point_tuple>::print(args1);
    //    INFO_NONEWLINE("+");  __tuple_print<std::tuple<ArgTypes...>>::print(shift_args);
    //    INFO_NONEWLINE("-->");__tuple_print<point_tuple>::print(out_args);
        };
    out.fill(ShiftFunction);

    static std::function<value_type(arg_tuple)> ShiftAnalyticF;
    ShiftAnalyticF = [this, shift_args](const arg_tuple& in)->value_type {
        arg_tuple out_args = trs::shift(in,shift_args,grids_);
        return tuple_tools::unfold_tuple(this->tail_, out_args);
    };

    function_type tailF = tools::extract_tuple_f(ShiftAnalyticF);
    out.tail_ = tailF;

    return out;
}

//
// IO
//

namespace extra {
template <size_t D>
inline std::array<size_t, D> enumerate_indices_(const size_t index, const std::array<size_t, D> dims)
{
    std::array<size_t, D> indices;
    size_t t = index;
    for (int i=D-1; i>=0; i--) {
        indices[i]=t%dims[i];
        t-=indices[i];
        t/=dims[i];
        }
    return indices;
}
} // end of namespace extra

template <typename ContainerType, typename ...GridTypes>
std::ostream& operator<<(std::ostream& lhs, const grid_object_base<ContainerType,GridTypes...> &in)
{
    lhs << (in.data_);
    return lhs;
}



template <typename ContainerType, typename ...GridTypes>
void grid_object_base<ContainerType,GridTypes...>::savetxt(const std::string& fname) const
{
    INFO("Saving " << demangled_name(typeid(*this)) << " to " << fname);
    std::ofstream out;
    out.open(fname.c_str());
    size_t total_size = this->size();
    size_t last_grid_size = std::get<N-1>(grids_).size();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = extra::enumerate_indices_(i, dims_);
        //arg_tuple args = this->getArgsFromIndices(pts_index);
        point_tuple pts = this->points(pts_index);
        auto val = data_(pts_index);
        out << std::scientific << tuple_tools::serialize_tuple<arg_tuple>(pts) << "    " << num_io<value_type>(val) << std::endl;
        if (N > 1 && i && (i+1)%last_grid_size==0) out << std::endl;
        };
    out.close();
}

template <typename ContainerType, typename ...GridTypes>
void grid_object_base<ContainerType,GridTypes...>::loadtxt(const std::string& fname, real_type tol)
{
    INFO("Loading " << demangled_name(typeid(*this)) << " from " << fname);
    std::ifstream in;
    in.open(fname.c_str());
    if (in.fail()) { ERROR("Couldn't open file " << fname); throw exIOProblem(); };
    size_t total_size = this->size();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = extra::enumerate_indices_(i, dims_);
        point_tuple pts = this->points(pts_index);
        arg_tuple pts2 = tuple_tools::read_tuple<arg_tuple>(in); // ensure serialize_tuple in savetxt has the same type. Dropping here indices - they're wrong anyway.
        if (!tools::is_float_equal<arg_tuple>(pts,pts2,tol)) throw (exIOProblem());

        num_io<value_type> tmp2(data_(pts_index));
        in >> tmp2;
        data_(pts_index) = tmp2();
        };


    in.close();
}

//
// Fill values
//

template <typename ContainerType, typename ...GridTypes>
void grid_object_base<ContainerType,GridTypes...>::fill_point_function(const typename grid_object_base<ContainerType,GridTypes...>::point_function_type& in)
{
    size_t total_size = this->size();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = extra::enumerate_indices_(i, dims_);
        point_tuple args = trs::points(pts_index,grids_);
        auto val = tuple_tools::unfold_tuple(in, args);
        data_(pts_index) = val;
        };
}

template <typename ContainerType, typename ...GridTypes>
void grid_object_base<ContainerType,GridTypes...>::fill_function(const typename grid_object_base<ContainerType,GridTypes...>::function_type& in)
{
    size_t total_size = this->size();
    for (size_t i=0; i<total_size; ++i) {
        auto pts_index = extra::enumerate_indices_(i, dims_);
        arg_tuple args = trs::get_args(pts_index,grids_);
        auto val = tuple_tools::unfold_tuple(in, args);
        data_(pts_index) = val;
        };
    tail_ = in;
}

template <typename ContainerType, typename ...GridTypes>
template <typename CType2>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::copy_interpolate (
    const grid_object_base<CType2,GridTypes...>& rhs)
{
    tail_ = rhs.tail_;
    const std::function<value_type(arg_tuple)> bindf = [&](arg_tuple in){return rhs(in);};
    this->fill(bindf);
    return (*this);
}




//
// Math (obsolete)
//
/*
template <typename ContainerType, typename ...GridTypes>
template <typename ...ArgTypes>
inline grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator= (
    const std::function<value_type(ArgTypes...)> & in)
{
    this->fill(in);
    return *this;
}
*/

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator+= (
    const grid_object_base<ContainerType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_+=rhs.data_;
    //_f=tools::fun_traits<function_type>::add(_f, rhs._f);
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator+= (
    const value_type & rhs)
{
    data_+=rhs;
    //tail_=tools::fun_traits<function_type>::add(tail_, tools::fun_traits<function_type>::constant(rhs));
    return *this;
}


template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator*= (
    const grid_object_base<ContainerType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_*=rhs.data_;
    //tail_=tools::fun_traits<function_type>::multiply(tail_, rhs.tail_);
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator*= (
    const value_type & rhs)
{
    data_*=rhs;
    //tail_=tools::fun_traits<function_type>::multiply(tail_, tools::fun_traits<function_type>::constant(rhs));
    return *this;
}


template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator/= (
    const grid_object_base<ContainerType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_/=rhs.data_;
    //tail_=tools::fun_traits<function_type>::divide(tail_, rhs.tail_);
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator/= (
    const value_type & rhs)
{
    data_/=rhs;
    //tail_=tools::fun_traits<function_type>::divide(tail_, tools::fun_traits<function_type>::constant(rhs));
    return *this;
}


template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator-= (
    const grid_object_base<ContainerType,GridTypes...>& rhs)
{
    //static_assert(rhs.grids_ == grids_, "Grid mismatch");
    data_-=rhs.data_;
    //tail_=tools::fun_traits<function_type>::subtract(tail_, rhs.tail_);
    return *this;
}

template <typename ContainerType, typename ...GridTypes>
grid_object_base<ContainerType,GridTypes...>& grid_object_base<ContainerType,GridTypes...>::operator-= (
    const value_type & rhs)
{
    data_-=rhs;
    //tail_=tools::fun_traits<function_type>::subtract(tail_, tools::fun_traits<function_type>::constant(rhs));
    return *this;
}


namespace extra {

    template <typename GridTupleType, typename OrigGridTupleType, int N>
    struct GridTupleFromVector
    {
        typedef std::vector<typename tools::grid_tuple_traits<OrigGridTupleType>::arg_tuple> vector_type;
        static constexpr int TupleSize = std::tuple_size<OrigGridTupleType>::value;
        typedef typename std::tuple_element<TupleSize - N, OrigGridTupleType>::type grid_type;
        typedef typename grid_type::value_type value_type;
        typedef typename tuple_tools::extra::tuple_io<GridTupleType>::tail_type grid_tuple_tail;
        static GridTupleType extract(vector_type const& in) {
            std::unordered_set<value_type, boost::hash<value_type>> vals_set;
            std::vector<value_type> vals_vec;
            for (size_t i=0; i<in.size(); i++) {
                value_type v = std::get<TupleSize - N>(in[i]);
                auto x = vals_set.insert(v);
                if (x.second) vals_vec.push_back(v);
                }
            return std::tuple_cat(
                std::forward_as_tuple(grid_type(vals_vec)),
                GridTupleFromVector<grid_tuple_tail, OrigGridTupleType, N-1>::extract(in)
                );
            }
    };

    template <typename GridTupleType, typename OrigGridTupleType>
    struct GridTupleFromVector<GridTupleType, OrigGridTupleType, 1>
    {
        typedef std::vector<typename tools::grid_tuple_traits<OrigGridTupleType>::arg_tuple> vector_type;
        static constexpr int TupleSize = std::tuple_size<OrigGridTupleType>::value;
        typedef typename std::tuple_element<TupleSize - 1, OrigGridTupleType>::type grid_type;
        typedef typename grid_type::value_type value_type;
        static constexpr int N = 1;
        static GridTupleType extract(vector_type const& in) {
            static_assert(std::tuple_size<GridTupleType>::value == 1, "wrong tuple size");
            std::unordered_set<value_type, boost::hash<value_type>> vals_set;
            std::vector<value_type> vals_vec;
            for (size_t i=0; i<in.size(); i++) {
                value_type v = std::get<TupleSize - N>(in[i]);
                auto x = vals_set.insert(v);
                if (x.second) vals_vec.push_back(v);
                }
            return std::forward_as_tuple(grid_type(vals_vec));
            }
    };



    template <typename GridTuple>
    GridTuple arg_vector_to_grid_tuple(std::vector<typename tools::grid_tuple_traits<GridTuple>::arg_tuple> const& in)
    {
        constexpr int ngrids =  std::tuple_size<GridTuple>::value;
        return GridTupleFromVector<GridTuple,GridTuple,ngrids>::extract(in);
    }
}

template <typename GridObjectType>
GridObjectType loadtxt(std::string const& fname, double tol)
{
// load some grid_object from file
    typedef typename GridObjectType::arg_tuple arg_tuple;
    typedef typename GridObjectType::value_type value_type;
    typedef typename GridObjectType::grid_tuple grid_tuple;

    INFO("Loading " << demangled_name(typeid(GridObjectType)) << " from " << fname);
    std::ifstream in;
    in.open(fname.c_str());
    if (in.fail()) { ERROR("Couldn't open file " << fname); throw typename GridObjectType::exIOProblem(); };

    std::vector<arg_tuple> grid_vals;
    std::vector<value_type> vals;

    while (!in.eof()) {
        arg_tuple pts = tuple_tools::read_tuple<arg_tuple>(in); // ensure serialize_tuple in savetxt has the same type. Dropping here indices - they're wrong anyway.
        if (in.eof()) break;
        value_type v;
        in >> num_io<value_type>(v);
        grid_vals.push_back(pts);
        vals.push_back(v);
        }

//    for (int i=0; i<grid_vals.size(); i++)
//        DEBUG(tuple_tools::print_tuple(grid_vals[i]));

    grid_tuple grids = extra::arg_vector_to_grid_tuple<grid_tuple>(grid_vals);
    GridObjectType out(grids);
    std::copy(vals.begin(), vals.end(), out.data().data());
    return out;
}


} // end of namespace GFTools
