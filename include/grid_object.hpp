#pragma once
#ifndef ___GFTOOLS_GRID_OBJECT_HPP___
#define ___GFTOOLS_GRID_OBJECT_HPP___

#include "defaults.hpp"
#include "tools.hpp"
#include "grid_base.hpp"
#include "grid_tools.hpp"
#include "container.hpp"

namespace gftools {

/** A grid_object is a wrapper over container class, that stores data,
 * defined on multiple different grids.
 */
template< typename ValueType, typename ...GridTypes> 
class grid_object 
{
public:
    /** A typedef for a function that gives the analytical value of the object, when it's not stored. */
    typedef typename GridArgTypeExtractor<ValueType, std::tuple<GridTypes...> >::arg_type FunctionType;
    /** A typedef for a function that gives the analytical value of the object, when it's not stored. */
    typedef typename GridPointExtractor<ValueType, std::tuple<GridTypes...> >::point_type PointFunctionType;
    /** A typedef for a tuple of grids. */
    typedef std::tuple<GridTypes...> GridTupleType;
    /** A typedef for a tuple of grid points. */
    typedef typename GridPointExtractor<ValueType, std::tuple<GridTypes...> >::arg_tuple_type PointTupleType;
    /** A typedef for a tuple of grid point values. */
    typedef typename GridArgTypeExtractor<ValueType, std::tuple<GridTypes...> >::arg_tuple_type ArgTupleType;
    /** A typedef for the values stored in the container. */
    typedef ValueType value_type;
    static constexpr size_t N = sizeof...(GridTypes);
    typedef std::array<size_t, N> PointIndices;
    class exPointMismatch : public std::exception { virtual const char* what() const throw() { return "Index mismatch."; }; };
protected:
    /** Grids on which the data is defined. */
    mutable std::tuple<GridTypes...> grids_;
public:
    /** The dimensions of the container - deduced from grids. */
    PointIndices dims_;
    /** A pointer to the container. A pointer is used as there exist no default 
     * constructor for the container.
     */
    //std::unique_ptr<container<ValueType, N>> data_;
    mutable container<ValueType, N> data_;

protected:
    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> PointTupleType getPointsFromIndices(PointIndices in) const;
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> PointTupleType getPointsFromIndices(PointIndices in) const;
    PointTupleType getPointsFromIndices(PointIndices in) const;

    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> ArgTupleType getArgsFromIndices(PointIndices in) const;
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> ArgTupleType getArgsFromIndices(PointIndices in) const;
    ArgTupleType getArgsFromIndices(PointIndices in) const;

    template <int M = N-1, typename std::enable_if<M >= 1, bool>::type = 0> PointIndices getIndicesFromPoints(PointTupleType in) const;
    template <int M = 0,   typename std::enable_if<M == 0, bool>::type = 0> PointIndices getIndicesFromPoints(PointTupleType in) const;
    PointIndices getIndicesFromPoints(PointTupleType in) const;


    /** A helper recursive template utility to extract and set data from the container. */
    template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes> struct containerExtractor { 
        /** Gets the data by values. */
        static ValueType get(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
        static ValueType& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
             };
    /** Specialization of containerExtractor for 1-dim container. */
    template <typename CT, typename ArgType1> struct containerExtractor<1,CT,ArgType1> {
        static ValueType get(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static ValueType get(CT &data, const std::tuple<GridTypes...> &grids, const std::tuple<ArgType1>& arg1); 

        static ValueType& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1); 
        static ValueType& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const std::tuple<ArgType1>& arg1); 
        };

    /** Returns _f(in). */
    template <typename ...ArgTypes> ValueType __get_f(const std::tuple<ArgTypes...>& in) const;
    PointIndices _getPointsIndices(const size_t index) const;

public:
    /** This function returns the value of the object when the point is not in container. */
    FunctionType _f;
    /** Constructs a grid object out of a tuple containing various grids. */
    grid_object( const std::tuple<GridTypes...> &grids);
    template <int M = N, typename std::enable_if<M == 1, bool>::type = 0>  grid_object( GridTypes... grids):grid_object(std::forward_as_tuple(grids...)){};
    /** Constructor of grids and data. */
    grid_object( const std::tuple<GridTypes...> &grids, const container<ValueType, sizeof...(GridTypes)>& data);
    /** Copy constructor. */
    grid_object( const grid_object<ValueType, GridTypes...>& rhs);
    /** Move constructor. */
    grid_object( grid_object<ValueType, GridTypes...>&& rhs);

    const std::tuple<GridTypes...> grids() const;
    size_t size() const;
    /** Returns an Mth grid in grids_. */
    template<size_t M = 0> 
        auto grid() const -> const typename std::add_lvalue_restd::tuple_element<M, std::tuple<GridTypes...>>::type& ;
    /** Returns element number i, which corresponds to (*_grid)[i]. */
    auto operator[](size_t i)->decltype(data_[i]);
    /** Const operator[]. */
    //auto operator[](size_t i) const -> const decltype(data_[i]);
    //template <size_t M> ValueType& operator[](const std::array<size_t,M>& in);
    /** Returns the data_ container. */
    container<ValueType, sizeof...(GridTypes)>& getData(){return data_;};
    /** Fills the container with a provided function. */
    template <typename ...ArgTypes> void fill(const std::function<ValueType(ArgTypes...)> &);
    void fill(const FunctionType &in);
    void fill(const PointFunctionType &in);
    //template <typename ...ArgTypes> void fill_tuple(const std::function<ValueType(const std::tuple<ArgTypes...>)> &);
    void fill_tuple(const std::function<ValueType(ArgTupleType)>& in);
    void fill_tuple(const std::function<ValueType(PointTupleType)>& in);
    /** Fills the container with any proper class with call operator. Untested */
    //template <template <typename, class> class Filler, typename ...ArgTypes> void fill(const Filler<ValueType,ArgTypes...> &);

    /** Return the value by grid values. */
    template <typename ...ArgTypes> ValueType& get(const ArgTypes&... in);
    template <typename ...ArgTypes> ValueType& get(const std::tuple<ArgTypes...>& in);
    ValueType& get(const PointTupleType& in);
    template <typename ...ArgTypes> ValueType operator()(const ArgTypes&... in) const;
    template <typename ...ArgTypes> ValueType operator()(const std::tuple<ArgTypes...>& in) const;
    ValueType operator()(const PointTupleType& in) const;
    //ValueType operator()(const ArgTupleType& in) const;
    //template <typename ...ArgTypes> auto operator()(const ArgType1& in)->decltype() const;

    /** A shortcut for fill method. */
    //template <typename ...ArgTypes> grid_object& operator= (const std::function<ValueType(ArgTypes...)> &);
    /** Same as operator=, but allows for non-equal grids. Slow. Uses analytic function to provide missing values. */
    grid_object& copyInterpolate(const grid_object &rhs);
    /** Algebraic operators. */
    grid_object& operator= (const grid_object & rhs);
    grid_object& operator= (const ValueType & rhs);
    grid_object& operator*= (const grid_object & rhs);
    grid_object& operator*= (const ValueType& rhs);
    grid_object operator* (const grid_object & rhs) const;
    grid_object operator* (const ValueType & rhs) const;
    grid_object& operator+= (const grid_object & rhs);
    grid_object& operator+= (const ValueType& rhs);
    grid_object operator+ (const grid_object & rhs) const;
    grid_object operator+ (const ValueType & rhs) const;
    grid_object& operator-= (const grid_object & rhs);
    grid_object& operator-= (const ValueType& rhs);
    grid_object operator- (const grid_object & rhs) const;
    grid_object operator- (const ValueType & rhs) const;
    grid_object& operator/= (const grid_object & rhs);
    grid_object& operator/= (const ValueType& rhs);
    grid_object operator/ (const grid_object & rhs) const;
    grid_object operator/ (const ValueType & rhs) const;
    friend grid_object operator* (const ValueType & lhs, const grid_object & rhs) {return rhs*lhs;};
    friend grid_object operator+ (const ValueType & lhs, const grid_object & rhs) {return rhs+lhs;};
    friend grid_object operator- (const ValueType & lhs, const grid_object & rhs) {return rhs*(-1.0)+lhs;};
    friend grid_object operator/ (const ValueType & lhs, const grid_object & rhs) {grid_object out(rhs); out=lhs; return out/rhs;};

    /** Returns the complex conjugate of this object, if it's complex valued. */
    template <typename U = ValueType, typename std::enable_if<std::is_convertible<U, complex_type>::value, int>::type=0>
        grid_object conj();
    /** Returns a norm of difference between two objects. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, complex_type>::value, int>::type = 0>
    real_type diff(const grid_object &rhs) const; 
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, real_type>::value, int>::type = 0>
    real_type diff(const grid_object &rhs) const;
    /** Returns the sum of all elements in the container. */
    ValueType sum();

    /** Returns an object with arguments, shifted by the given values.
     * \param[in] args A pack of arguments to shift the object
     */
    template <typename ...ArgTypes> grid_object shift(ArgTypes... args) const;
    template <typename ...ArgTypes> grid_object shift(const std::tuple<ArgTypes...>& arg_tuple) const;
    /** Save the data to the txt file. */
    void savetxt(const std::string& fname) const;
    /** Loads the data to the txt file. */
    void loadtxt(const std::string& fname, real_type tol = 1e-8);
    /** Dumps the object to the stream. */
    template <typename ValType, class ...GridTypes2> friend std::ostream& operator<<(std::ostream& lhs, const grid_object<ValType,GridTypes2...> &in);
    
    class exIOProblem : public std::exception { virtual const char* what() const throw(){return "IO problem.";} }; 
    class ex_wrong_index : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
    
/** Returns a tuple of input args shifted by values from another tuple. */
    template <typename OrigArg1, typename ...OrigArgs, typename ArgType1, typename ...ArgTypes, 
        typename std::enable_if<sizeof...(OrigArgs)==sizeof...(ArgTypes), int>::type = 0,  
        typename std::enable_if<sizeof...(OrigArgs)!=0, int>::type = 0 > 
        std::tuple<OrigArg1, OrigArgs...> _shiftArgs(const std::tuple<OrigArg1, OrigArgs...>&in, const std::tuple<ArgType1, ArgTypes...>& shift_args) const {
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-sizeof...(ArgTypes)-1>(grids_).shift(arg1,shift_arg1); 
    return std::tuple_cat(std::forward_as_tuple(out1),this->_shiftArgs(tuple_tail(in),tuple_tail(shift_args)));
}

    /** Specialization of _shiftArgs for a tuple of 1 element. */
    template <typename OrigArg1, typename ArgType1>
        std::tuple<OrigArg1> _shiftArgs(const std::tuple<OrigArg1>&in, const std::tuple<ArgType1>& shift_args) const;

};

} // end of namespace gftools
#endif // endif::ifndef ___GFTOOLS_GRID_OBJECT_H___

#include "grid_object.hxx"

