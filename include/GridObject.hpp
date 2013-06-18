#pragma once
#ifndef ___GFTOOLS_GRID_OBJECT_HPP___
#define ___GFTOOLS_GRID_OBJECT_HPP___

#include "Defaults.hpp"
#include "Tools.hpp"
//#include "Grid.h"
#include "Container.hpp"

namespace GFTools {

/** A GridObject is a wrapper over Container class, that stores data,
 * defined on multiple different grids.
 */
template< typename ValueType, typename ...GridTypes> 
class GridObject 
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
    typedef ValueType ValType;
    static constexpr size_t N = sizeof...(GridTypes);
    typedef std::array<size_t, N> PointIndices;
    class exPointMismatch : public std::exception { virtual const char* what() const throw() { return "Index mismatch."; }; };
protected:
    /** Grids on which the data is defined. */
    mutable std::tuple<GridTypes...> _grids;
public:
    /** The dimensions of the Container - deduced from grids. */
    PointIndices _dims;
protected:
    /** A pointer to the Container. A pointer is used as there exist no default 
     * constructor for the Container.
     */
    //std::unique_ptr<Container<ValueType, N>> _data;
    mutable Container<ValueType, N> _data;

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
    template <size_t Nc, typename CT, typename ArgType1, typename ...ArgTypes> struct ContainerExtractor { 
        /** Gets the data by values. */
        static ValueType get(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
        static ValueType& get_ref(CT &data, const std::tuple<GridTypes...> &grids, const ArgType1& arg1, const ArgTypes&... args);
             };
    /** Specialization of ContainerExtractor for 1-dim container. */
    template <typename CT, typename ArgType1> struct ContainerExtractor<1,CT,ArgType1> {
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
    GridObject( const std::tuple<GridTypes...> &grids);
    /** Constructor of grids and data. */
    GridObject( const std::tuple<GridTypes...> &grids, const Container<ValueType, sizeof...(GridTypes)>& data);
    /** Copy constructor. */
    GridObject( const GridObject<ValueType, GridTypes...>& rhs);
    /** Move constructor. */
    GridObject( GridObject<ValueType, GridTypes...>&& rhs);

    const std::tuple<GridTypes...> getGrids() const;
    const size_t getTotalContainerSize() const;
    /** Returns an Mth grid in _grids. */
    template<size_t M> auto getGrid() const -> const decltype(std::get<M>(_grids));
    /** Returns the top level grid. */
    auto getGrid() const -> const decltype(std::get<0>(_grids));
    /** Returns element number i, which corresponds to (*_grid)[i]. */
    auto operator[](size_t i)->decltype(_data[i]);
    /** Const operator[]. */
    auto operator[](size_t i) const ->decltype(_data[i]) const;
    //template <size_t M> ValueType& operator[](const std::array<size_t,M>& in);
    /** Returns the _data Container. */
    Container<ValueType, sizeof...(GridTypes)>& getData(){return _data;};
    /** Fills the Container with a provided function. */
    template <typename ...ArgTypes> void fill(const std::function<ValueType(ArgTypes...)> &);
    void fill(const FunctionType &in);
    void fill(const PointFunctionType &in);
    //template <typename ...ArgTypes> void fill_tuple(const std::function<ValueType(const std::tuple<ArgTypes...>)> &);
    void fill_tuple(const std::function<ValueType(ArgTupleType)>& in);
    void fill_tuple(const std::function<ValueType(PointTupleType)>& in);
    /** Fills the Container with any proper class with call operator. Untested */
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
    //template <typename ...ArgTypes> GridObject& operator= (const std::function<ValueType(ArgTypes...)> &);
    /** Same as operator=, but allows for non-equal grids. Slow. Uses analytic function to provide missing values. */
    GridObject& copyInterpolate(const GridObject &rhs);
    /** Algebraic operators. */
    GridObject& operator= (const GridObject & rhs);
    GridObject& operator= (const ValueType & rhs);
    GridObject& operator*= (const GridObject & rhs);
    GridObject& operator*= (const ValueType& rhs);
    GridObject operator* (const GridObject & rhs) const;
    GridObject operator* (const ValueType & rhs) const;
    GridObject& operator+= (const GridObject & rhs);
    GridObject& operator+= (const ValueType& rhs);
    GridObject operator+ (const GridObject & rhs) const;
    GridObject operator+ (const ValueType & rhs) const;
    GridObject& operator-= (const GridObject & rhs);
    GridObject& operator-= (const ValueType& rhs);
    GridObject operator- (const GridObject & rhs) const;
    GridObject operator- (const ValueType & rhs) const;
    GridObject& operator/= (const GridObject & rhs);
    GridObject& operator/= (const ValueType& rhs);
    GridObject operator/ (const GridObject & rhs) const;
    GridObject operator/ (const ValueType & rhs) const;
    friend GridObject operator* (const ValueType & lhs, const GridObject & rhs) {return rhs*lhs;};
    friend GridObject operator+ (const ValueType & lhs, const GridObject & rhs) {return rhs+lhs;};
    friend GridObject operator- (const ValueType & lhs, const GridObject & rhs) {return rhs*(-1.0)+lhs;};
    friend GridObject operator/ (const ValueType & lhs, const GridObject & rhs) {GridObject out(rhs); out=lhs; return out/rhs;};

    /** Returns the complex conjugate of this object, if it's complex valued. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type=0>
        GridObject conj();
    /** Returns a norm of difference between two objects. */
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, ComplexType>::value, int>::type = 0>
    RealType diff(const GridObject &rhs) const; 
    template <typename U = ValueType, typename std::enable_if<std::is_same<U, RealType>::value, int>::type = 0>
    RealType diff(const GridObject &rhs) const;
    /** Returns the sum of all elements in the container. */
    ValueType sum();

    /** Returns an object with arguments, shifted by the given values.
     * \param[in] args A pack of arguments to shift the object
     */
    template <typename ...ArgTypes> GridObject shift(ArgTypes... args) const;
    template <typename ...ArgTypes> GridObject shift(const std::tuple<ArgTypes...>& arg_tuple) const;
    /** Save the data to the txt file. */
    template<typename U = typename std::enable_if<sizeof...(GridTypes)==1>::type> 
    void savetxt(const std::string& fname) const;
    /** Loads the data to the txt file. */
    template<typename U = typename std::enable_if<sizeof...(GridTypes)==1>::type> 
    void loadtxt(const std::string& fname);
    /** Dumps the object to the stream. */
    template <typename ValType, class ...GridTypes2> friend std::ostream& operator<<(std::ostream& lhs, const GridObject<ValType,GridTypes2...> &in);
    
    class exIOProblem : public std::exception { virtual const char* what() const throw(){return "IO problem.";} }; 
    class exWrongIndex : public std::exception { virtual const char* what() const throw(){return "Index out of bounds";}}; 
    
/** Returns a tuple of input args shifted by values from another tuple. */
    template <typename OrigArg1, typename ...OrigArgs, typename ArgType1, typename ...ArgTypes, 
        typename std::enable_if<sizeof...(OrigArgs)==sizeof...(ArgTypes), int>::type = 0,  
        typename std::enable_if<sizeof...(OrigArgs)!=0, int>::type = 0 > 
        std::tuple<OrigArg1, OrigArgs...> _shiftArgs(const std::tuple<OrigArg1, OrigArgs...>&in, const std::tuple<ArgType1, ArgTypes...>& shift_args) const {
    OrigArg1 arg1 = std::get<0>(in);
    ArgType1 shift_arg1 = std::get<0>(shift_args);
    OrigArg1 out1 = std::get<sizeof...(GridTypes)-sizeof...(ArgTypes)-1>(_grids).shift(arg1,shift_arg1); 
    return std::tuple_cat(std::forward_as_tuple(out1),this->_shiftArgs(__tuple_tail(in),__tuple_tail(shift_args)));
}

    /** Specialization of _shiftArgs for a tuple of 1 element. */
    template <typename OrigArg1, typename ArgType1>
        std::tuple<OrigArg1> _shiftArgs(const std::tuple<OrigArg1>&in, const std::tuple<ArgType1>& shift_args) const;

};

} // end of namespace GFTools
#endif // endif::ifndef ___GFTOOLS_GRID_OBJECT_H___

#include "GridObject.hxx"

