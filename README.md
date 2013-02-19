This is a set of tools to work with numeric condmat problems. It consists of three basic objects:

- Grid<ValueType>. A class that represents a discrete set of points of type ValueType and allows to perform operations with objects defined on the
manifold of ValueType values. ValueType is normally a number (int, double, complex<double>)
There exists a set of implementation of grids for common purposes:
    - MatsubaraGrid< b >. A grid of Matsubara frequencies. "b" corresponds to a choice of fermionic(1) or bosonic(0) grids.
    - RealGrid. A grid of floats, e.g. frequencies. 
    - KMesh. An equidistant set of points from 0 to 2pi.
- Container<N, ValueType>. A recursive container of depth <N> to store ValueType objects.
- GridObject<Grids..., ValueType>. An multidimensional object that is defined on a set of grids, Grids... that contains ValueType values.
Example :
    - GridObject<MatsubaraGrid<1>,ComplexType> - Green's function in imaginary frequencies
    - GridObject<MatsubaraGrid<1>,KMesh,ComplexType> - Green's function in imaginary frequencies and reciprocal 1d space.
    - GridObject<RealGrid,ComplexType> - Retarded or Advanced Green's function, etc...

The code is provided as a header-only library with a set of examples. Requires C++11 and >=Eigen3.1 library.
