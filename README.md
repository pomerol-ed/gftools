This is a set of tools to work in condmat problems. It consists of three basic objects:

- Grid<ValueType>. A class to represent a set of discrete points of type ValueType and to perform operations with objects defined on the
manifold of ValueType values. ValueType is normally a number (int, double, complex<double>)
There exists a set of implementation of grids for common purposes:
    - MatsubaraGrid< b >. A grid of Matsubara frequencies. "b" corresponds to a choice of fermionic(1) or bosonic(0) grids.
    - RealGrid. A grid of floats, e.g. frequencies. 
    - KMesh. An equidistant set of points from 0 to 2pi.
- Container<N, ValueType>. A recursive container of depth <N> to store ValueType objects.
- GridObject<Grids..., ValueType>. An multidimensional object that is defined on grids Grids... that contains ValueType values.
Example : Green's function.

The code is provided as a header-only library with a set of examples. Requires C++11 and >=Eigen3.1 library.
