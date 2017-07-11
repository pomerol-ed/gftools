##### Overview
gftools is a set of tools to work with numerical condmat problems. It provides three basic classes of objects:

- grids. 
There exists a set of implementation of grids for common purposes:
    - `matsubara_grid<b>`. A grid of Matsubara frequencies. "b" corresponds to a choice of fermionic(1) or bosonic(0) grids.
    - `real_grid `. A grid of floats, e.g. frequencies. 
    - `kmesh`. An equidistant periodic set of points (for sampling Brillouine zones).
    - `enum_grid`. A grid of integers.
- `container<ValueType, N>`. A recursive multidimensional container of depth <N> to store ValueType objects.
- `grid_object<ValueType, Grids...>`. Green's function. In a broader scope - a multidimensional object that is defined on a set of grids (Grids...) that contains ValueType values.
Examples :
    - `grid_object<complex_type,matsubara_grid<1>>` - Green's function in imaginary frequencies
    - `grid_object<complex_type,matsubara_grid<1>,kmesh>` - Green's function in imaginary frequencies and reciprocal 1d space.
    - `grid_object<complex_type,real_grid>` - Retarded or Advanced Green's function, etc...

Look for examples in "example/" directory for detailed information.

##### Installation ###
The code is is provided as a header-only library with a set of examples and tests.
The `gftools.hpp` in the repo root can be included in any derivative projects.
To compile examples and tests create a build directory and run 

1. `cmake -DExamples=ON -DTesting=ON {path_to_gftools}`
2. `make`
3. `make test` (for running tests)
4. example will be build in example subdirectory
5. `make doc` to generate documentation

##### Dependencies 
- c++11-compatible compiler (tested with clang >= 3.1, gcc >= 4.8.2, icpc >= 14.0.2)  
- *Boost* headers 
- *Eigen* >=3.1 
- *git* to fetch the code 
- *cmake* to build tests and examples (optional)
- *doxygen* for documentation (optional)

##### Extra features
- FFT support via FFTW
- HDF5 support via alpscore (http://www.alpscore.org)

##### Author
- Andrey Antipov, *Andrey.E.Antipov[at]gmail.com*, 2013-now.

##### Distribution
Open-source under GPLv2 license.

##### Academic usage
Please attribute this work by a citation to arXiv:1705.00024 (2017).
