
/* Here's a small example of the code, written with usage of some basic
 * objects of the library. It covers initialization of the grids and objects,
 * defined on the grid. Here's what's done:
 */

// This includes the grid of real values.
#include "RealGrid.hpp"
// This includes the definition of a container class - a tool to store data.
#include "Container.hpp"
// This include the GridObject class - an object to store data, defined on a grid.
#include "GridObject.hpp"

/* Some more basic definitions, that are used:
 * RealType is double
 * ComplexType is complex<double>
 * INFO( something ) is the same as std::cout << something << std::endl;
 * DEBUG( something ) is the same as INFO, but it also prints the number of the
 * line in the code, where it is called from.
 * ERROR( something ) is the same as INFO but is used for outputting the errors,
 * it sends the text to the error stream (stderr).
 */

// This is a shorthand for not writing GFTools::RealGrid, GFTools::GridObject etc
using namespace GFTools;

int main()
{
    // Let's start with grids.
    // This generates a uniform grid of 20 real numbers from -5 to 5
    RealGrid grid1(-5,5,20);
    // And prints it
    INFO(grid1);
    
    // New c++ feature - store a function
    std::function<RealType(int)> f1;
    /* The right hand side of next line is called a "lambda" construction.
     * We construct some function of integer and assign it to f1. */
    f1 = [](int n){return std::pow(3,n);};
    // f1 is now 3^n. Let's now generate a non-uniform grid and print it.
    RealGrid grid2(0,10,f1);
    INFO(grid2)
    
    // Another option is to generate a set of values and generate a grid.
    std::vector<RealType> v1;
    for (RealType a=-3;a<=3;a+=0.01) v1.push_back(a);
    RealGrid grid3(v1);
    INFO(grid3);
    
    // Now let's make an integration over a grid
    // Let's create another function
    std::function<RealType(RealType)> sinF;
    sinF = [](RealType x){return sin(x);};
    // sinF is now sin(x). We integrate an even function and generate zero.
    INFO(grid1.integrate(sinF));
    
    // Ok, that's all for grids.
    // Now, let's define objects on a grid.
    // First, make a shorthand
    typedef GridObject<RealType,RealGrid> GFRetarded;
    // Now in order to construct an object on a grid, we need to provide a grid.
    // Let's use the one, already created.
    GFRetarded G1(grid1);
    // After construction, GridObject stores only zeros.
    // Let's fill it with something more interesting, say sin(x).
    G1.fill(sinF);
    INFO(G1)
    // And now let's integrate it. It should give zero once again.
    INFO(grid1.integrate(G1));
    
    // In order to have an access to the data points, one can just use [] operator.
    // Notice the values are counted from zero.
    INFO(G1[2]);
    
    //GridObjects also support mathematical operations. Let's create another
    //GridObject first.
    GFRetarded G2(grid1);
    std::function<RealType(RealType)> cosF;
    cosF = [](RealType x){return cos(x);};
    G2.fill(cosF);
    INFO(G2);
    // Now let's make some operations. sin^2(x) + cos^2(x) is trivial.
    GFRetarded G3 = G2*G2+G1*G1;
    INFO(G3);
    
    // Some advanced features. Grids allow to interpolate values.
    // ERROR, this doesn't work well now.
    INFO(grid1.getValue(G1,4.1));
    
    /*GridObjects allow for the extrapolation of the values, if they are out 
     of grid bounds. For this we need to assign the extrapolation
     function. */
    G1._f = sinF;
    INFO(G1(2*PI));
}
