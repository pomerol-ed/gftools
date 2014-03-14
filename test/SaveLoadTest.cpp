#include <numeric>

#include "Defaults.hpp"
#include "Tools.hpp"
#include "MatsubaraGrid.hpp"
#include "Container.hpp"
#include "GridObject.hpp"
#include "KMesh.hpp"
#include "EnumerateGrid.hpp"

#include <iostream>
#include <ctime>
#include <array>

using namespace GFTools;

int main()
{
    FMatsubaraGrid fgrid(0,2,10);
    KMesh qgrid(5);

    GridObject<ComplexType,FMatsubaraGrid,KMesh> D2(std::make_tuple(fgrid,qgrid));
    GridObject<ComplexType,FMatsubaraGrid,KMesh> D3(std::make_tuple(fgrid,qgrid));
    
    D2.getData()[0][1]=4.0;
    D2.getData()[1][2]=3.1;
    D2.savetxt("data_in.dat");
    INFO("Saved" << D2);

    try {
    D3.loadtxt("data_in.dat");
        }
    catch (...)
    {
        ERROR("Caught error");
        return EXIT_FAILURE;
    };
    INFO("Loaded" << D2);
    D3.savetxt("data_out.dat");
    if (D3.diff(D2)>1e-5) return EXIT_FAILURE;

    GridObject<ComplexType,FMatsubaraGrid,KMesh,KMesh> D4(std::make_tuple(fgrid,qgrid,qgrid));
    GridObject<ComplexType,FMatsubaraGrid,KMesh,KMesh> D5(std::make_tuple(fgrid,qgrid,qgrid));

    GridObject<ComplexType,FMatsubaraGrid,KMesh,KMesh>::FunctionType f1 = [](ComplexType w, RealType k1, RealType k2){return 1.0/(w + cos(k1)+cos(k2));};
    D4.fill(f1);
    D4.savetxt("data_in2.dat");
    D5.loadtxt("data_in2.dat");
    if (D4.diff(D5)>1e-5) return EXIT_FAILURE;
    D5.savetxt("data_out2.dat");
    
    return EXIT_SUCCESS;
}
