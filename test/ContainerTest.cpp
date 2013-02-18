#include <numeric>

#include "Defaults.hpp"
#include "Container.hpp"
#include "Logger.h"

#include <iostream>
#include <ctime>
#include <array>

using namespace GFTools;

int main()
{
    Log.setDebugging(true);
    std::cout << "Hi!" << std::endl;

    std::array<size_t,2> Ar1 {{1,3}};
    std::array<size_t,3> Ar2 {{1,1,3}};
    Container<2,ComplexType> B(Ar1);
    Container<2,ComplexType> C(Ar1);
    Container<3,ComplexType> D(Ar2);
    Container<3,ComplexType> E(D);
    B[0][2]=3.0;
    C[0][2]=-2.0;
    DEBUG(B+C*2);
    DEBUG(B[0]+C[0]);
    DEBUG((B+C)[0]);
    DEBUG((B[0]*3.0));
    DEBUG((B-C)[0]);
    D[0][0][0]=-1.0;
    D[0][0][1]=1.0;
    E*=(-1);
    DEBUG(D);
    DEBUG(E);
    DEBUG(D+E);
    DEBUG(D*5+2.0);
    DEBUG(D[0]+B);

    INFO("Matrix test");
    INFO("2-dim");
    std::array<size_t,2> dim1 {{4,5}};
    Container<2,RealType> Vals (dim1);
    Vals[1][0] = 1.0;
    Vals[0][3] = -0.4;
    Vals[2][1] = 1.5;
    INFO(Vals);
    INFO(Vals.getAsMatrix());
    decltype(Vals) Vals_2(Vals.getAsMatrix());
    INFO(Vals_2);

    decltype(Vals) Vals_3(dim1);
    Vals_3 = Vals.getAsMatrix();
    INFO(Vals_3);

    if ((Vals-Vals_2).sum()!=0 || (Vals-Vals_3).sum()!=0) return EXIT_FAILURE;
    INFO("1-dim test");
    Container<1,RealType> Vals_1d(5);
    Vals_1d[2]=1.6;
    Vals_1d[4]=3;
    INFO(Vals_1d);
    INFO(Vals_1d.getAsVector());
    INFO(Vals_1d.getAsDiagonalMatrix());
    

    return EXIT_SUCCESS;
}
