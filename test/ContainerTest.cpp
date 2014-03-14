#include <numeric>

#include "Defaults.hpp"
#include "Container.hpp"

#include <iostream>
#include <ctime>
#include <array>

#include <boost/multi_array.hpp>

bool equal(const double& a, const double& b)
{
  return a == b;
}

template <typename ArrayA, typename ArrayB>
bool equal(const ArrayA& A, const ArrayB& B)
{
  typename ArrayA::const_iterator ia;
  typename ArrayB::const_iterator ib = B.begin();
  for (ia = A.begin(); ia != A.end(); ++ia, ++ib)
    if (!equal(*ia, *ib))
      return false;
  return true;
}

using namespace GFTools;

int main()
{
    std::cout << "Hi!" << std::endl;
    std::array<size_t,1> Ar1 {{5}};
    std::array<size_t,2> Ar2 {{3,3}};
    std::array<size_t,3> Ar3 {{2,3,3}};


    INFO("Container")

    Container<ComplexType,2> A1(Ar2), A2(Ar2), A3(Ar2);
    A1 = 1.;
    A2 = 2.;
    DEBUG(A1);
    A1+=A2;
    DEBUG(A1);
    A3 = 3.0*(A1+0.5)*2+A2;
    DEBUG(A3);
    DEBUG(A1);
    auto a01 = 3*(A1[1]+A1[0])*2;
    DEBUG(a01);
    DEBUG(A1);

    Container<ComplexType,1> E1(A1[1]), E2(A1[2]);
    E1*=3;
    DEBUG(E1);
    DEBUG(A1);
    

    Container<ComplexType,2> B(Ar2);
    Container<ComplexType,2> C(Ar2);
    Container<ComplexType,3> D(Ar3);
    
    B[0][2]=3.0;
    C[0][2]=-2.0;

    std::array<size_t, 2> coord1 {{0,1}};
    DEBUG(B._data(coord1));

    
    Container<ComplexType,2>::iterator t1;
    auto it2 = B.begin();
    auto it3 = B[0].begin();
    
    DEBUG(B[0][2]);
    DEBUG(B[0][0]);

    //DEBUG(B._data.shape()[0]);
    //DEBUG(B[0]._data.shape()[0]);

    DEBUG(B);
    DEBUG(B[0]);

    D[0][0][0]=-1.0;
    D[0][0][1]=1.0;
    DEBUG(D);
    DEBUG(D._data.num_elements());

    DEBUG(C);
    C=B;
    DEBUG(C);

    C[0][0]=1.0;
    C[0][2]=-5;
    DEBUG(C);
    C+=B;
    DEBUG(C);
    C+=1.0;
    DEBUG(C);
    DEBUG(C+B);
    Container<ComplexType,3> E(D);

    DEBUG(B);
    B+=2.0;
    DEBUG(B);
    DEBUG(B[0]+C[0]);
    DEBUG((B+C)[0]);
    DEBUG((B-C)[0]);

    B*=4.;
    DEBUG(B);
    DEBUG(B[0]*2.0);
    DEBUG(D+2.0);
    DEBUG(D[0]+B);
    //DEBUG((B[0]*3.0));
    E*=(-1);
    DEBUG(D);
    DEBUG(E);
    DEBUG(D+E);
    DEBUG(D*5+2.0);


    

    INFO("Matrix test");
    INFO("2-dim");
    std::array<size_t,2> dim1 {{4,5}};
    Container<RealType,2> Vals (dim1);
    Vals[1][0] = 1.0;
    Vals[0][3] = -0.4;
    Vals[2][1] = 1.5;
    INFO(Vals);
    INFO(Vals.getAsMatrix());
    DEBUG("!");
    Container<RealType,2> Vals_2(Vals.getAsMatrix());
    INFO(Vals_2);

    decltype(Vals) Vals_3(dim1);
    Vals_3 = Vals.getAsMatrix();
    INFO(Vals_3);

    if ((Vals-Vals_2).sum()!=0 || (Vals-Vals_3).sum()!=0) return EXIT_FAILURE;
    INFO("1-dim test");
    Container<RealType,1> Vals_1d(std::array<size_t,1>({{5}}));
    Vals_1d[2]=1.6;
    Vals_1d[4]=3;
    INFO(Vals_1d);
    INFO(Vals_1d.getAsVector());
    INFO(Vals_1d.getAsDiagonalMatrix());

    INFO("View test");
    Container<RealType, 3> C1(std::array<size_t,3>({{2,1,3}}));
    C1[0][0][1] = 1.0; C1[0][0][2] = 2.0; C1[1][0][0] = 3.0; C1[1][0][2] = 4.0;
    INFO_NONEWLINE("Data in memory: "); for (size_t i = 0; i<C1._data.num_elements(); ++i) { INFO_NONEWLINE(*(C1._data.origin()+i) << " "); }; INFO("");
    INFO("Container in : " << C1);
    typedef boost::multi_array<RealType, 3>::array_view<3>::type myview_3_t;
    typedef boost::multi_array_types::index_range range;
    myview_3_t boost_view1 = C1._data[boost::indices[range(0,1)][range(0,1)][range(0,2)]]; 
    ContainerView<RealType, 3, 3> C1_v2(myview_3_t(C1._data[boost::indices[range().stride(2)][range()][range().stride(1)]]));
    ContainerView<RealType, 3, 3> C1_v3(C1_v2);
    INFO(C1_v2);
    DEBUG(C1_v2*C1_v3);
    DEBUG(C1);
    C1_v2 *= C1_v3;
    DEBUG(C1_v2);
    DEBUG(C1);
    ContainerView<RealType, 3> C1_v(C1._data[boost::indices[range()][range()][range().stride(3)]]);
    INFO("Container view : " << C1_v);

    return EXIT_SUCCESS;
}
