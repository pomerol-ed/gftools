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
/*
    std::array<size_t,3> sizes = { { 3, 3, 3 } };
    std::array<size_t,2> size2 = { { 3 , 3 } };
    std::array<size_t,1> size1 = { { 3 } };
    typedef boost::multi_array<double,3> array;
    typedef array::size_type size_type;
    size_type num_elements = 27;
    std::vector<double> vals(num_elements, 4.5);

    boost::multi_array<double, 3> A1(sizes);
    A1.assign(vals.begin(),vals.end());

    auto D1 = A1[1][1];
    array::subarray<1>::type D2 = D1;
    A1[1][1][0]=3.3;
    DEBUG(D2[0]);

    //boost::multi_array_ref<double, 1> array2(&(*D2.begin()), size1);
    boost::multi_array_ref<double, 1> array2(A1[0].origin(), size2);

    typedef array::subarray<2>::type subarray;
    subarray B1 = A1[1];
    subarray::value_type C1 = B1[0];

    DEBUG(equal(A1[1][0],C1));
    DEBUG(equal(B1[0],C1));

*/
    std::cout << "Hi!" << std::endl;
    std::array<size_t,2> Ar1 {{1,3}};
    std::array<size_t,3> Ar2 {{2,1,3}};


    INFO("Container")

    Container<ComplexType,2> B(Ar1);
    Container<ComplexType,2> C(Ar1);
    Container<ComplexType,3> D(Ar2);
    
    B[0][2]=3.0;
    C[0][2]=-2.0;
    DEBUG(B[0][2]);
    DEBUG(B[0][0]);

    DEBUG(B._data.shape()[0]);
    DEBUG(B[0]._data.shape()[0]);

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
    C+=B;
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
    decltype(Vals) Vals_2(Vals.getAsMatrix());
    INFO(Vals_2);

    decltype(Vals) Vals_3(dim1);
    Vals_3 = Vals.getAsMatrix();
    INFO(Vals_3);

    //if ((Vals-Vals_2).sum()!=0 || (Vals-Vals_3).sum()!=0) return EXIT_FAILURE;
/*
    INFO("1-dim test");
    Container<RealType,1> Vals_1d(5);
    Vals_1d[2]=1.6;
    Vals_1d[4]=3;
    INFO(Vals_1d);
    INFO(Vals_1d.getAsVector());
    INFO(Vals_1d.getAsDiagonalMatrix());
  */  

    return EXIT_SUCCESS;
}
