#ifndef _RCAC_HPP_
#define _RCAC_HPP_

#include <iostream>
#include "RCAC_V4.hpp"
// #include "Eigen/Core"
using namespace std;
// using namespace Eigen;
int main()
{
    
    matrix::Matrix<float, 1, 2> filtNu;
    filtNu(0,0) = 0;
    filtNu(0,1) = 1;

    RCAC myRCAC(5.0, 1.0, 2, filtNu);
    myRCAC.init_RCAC();
    myRCAC.buildRegressor(0,0,0);
    myRCAC.set_RCAC_data(1,2);
    myRCAC.filter_data();
    myRCAC.update_theta();
    myRCAC.compute_uk(0,0,0,0);
    // uout = myRCAC.compute_uk(z1, g1, zd1,u(0,0));
    std::cout << filtNu(0,0);
}
#endif