#ifndef _RCAC_HPP_
#define _RCAC_HPP_

#include <iostream>
#include "RCAC.h"
#include "matrix/math.hpp"
using namespace std;

int main()
{
    //**************Simulation Parameters***************
    //Initialize the plant
    // int lx = 2;
    // int lu = 1;
    // int ly = 1;
    // int lz = ly;
    matrix::Matrix<float, 2, 2> A;
    matrix::Matrix<float, 2, 1> B;
    matrix::Matrix<float, 1, 2> C;
    A(0,0) = 1.7;
    A(0,1) = -0.72;
    A(1,0) = 1.0;
    A(1,1) = 0.0;
    
    B(0,0) = 1.0;
    B(1,0) = 0.0;

    C(0,0) = 1.0;
    C(0,1) = -0.85;
    

    matrix::Matrix<float, 2, 1> x;
    matrix::Matrix<float, 1, 1> y;
    matrix::Matrix<float, 1, 1> z;
    matrix::Matrix<float, 1, 1> u;
    
    x(0,0) = 1.0;
    x(1,0) = 1.0;
    //End time of the simulation
    int kend = 18;
    //*************************************************
    
    // RCAC RCAC2;
    RCAC myRCAC;
    // myRCAC.set_RCAC_parameters(5,1,2,1);
    myRCAC.init_RCAC(5,1,2);

    // u = MatrixXd::Zero(lu, 1);

    double z1, z1km1, g1, zd1;
    float uout;
    z1km1 = 0;
    g1 = 0;


    for (int k = 0; k < kend; k++)
    {
        y = C * x;
        z = y - 1;
        x = A * x + B * u;
        z1 = z(0, 0);
        g1 = g1 + z(0, 0);
        zd1 = z1 - z1km1;
        z1km1 = z1;
        uout = myRCAC.compute_uk(z1, g1, zd1, u(0,0));
        u(0, 0) = uout;
        cout << myRCAC.getkk() << "\t" << z(0, 0) << endl;
        // std::cout << " y: " << y(0) << ", "  << y(1) << "\n z: " << z(0) << ", "  << z(1) << "\n u: " << u(0) << ", "  << u(1) << std::endl;
    }
    // ****************************************************
}

#endif