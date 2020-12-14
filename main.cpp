#ifndef _RCAC_HPP_
#define _RCAC_HPP_

#include <iostream>
#include "RCAC_V4.hpp"
#include "Eigen/Core"
using namespace std;
using namespace Eigen;
int main()
{
    //**************Simulation Parameters***************
    //Initialize the plant
    int lx = 2;
    int lu = 1;
    int ly = 1;
    int lz = ly;
    MatrixXd A(lx, lx);
    A << 1.7, -0.72,
        1.0, 0.0;
    MatrixXd B(lx, lu);
    B << 1.0,
        0.0;
    MatrixXd C(ly, lx);
    C << 1.0, -.85;
    MatrixXd x(lx, 1);
    x << 1.0,
        1.0;
    //End time of the simulation
    int kend = 27;
    //*************************************************
    matrix::Matrix<float, 1, 2> filtNu;
    filtNu(0,0) = 0;
    filtNu(0,1) = 1;

    RCAC myRCAC(5.0, 1.0, 2, filtNu);

    // RCAC myRCAC(5.0, 1.0, nf, filtNu);
    myRCAC.init_RCAC();

    //Run Simulation
    VectorXd y(ly, 1);
    VectorXd z(lz, 1);
    VectorXd u(lu, 1);
    u = MatrixXd::Zero(lu, 1);

    double z1, z1km1, g1, zd1;
    float uout;
    z1km1 = 0;
    g1 = 0;


    for (int k = 0; k < kend; k++)
    {
        y = C * x;
        z = y - VectorXd::Ones(ly);
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