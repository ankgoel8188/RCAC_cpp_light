#ifndef _RCAC_V4_HPP_
#define _RCAC_V4_HPP_

// #include "Eigen/Dense"
// #include "Eigen/Core"
#include "matrix/math.hpp"
// #include "Eigen/Sparse"
// #include <iostream>

using namespace matrix;

/**
 * The parent RCAC class. This class handles all the low level computation of RCAC
 * such as the filtering, coefficient updates, and keeping track of the regressors.
 *
 * Almost all the methods are polymorphic and can be modified by child classes
 * to create RCAC algorithms with more complex filtering.
 */
class RCAC
{
    float P0;
    float lambda;
    int nf;
    matrix::Matrix<float, 1, 2> filtNu;

public:

    // RCAC();
    RCAC( double, double, int, matrix::Matrix<float, 1,2>);

    /**
         * Returns RCAC's computed value for the control. Must run oneStep at least once.
         */
    //Function getControl: Get the computed control input
    float get_uk()
    {
        return u_k;
    };

    /**
         * Returns a vector of the current RCAC coefficients.
         */
    //Function getCoeff: Get the RCAC coefficients
    matrix::Matrix<float, 3,1> get_thetak()
    {
        return theta;
    };


    /**
         * Returns the timestep of RCAC
         */
    //Function getkk: Get timestep
    int getkk()
    {
        return kk;
    };

    double getP()
    {
        return P(0, 0);
    };

    matrix::Matrix<float, 1,3> getPhi()
    {
        return Phi_k;
    };

    // void buildRegressor(Eigen::VectorXd &zIn);

    void init_RCAC();
    void set_RCAC_data(double , double);
    void buildRegressor(double zkm1, double zkm1_int, double zkm1_diff);
    void filter_data();
    void update_theta();
    float compute_uk(double, double, double, double);

protected:
    //RCAC Working variables
    matrix::Matrix<float, 3,3> P;
    matrix::Matrix<float, 3,1> theta;

    float u_k, u_km1, u_filt;
    float z_km1, z_filt;
    matrix::Matrix<float, 1,3> Phi_k, Phi_filt;

    matrix::Matrix<float, 2,1> ubar;        // Size nf by 1
    matrix::Matrix<float, 3,3> Phibar;      // Size nf+1 by 1
    matrix::Matrix<float, 2,3> PhibarBlock;      // Size nf by 1

    float Gamma;
    float Idty_lz;
    matrix::Matrix<float, 1,1> one_matrix;
    matrix::Matrix<float, 1,1> dummy;

    int kk = 1;
};


#endif
