#ifndef _RCAC_V4_HPP_
#define _RCAC_V4_HPP_

// #include "Eigen/Dense"
// #include "Eigen/Core"
#include "matrix/math.hpp"
// #include "Eigen/Sparse"
// #include <iostream>

using namespace matrix;
using namespace std;
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
    RCAC(double, double, int, matrix::Matrix<float, 1, 2>);

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
    matrix::Matrix<float, 3, 1> get_thetak()
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

    matrix::Matrix<float, 1, 3> getPhi()
    {
        return Phi_k;
    };

    // void buildRegressor(Eigen::VectorXd &zIn);

    void init_RCAC();
    void set_RCAC_data(double, double);
    void buildRegressor(double zkm1, double zkm1_int, double zkm1_diff);
    void filter_data();
    void update_theta();
    float compute_uk(double, double, double, double);

protected:
    //RCAC Working variables
    matrix::Matrix<float, 3, 3> P;
    matrix::Matrix<float, 3, 1> theta;

    float u_k, u_km1, u_filt;
    float z_km1, z_filt;
    matrix::Matrix<float, 1, 3> Phi_k, Phi_filt;

    matrix::Matrix<float, 2, 1> ubar;        // Size nf by 1
    matrix::Matrix<float, 3, 3> Phibar;      // Size nf+1 by 1
    matrix::Matrix<float, 2, 3> PhibarBlock; // Size nf by 1

    float Gamma;
    float Idty_lz;
    matrix::Matrix<float, 1, 1> one_matrix;
    matrix::Matrix<float, 1, 1> dummy;

    int kk = 1;
};

RCAC::RCAC(double P0_val, double lambda_val, int nf_val, matrix::Matrix<float, 1, 2> filtNu_val)
{
    P0 = P0_val;
    lambda = lambda_val;
    nf = nf_val;
    filtNu = filtNu_val;
}

void RCAC::init_RCAC()
{
    P = eye<float, 3>() * P0;
    theta.setZero();

    u_km1 = 0;
    z_km1 = 0;

    Phi_k.setZero();

    ubar.setZero();
    Phibar.setZero();

    Phi_filt.setZero();

    one_matrix = eye<float, 1>();
}

void RCAC::set_RCAC_data(double zkm1, double ukm1)
{
    z_km1 = zkm1;
    u_km1 = ukm1;
}

void RCAC::buildRegressor(double z, double z_int, double z_diff)
{
    Phi_k(0, 0) = z;
    Phi_k(0, 1) = z_int;
    Phi_k(0, 2) = z_diff;
}

void RCAC::filter_data()
{
    for (int ii = nf - 1; ii > 0; ii--)
    {
        ubar(ii, 0) = ubar(ii - 1, 0);
    }
    ubar(0, 0) = u_km1;

    for (int ii = nf; ii > 0; ii--)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            Phibar(ii, jj) = Phibar(ii - 1, jj);
        }
    }
    for (int jj = 0; jj < 3; jj++)
    {
        Phibar(0, jj) = Phi_k(0, jj);
    }

    // UbarBlock = ubar.block(0, 0, nf - 1, 1);
    // ubar.block(1, 0, nf - 1, 1) = UbarBlock;
    // ubar(0, 0) = u_km1(0, 0);

    // PhibarBlock = Phibar.block(0, 0, nf, 3);
    // Phibar.block(1, 0, nf, 3) = PhibarBlock;
    // Phibar.block(0, 0, 1, 3) = Phi_k;

    z_filt = z_km1;

    dummy = filtNu * ubar;
    u_filt = dummy(0, 0);

    for (int ii = 0; ii < nf; ii++)
    {
        for (int jj = 0; jj < 3; jj++)
        {
            PhibarBlock(ii, jj) = Phibar(ii + 1, jj);
        }
    }
    Phi_filt = filtNu * PhibarBlock;
}

void RCAC::update_theta()
{
    if (kk > 3)
    {
        dummy = Phi_filt * P * Phi_filt.transpose();
        Gamma = lambda + dummy(0, 0);
        P = P - P * Phi_filt.transpose() * 1 / Gamma * Phi_filt * P;
        P = P / lambda;

        theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix + Phi_filt * theta - u_filt);
    }
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << "\t" << theta.transpose() << endl;
}

float RCAC::compute_uk(double z, double z_int, double z_diff, double u)
{
    set_RCAC_data(z, u);
    buildRegressor(z, z_int, z_diff);
    filter_data();
    update_theta();
    dummy = Phi_k * theta;
    u_k = dummy(0, 0);
    // shift_data();
    kk = kk + 1;
    return u_k;
}

#endif