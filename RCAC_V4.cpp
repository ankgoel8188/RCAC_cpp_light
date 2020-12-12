#include "RCAC_V4.hpp"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Sparse"

// using namespace std;

RCAC::RCAC(double P0_val, double N1_val, double lambda_val)
{
    P0 = P0_val;
    lambda = lambda_val;
    N1 = N1_val;
}

void RCAC::init_RCAC()
{
    int lu = 1;
    int lz = 1;
    P.setIdentity(3 * lu, 3 * lu);
    P = P0 * P;
    // cout << kk << "\n"
    //      << P << endl;
    theta.setZero(3 * lu);

    u_k.setZero(lu);
    u_km1.setZero(lu);

    z_k.setZero(lz);
    z_km1.setZero(lz);

    Phi_k = Eigen::MatrixXd::Zero(lu, 3 * lu);
    Phi_km1 = Phi_k;
    Phi_km2 = Phi_k;
    Idty_lz.setIdentity(lz, lz);
}

void RCAC::buildRegressor(double zk, double zk_int, double zk_diff)
{
    Phi_k(0, 0) = zk;
    Phi_k(0, 1) = zk_int;
    Phi_k(0, 2) = zk_diff;
}

void RCAC::set_performance(double zk)
{
    z_k(0, 0) = zk;
}

void RCAC::update_theta()
{
    z_filt = z_k;
    u_filt = N1*u_km2;
    Phi_filt = N1*Phi_km2;

    if (kk > 3)
    {
        Gamma = lambda * Idty_lz + Phi_filt * P * Phi_filt.transpose() * N1;
        P = P - P * Phi_filt.transpose() * Gamma.inverse() * Phi_filt * P;
        P = P / lambda;
        theta = theta - P * Phi_filt.transpose() * (z_filt + Phi_filt * theta - u_filt);
    }
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << "\t" << theta.transpose() << endl;   
}

Eigen::VectorXd RCAC::compute_uk(double zk, double zk_int, double zk_diff)
{
    buildRegressor(zk, zk_int, zk_diff);
    set_performance(zk);
    update_theta();
    u_k = Phi_k * theta;
    shift_data();
    return u_k;
}

void RCAC::shift_data()
{
    u_km2 = u_km1;
    u_km1 = u_k;
    Phi_km2 = Phi_km1;
    Phi_km1 = Phi_k;
    kk = kk + 1;
}