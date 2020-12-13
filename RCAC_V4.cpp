#include "RCAC_V4.hpp"
#include <iostream>

using namespace std;

RCAC::RCAC(double P0_val, double lambda_val, int nf_val, Eigen::MatrixXd filtNu_val)
{
    P0 = P0_val;
    lambda = lambda_val;
    nf = nf_val;
    filtNu = filtNu_val;
}

void RCAC::init_RCAC()
{
    int lu = 1;
    int lz = 1;
    P.setIdentity(3 * lu, 3 * lu);
    P = P0 * P;
    theta.setZero(3 * lu);

    u_km1.setZero(lu);

    z_km1.setZero(lz);

    Phi_k = Eigen::MatrixXd::Zero(lu, 3 * lu);

    Idty_lz.setIdentity(lz, lz);

    ubar.setZero(nf, 1);
    Phibar.setZero(nf + 1, 3);

    PhibarBlock.setZero(nf, 3);
    UbarBlock.setZero(nf-1,1);
}

void RCAC::set_RCAC_data(double zkm1, double ukm1)
{
    z_km1(0, 0) = zkm1;
    u_km1(0, 0) = ukm1;
}

void RCAC::buildRegressor(double z, double z_int, double z_diff)
{
    Phi_k(0, 0) = z;
    Phi_k(0, 1) = z_int;
    Phi_k(0, 2) = z_diff;
}

void RCAC::filter_data()
{
    UbarBlock = ubar.block(0, 0, nf - 1, 1);
    ubar.block(1, 0, nf - 1, 1) = UbarBlock;
    ubar(0, 0) = u_km1(0, 0);

    PhibarBlock = Phibar.block(0, 0, nf, 3);
    Phibar.block(1, 0, nf, 3) = PhibarBlock;
    Phibar.block(0, 0, 1, 3) = Phi_k;

    z_filt = z_km1;
    u_filt = filtNu * ubar;
    Phi_filt = filtNu * Phibar.block(1, 0, nf, 3);
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << endl;
}

void RCAC::update_theta()
{
    if (kk > 3)
    {
        Gamma = lambda * Idty_lz + Phi_filt * P * Phi_filt.transpose();
        P = P - P * Phi_filt.transpose() * Gamma.inverse() * Phi_filt * P;
        P = P / lambda;
        theta = theta - P * Phi_filt.transpose() * (z_filt + Phi_filt * theta - u_filt);
    }
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << "\t" << theta.transpose() << endl;
}

Eigen::VectorXd RCAC::compute_uk(double z, double z_int, double z_diff, double u)
{
    set_RCAC_data(z, u);
    buildRegressor(z, z_int, z_diff);
    filter_data();
    update_theta();
    u_k = Phi_k * theta;
    // shift_data();
    kk = kk + 1;
    return u_k;
}