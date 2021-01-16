#include "RCAC.h"
// #include <iostream>

// using namespace std;

// RCAC::RCAC(float P0_val, float lambda_val, int nf_val, float N1_val)
// {
//     P0 = P0_val;
//     lambda = lambda_val;
//     nf = nf_val;
// //     filtNu = filtNu_val;
//     N1 = N1_val;
// }

void RCAC::set_RCAC_parameters(float P0_val, float lambda_val, int nf_val, float N1_val)
{
    P0 = P0_val;
    lambda = lambda_val;
    nf = nf_val;
//     filtNu = filtNu_val;
N1 = N1_val;
}

void RCAC::init_RCAC()
{
    P = eye<float, 3>() * P0;
    theta.setZero();
    filtNu(0,0)=0;
    filtNu(0,1)=N1;
    u_km1 = 0;
    z_km1 = 0;

    Phi_k.setZero();

    ubar.setZero();
    Phibar.setZero();

    Phi_filt.setZero();

    one_matrix = eye<float, 1>();
}

void RCAC::set_RCAC_data(float zkm1, float ukm1)
{
    z_km1 = zkm1;
    u_km1 = ukm1;
}

void RCAC::buildRegressor(float z, float z_int, float z_diff)
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

float RCAC::compute_uk(float z, float z_int, float z_diff, float u)
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
