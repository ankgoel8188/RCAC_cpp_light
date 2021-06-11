#pragma once
#include "matrix/math.hpp"

// using namespace matrix;
// using namespace std;

/**
 * The parent RCAC class. This class handles all the low level computation of RCAC
 * such as the filtering, coefficient updates, and keeping track of the regressors.
 *
 * The notation here follows the JGCD 2019 implementation
 */
template<size_t MatDim>
class RCAC
{
    float P0;
    float lambda;
    float N_nf;

public:
    RCAC();
    RCAC(float P0_val);
    RCAC(float P0_val, float lambda_val, float N_nf_val);

    ~RCAC() = default;
    RCAC(const RCAC & obj);
    RCAC& operator=(const RCAC & obj);

    int   getkk() {return kk;};
    float get_rcac_uk() {return u_k;};
    float get_rcac_zk() {return z_k;};

    //TODO: FIX THESE GET FUNCTIONS, NOT FOOLPROOF
    float get_rcac_theta(int i) {return theta(i,0);}
    float get_rcac_P(int i, int j){return P(i, j);};
    float get_rcac_Phi(int i) {return Phi_k(0,i);}

    void set_RCAC_data(float, float);
    void buildRegressor(float zkm1, float zkm1_int, float zkm1_diff);
    void filter_data();
    void update_theta();

    float compute_uk(float _z_in, matrix::Matrix<float, 1, MatDim> _phi_in, float _u_km1_in);

protected:
    const int nf = 2;
    matrix::Matrix<float, 1, 2> filtNu;         // 1st order filter. Gf = filtNu(0) + filtNu(1)/q
                                                // In most cases, filtNu(0) = 0 and filtNu(1) = +-1

    //RCAC internal variables
    matrix::Matrix<float, MatDim, MatDim> P;
    matrix::Matrix<float, MatDim, 1> theta;          // Kp = theta(0,0)
                                                // Ki = theta(1,0)
                                                // Kd = theta(2,0)

    float u_k, u_km1, u_filt;
    float z_k, z_filt;
    matrix::Matrix<float, 1, MatDim> Phi_k, Phi_filt;

    matrix::Matrix<float, 2, 1> ubar;        // Size nf by 1
    matrix::Matrix<float, 3, MatDim> Phibar;      // Size nf+1 by 1
    matrix::Matrix<float, 2, MatDim> PhibarBlock; // Size nf by 1

    float Gamma;
    float Idty_lz;
    matrix::Matrix<float, 1, 1> one_matrix;
    matrix::Matrix<float, 1, 1> dummy;

    int kk = 0;
};


// TEST: Template is not working correctly, so temp fix

template<size_t MatDim>
RCAC<MatDim>::RCAC() : RCAC(0.1, 1.0, 1.0) {}

template<size_t MatDim>
RCAC<MatDim>::RCAC(float P0_val) : RCAC(P0_val, 1.0, 1.0) {}

template<size_t MatDim>
RCAC<MatDim>::RCAC(float P0_val, float lambda_val, float N_nf_val) :
    P0(P0_val), lambda(lambda_val), N_nf(N_nf_val)
{
    // Initialize interal RCAC variables
    P = matrix::eye<float, MatDim>() * P0;
    theta.setZero();
    filtNu.setZero();
    filtNu(0,nf-1)=N_nf;
    u_k = 0;
    u_km1 = 0;
    u_filt = 0;
    z_k = 0;
    z_filt = 0;

    Gamma = 0;
    Idty_lz = 0;

    Phi_k.setZero();
    Phi_filt.setZero();

    ubar.setZero();
    Phibar.setZero();
    PhibarBlock.setZero();

    one_matrix = matrix::eye<float, 1>();
    dummy.setZero();
}

template<size_t MatDim>
RCAC<MatDim>::RCAC(const RCAC & obj)
{
    P0 = obj.P0;
    filtNu = obj.filtNu;
    P = obj.P;
    theta = obj.theta;
    u_k = obj.u_k;
    u_km1 = obj.u_km1;
    u_filt = obj.u_filt;
    z_k = obj.z_k;
    z_filt = obj.z_filt;
    Phi_k = obj.Phi_k;
    Phi_filt = obj.Phi_filt;
    ubar = obj.ubar;
    Phibar = obj.Phibar;
    PhibarBlock = obj.PhibarBlock;
    Gamma = obj.Gamma;
    Idty_lz = obj.Idty_lz;
    one_matrix = obj.one_matrix;
    dummy = obj.dummy;
}

template<size_t MatDim>
RCAC<MatDim>& RCAC<MatDim>::operator=(const RCAC & obj)
{
    P0 = obj.P0;
    filtNu = obj.filtNu;
    P = obj.P;
    theta = obj.theta;
    u_k = obj.u_k;
    u_km1 = obj.u_km1;
    u_filt = obj.u_filt;
    z_k = obj.z_k;
    z_filt = obj.z_filt;
    Phi_k = obj.Phi_k;
    Phi_filt = obj.Phi_filt;
    ubar = obj.ubar;
    Phibar = obj.Phibar;
    PhibarBlock = obj.PhibarBlock;
    Gamma = obj.Gamma;
    Idty_lz = obj.Idty_lz;
    one_matrix = obj.one_matrix;
    dummy = obj.dummy;
    return *this;
}

template<size_t MatDim>
void RCAC<MatDim>::set_RCAC_data(float z_k_val, float u_km1_val)
{
    z_k = z_k_val;
    u_km1 = u_km1_val;
}

template<size_t MatDim>
void RCAC<MatDim>::buildRegressor(float z, float z_int, float z_diff)
{
    Phi_k(0, 0) = z;
    Phi_k(0, 1) = z_int;
    Phi_k(0, 2) = z_diff;
    // for (size_t i = 0; i < MatDim; ++i)
    // {
    //     Phi_k(0, i) =
    // }
}

template<size_t MatDim>
void RCAC<MatDim>::filter_data()
{
    //TODO: CHECK LOGIC HERE.
    for (int ii = nf - 1; ii > 0; ii--)
    {
        ubar(ii, 0) = ubar(ii - 1, 0);
    }
    ubar(0, 0) = u_km1;

    for (int ii = nf; ii > 0; ii--)
    {
        for (int jj = 0; jj < (int)MatDim; jj++)
        {
            Phibar(ii, jj) = Phibar(ii - 1, jj);
        }
    }
    for (int jj = 0; jj < (int)MatDim; jj++)
    {
        Phibar(0, jj) = Phi_k(0, jj);
    }

    z_filt = z_k;

    dummy = filtNu * ubar;
    u_filt = dummy(0, 0);

    for (int ii = 0; ii < nf; ii++)
    {
        for (int jj = 0; jj < (int)MatDim; jj++)
        {
            PhibarBlock(ii, jj) = Phibar(ii + 1, jj);
        }
    }
    Phi_filt = filtNu * PhibarBlock;
}


template<size_t MatDim>
void RCAC<MatDim>::update_theta()
{
    if (kk > 3)
    {
        // Phi_filt and P incompatible
        dummy = Phi_filt * P * Phi_filt.transpose();
        Gamma = lambda + dummy(0, 0);
        P = P - P * Phi_filt.transpose() * 1 / Gamma * Phi_filt * P;
        P = P / lambda;

        theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix + Phi_filt * theta - u_filt);
    }
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << "\t" << theta.transpose() << endl;
}

template<size_t MatDim>
float RCAC<MatDim>::compute_uk(float _z_in, matrix::Matrix<float, 1, MatDim> _phi_in, float _u_km1_in)
{
    // std::cout << one_matrix(0, 0) << std::endl;
    // set_RCAC_data(z, u);
    // buildRegressor(z, z_int, z_diff);
    Phi_k = _phi_in;
    z_k = _z_in;
    u_km1 = _u_km1_in;

    filter_data();
    update_theta();
    dummy = Phi_k * theta;
    u_k = dummy(0, 0);
    kk = kk + 1;
    return u_k;
}
