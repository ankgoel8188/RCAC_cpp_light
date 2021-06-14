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
template<size_t l_theta, size_t l_Rblock>
class RCAC
{
    //TO DO: Implement R as variable vector to allow user to specify either one or two values for diagonal elements in matrix containing Rz, Ru
    float P0;
    float lambda;
    float Rz;
    float Ru;
    float N_nf;
    int e_fun_num;

public:
    RCAC();
    RCAC(float P0_val);
    RCAC(float P0_val, float lambda_val, float Rz_val, float Ru_val, float N_nf_val);

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

    // void set_RCAC_data(float, float);
    //void buildRegressor(float zkm1, float zkm1_int, float zkm1_diff);
    void normalize_e();
    void filter_data();
    void build_Rblock();
    void update_theta();

    float compute_uk(float _z_in, matrix::Matrix<float, 1, l_theta> _phi_in, float _u_km1_in, int e_fun_num_in);

protected:
    const int nf = 2;
    float mu = 1.0;
    float nu = 1.0;

    matrix::Matrix<float, 1, 2> filtNu;         // 1st order filter. Gf = filtNu(0) + filtNu(1)/q
                                                // In most cases, filtNu(0) = 0 and filtNu(1) = +-1

    //RCAC internal variables
    matrix::Matrix<float, l_theta, l_theta> P;
    matrix::Matrix<float, l_theta, 1> theta;          // Kp = theta(0,0)
                                                // Ki = theta(1,0)
                                                // Kd = theta(2,0)

    float u_k, u_km1, u_filt;
    float z_k, z_filt;
    matrix::Matrix<float, 1, l_theta> Phi_k, Phi_filt;

    matrix::Matrix<float, 2, 1> ubar;        // Size nf by 1
    matrix::Matrix<float, 3, l_theta> Phibar;      // Size nf+1 by 1
    matrix::Matrix<float, 2, l_theta> PhibarBlock; // Size nf by 1

    float Gamma;
    float Idty_lz;
    matrix::Matrix<float, 1, 1> one_matrix;
    matrix::Matrix<float, 1, 1> dummy;
    matrix::Matrix<float, l_Rblock, l_theta> Phiblock;
    matrix::Matrix<float, l_Rblock, l_Rblock> Rblock, PhiB_P_PhiB_t;

    int kk = 0;
};


// TEST: Template is not working correctly, so temp fix

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC() : RCAC(0.1, 1.0, 1.0, 1.0, 1.0) {}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(float P0_val) : RCAC(P0_val, 1.0, 1.0, 1.0, 1.0) {}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(float P0_val, float lambda_val, float Rz_val, float Ru_val, float N_nf_val) :
    P0(P0_val), lambda(lambda_val), Rz(Rz_val), Ru(Ru_val), N_nf(N_nf_val)
{
    // Initialize interal RCAC variables
    P = matrix::eye<float, l_theta>() * P0;
    theta.setZero();
    filtNu.setZero();
    filtNu(0,nf-1) = N_nf;
    u_k = 0;
    u_km1 = 0;
    u_filt = 0;
    z_k = 0;
    z_filt = 0;

    Gamma = 0;
    Idty_lz = 0;

    // set Ru, Rz to 1 for testing
    // Ru = 1;
    // Rz = 1;

    Phi_k.setZero();
    Phi_filt.setZero();

    ubar.setZero();
    Phibar.setZero();
    PhibarBlock.setZero();

    one_matrix = matrix::eye<float, 1>();
    dummy.setZero();
    Phiblock.setZero();
    PhiB_P_PhiB_t.setZero();
    Rblock.setZero();
}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(const RCAC & obj)
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
    Phiblock = obj.Phiblock;
    PhiB_P_PhiB_t = obj.PhiB_P_PhiB_t;
    Rblock = obj.Rblock;
    Ru = obj.Ru;
    Rz = obj.Rz;
}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>& RCAC<l_theta, l_Rblock>::operator=(const RCAC & obj)
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
    Phiblock = obj.Phiblock;
    PhiB_P_PhiB_t = obj.PhiB_P_PhiB_t;
    Rblock = obj.Rblock;
    Ru = obj.Ru;
    Rz = obj.Rz;
    return *this;
}

// template<size_t l_theta, size_t l_Rblock>
// void RCAC<l_theta, l_Rblock>::set_RCAC_data(float z_k_val, float u_km1_val)
// {
//     z_k = z_k_val;
//     u_km1 = u_km1_val;
// }

// template<size_t l_theta, size_t l_Rblock>
// void RCAC<l_theta, l_Rblock>::buildRegressor(float z, float z_int, float z_diff)
// {
//     Phi_k(0, 0) = z;
//     Phi_k(0, 1) = z_int;
//     Phi_k(0, 2) = z_diff;
//     // for (size_t i = 0; i < l_theta; ++i)
//     // {
//     //     Phi_k(0, i) =
//     // }
// }

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::normalize_e()
{
    float pi = (float)M_PI;

    switch (e_fun_num)
    {
        case 1:
            z_k = (mu * nu * z_k) / (mu + nu * abs(z_k));    //norm_e_fun_1(Phi_k(0,i));
            break;
        case 2:
            z_k = ((2*nu)/pi)*(float)atan((pi*nu*z_k)/(2*mu));    //norm_e_fun_2(Phi_k(0,i));
            break;
        case 3:
            z_k = (mu*nu*z_k)/(float)sqrt(mu*mu + nu*nu*z_k*z_k);    //norm_e_fun_3(Phi_k(0,i));
            break;
        case 4:
            z_k = mu*(float)tanh((nu*z_k)/mu);    //norm_e_fun_4(Phi_k(0,i));
            break;
        case 5:
            // TO DO: double check
            z_k = mu * (float)erf(((float)sqrt(pi)*nu*z_k)/(2*mu));  //*integrate_rk4(
            break;
        default:
            break;
    }
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::filter_data()
{
    //TODO: CHECK LOGIC HERE.
    for (int ii = nf - 1; ii > 0; ii--)
    {
        ubar(ii, 0) = ubar(ii - 1, 0);
    }
    ubar(0, 0) = u_km1;

    for (int ii = nf; ii > 0; ii--)
    {
        for (int jj = 0; jj < (int)l_theta; jj++)
        {
            Phibar(ii, jj) = Phibar(ii - 1, jj);
        }
    }
    for (int jj = 0; jj < (int)l_theta; jj++)
    {
        Phibar(0, jj) = Phi_k(0, jj);
    }

    z_filt = z_k;

    dummy = filtNu * ubar;
    u_filt = dummy(0, 0);

    for (int ii = 0; ii < nf; ii++)
    {
        for (int jj = 0; jj < (int)l_theta; jj++)
        {
            PhibarBlock(ii, jj) = Phibar(ii + 1, jj);
        }
    }
    Phi_filt = filtNu * PhibarBlock;
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::build_Rblock()
{
    for (int i = 0; i < (int)l_theta ; i++)
    {
        Phiblock(0,i) = Phi_filt(0, i);
        Phiblock(1,i) = Phi_k(0,i);
    }

    Rblock(0,0) = Rz;
    Rblock(1,1) = Ru;
}


template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::update_theta()
{
    if (kk > 3)
    {
        if ((int)l_Rblock == 2)
        {
            PhiB_P_PhiB_t = Phiblock * P * Phiblock.transpose();
            PhiB_P_PhiB_t = geninv(Rblock) + PhiB_P_PhiB_t;
            P = P - P * Phiblock.transpose() * geninv(PhiB_P_PhiB_t) * Phiblock * P;

            theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix - u_filt + Phi_filt * theta) - P * Phi_k.transpose() * Ru * Phi_k * theta;
        }
        else
        {
            dummy = Phi_filt * P * Phi_filt.transpose();
            Gamma = lambda + dummy(0, 0);
            P = P - P * Phi_filt.transpose() * 1 / Gamma * Phi_filt * P;
            P = P / lambda;

            theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix + Phi_filt * theta - u_filt);
        }
    }
    // cout << kk << "\t" << z_filt << "\t" << u_filt << "\t" << Phi_filt << "\t" << theta.transpose() << endl;
}




template<size_t l_theta, size_t l_Rblock>
float RCAC<l_theta, l_Rblock>::compute_uk(float _z_in, matrix::Matrix<float, 1, l_theta> _phi_in, float _u_km1_in, int e_fun_num_in)
{
    // std::cout << one_matrix(0, 0) << std::endl;

    // buildRegressor(z, z_int, z_diff);
    Phi_k = _phi_in;
    z_k = _z_in;
    u_km1 = _u_km1_in;
    e_fun_num = e_fun_num_in;

    normalize_e();
    filter_data();
    if ((int)l_Rblock == 2)
    {
        build_Rblock();
    }
    update_theta();

    dummy = Phi_k * theta;
    u_k = dummy(0, 0);
    kk = kk + 1;

    return u_k;
}
