// #pragma once

// #include "matrix/math.hpp"
// #include <cfloat>
// #include <cmath>
// #include <algorithm>
#include "RCACParams.h"

// #include <iostream>

// using namespace matrix;
// using namespace std;

// NOTE: This Header is Unneccesary for uses outside of PX4.

#include <px4_platform_common/defines.h>

/**
 * The parent RCAC class. This class handles all the low level computation of RCAC
 * such as the filtering, coefficient updates, and keeping track of the regressors.
 *
 * The notation here follows the JGCD 2019 implementation
 */

template<size_t l_theta, size_t l_Rblock>
class RCAC
{
public:
    RCAC();
    RCAC(float P0_val);
    RCAC(const RCACParams & RCAC_Parameters_in);
    // RCAC(float P0_val, float lambda_val, float N_nf_val, int e_fun_num_val, float lim_int_val = FLT_MAX);
    // RCAC(float P0_val, float lambda_val, matrix::Matrix<float, l_Rblock, l_Rblock> Rblock_val, float N_nf_val, int e_fun_num_val, float lim_int_val = FLT_MAX);

    ~RCAC() = default;
    RCAC(const RCAC & obj);
    RCAC& operator=(const RCAC & obj);

    int   getkk() {return kk;};
    float get_rcac_uk() {return u_k;};
    float get_rcac_zk() {return z_k;};

    //TODO: FIX THESE GET FUNCTIONS, NOT FOOLPROOF
    float get_rcac_theta(int i) {return theta(i,0);}
    float get_rcac_P(int i, int j){return P(i, j);};
    // float get_rcac_Ru(){return Rblock(1,1);};
    float get_rcac_Phi(int i) {return Phi_k(0,i);}
    float get_rcac_integral() {return rcac_int;};
    // float get_rcac_N() {return _RCACParams.tuneParams.N_nf;}
    const RCACParams & get_rcacParams( return _RCACParams; )

    void set_lim_int(float lim_int_in);
    void normalize_e();

    void filter_data();
    void build_Phiblock();
    void update_theta_Rblock_ON();
    void update_theta_Rblock_OFF();
    void update_integral(const float rcac_error, const float dt);
    void reset_integral();
    void reset_kk();
    void init_var_helper();

    float compute_uk(float _z_in, matrix::Matrix<float, 1, l_theta> _phi_in, float _u_km1_in);

protected:
    matrix::Matrix<float, l_Rblock, l_Rblock> Rblock;
    const int nf = 2;
    float mu = 1.0;
    float nu = 1.0;
    int kk;

    RCACParams _RCACParams;

    matrix::Matrix<float, 1, 2> filtNu;         // 1st order filter. Gf = filtNu(0) + filtNu(1)/q
                                                // In most cases, filtNu(0) = 0 and filtNu(1) = +-1

    //RCAC internal variables
    matrix::Matrix<float, l_theta, l_theta> P;
    matrix::Matrix<float, l_theta, 1> theta;    // Kp = theta(0,0)
                                                // Ki = theta(1,0)
                                                // Kd = theta(2,0)

    float u_k, u_km1, u_filt;
    float z_k, z_filt;
    matrix::Matrix<float, 1, l_theta> Phi_k, Phi_filt;

    matrix::Matrix<float, 2, 1> ubar;               // Size nf by 1
    matrix::Matrix<float, 3, l_theta> Phibar;       // Size nf+1 by 1
    matrix::Matrix<float, 2, l_theta> PhibarBlock;  // Size nf by 1

    float Gamma;
    float Idty_lz;
    matrix::Matrix<float, 1, 1> one_matrix;
    matrix::Matrix<float, 1, 1> dummy;
    matrix::Matrix<float, l_Rblock, l_theta> Phiblock;
    matrix::Matrix<float, l_Rblock, l_Rblock> PhiB_P_PhiB_t;

    float rcac_int;

    //TODO: Initialize lim_int
};

// TEST: Template is not working correctly, so temp fix
template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC() : RCAC(0.1) {}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(float P0_val)
{
    init_var_helper();
    Rblock(0, 0) = _RCACParams.tuneParams.Ru;
    Rblock(1, 1) = _RCACParams.initParams.Rz;
    P = matrix::eye<float, l_theta>() * P0_val;
    filtNu(0,nf-1) = _RCACParams.tuneParams.N_nf;
}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(const RCACParams & RCAC_Parameters_in) : _RCACParams(RCAC_Parameters_in)
{
    init_var_helper();
    Rblock(0, 0) = RCAC_Parameters_in.tuneParams.Ru;
    Rblock(1, 1) = RCAC_Parameters_in.initParams.Rz;
    P = matrix::eye<float, l_theta>() * RCAC_Parameters_in.tuneParams.p0;
    filtNu(0,nf-1) = RCAC_Parameters_in.tuneParams.N_nf;
}


template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>::RCAC(const RCAC & obj)
{
    _RCACParams = obj._RCACParams;
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
    kk = obj.kk;
}

template<size_t l_theta, size_t l_Rblock>
RCAC<l_theta, l_Rblock>& RCAC<l_theta, l_Rblock>::operator=(const RCAC & obj)
{
    _RCACParams = obj._RCACParams;
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
    kk = obj.kk;
    return *this;
}


template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::init_var_helper()
{
    theta.setZero();
    filtNu.setZero();

    u_k = 0;
    u_km1 = 0;
    u_filt = 0;
    z_k = 0;
    z_filt = 0;
    kk = 0;

    Gamma = 0;
    Idty_lz = 0;

    Phi_k.setZero();
    Phi_filt.setZero();

    ubar.setZero();
    Phibar.setZero();
    PhibarBlock.setZero();

    one_matrix = matrix::eye<float, 1>();
    dummy.setZero();
    Phiblock.setZero();
    PhiB_P_PhiB_t.setZero();
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::normalize_e()
{
    float pi = (float)M_PI;

    // std::cout << "\nError normalization Update:\t" << e_fun_num;
    switch (_RCACParams.initParams.errorNormMode)
    {
        case 1:
            z_k = (mu * nu * z_k) / (mu + nu * abs(z_k));
            break;
        case 2:
            z_k = ((2*nu)/pi)*(float)atan((pi*nu*z_k)/(2*mu));
            break;
        case 3:
            z_k = (mu*nu*z_k)/(float)sqrt(mu*mu + nu*nu*z_k*z_k);
            break;
        case 4:
            z_k = mu*(float)tanh((nu*z_k)/mu);
            break;
        case 5:
            z_k = mu * (float)erf(((float)sqrt(pi)*nu*z_k)/(2*mu));
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
void RCAC<l_theta, l_Rblock>::build_Phiblock()
{
    for (int i = 0; i < (int)l_theta ; i++)
    {
        // std::cout << "\nPhiblock update - Rblock ON\n";
        Phiblock(0,i) = Phi_filt(0, i);
        Phiblock(1,i) = Phi_k(0,i);
    }
}


template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::update_theta_Rblock_ON()
{
    if (kk > 3)
    {
        // std::cout << "\nTheta update - Rblock ON\n";
        PhiB_P_PhiB_t = Phiblock * P * Phiblock.transpose();
        PhiB_P_PhiB_t = geninv(Rblock) + PhiB_P_PhiB_t;
        P = P - P * Phiblock.transpose() * geninv(PhiB_P_PhiB_t) * Phiblock * P;

        theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix - u_filt + Phi_filt * theta) - P * Phi_k.transpose() * Rblock(1,1) * Phi_k * theta;
    }
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::update_theta_Rblock_OFF()
{
    if (kk > 3)
    {
        // std::cout << "\nTheta update - Rblock OFF\n";
        dummy = Phi_filt * P * Phi_filt.transpose();
        Gamma = _RCACParams.initParams.lambda + dummy(0, 0);
        P = P - P * Phi_filt.transpose() * 1 / Gamma * Phi_filt * P;
        P = P / _RCACParams.initParams.lambda;

        theta = theta - P * Phi_filt.transpose() * (z_filt * one_matrix + Phi_filt * theta - u_filt);
    }
}

template<size_t l_theta, size_t l_Rblock>
float RCAC<l_theta, l_Rblock>::compute_uk(float _z_in, matrix::Matrix<float, 1, l_theta> _phi_in, float _u_km1_in)
{
    // Check RCAC_EN to see if the RCAC is enabled.
    if (!_RCACParams.tuneParams.RCAC_EN)
        return 0;

    u_km1 = _u_km1_in;
    Phi_k = _phi_in;
    z_k = _z_in;

    if (_RCACParams.initParams.errorNormMode > 0)
    {
        normalize_e();
    }

    filter_data();

    if (_RCACParams.initParams.RBlock_EN)
    {
        build_Phiblock();
        update_theta_Rblock_ON();
    }
    else
    {
        update_theta_Rblock_OFF();
    }

    dummy = Phi_k * theta;
    u_k = dummy(0, 0);
    kk++;

    return u_k;
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::update_integral(const float rcac_error, const float dt)
{
    // Update Jun 24, 2020 - Removed Anti-windup out of the class. Anti-windup needs to be implemented out of class.
    float rcac_i = rcac_int + rcac_error * dt;

    // do not propagate the result if out of range or invalid
    if (PX4_ISFINITE(rcac_i)) {
        rcac_int = math::constrain(rcac_i, -_RCACParams.initParams.lim_int, _RCACParams.initParams.lim_int);
    }

}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::reset_integral()
{
    rcac_int = 0;
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::set_lim_int(float lim_int_in)
{
    _RCACParams.initParams.lim_int = lim_int_in;
}

template<size_t l_theta, size_t l_Rblock>
void RCAC<l_theta, l_Rblock>::reset_kk()
{
    kk = 0;
}


template<size_t l_theta_IO, size_t l_RBlock_IO>
class RCAC_Public_IO {

	private:
	RCAC<l_theta_IO, l_RBlock_IO> * RCAC_ptr;

	public:
    RCAC_Public_IO() : RCAC_PUblic_IO(nullptr) {}

    RCAC_Public_IO(RCAC * RCAC_ptr_in) : RCAC_ptr(RCAC_ptr_in) {}

    float get_P11() {
        return RCAC_ptr->get_rcac_P(0, 0);
    }

    float get_uk() {
        return RCAC_ptr->get_rcac_uk();
    }

    float get_theta(size_t i) {
        return RCAC_ptr->get_rcac_theta(i);
    }

    size_t get_kk() {
        return RCAC_ptr->getkk();
    }


};
