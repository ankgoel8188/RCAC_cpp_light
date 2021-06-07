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
    float get_rcac_theta(int i) {return theta(i,0);}
    float get_rcac_P(int i, int j){return P(i, j);};
    float get_rcac_Phi(int i) {return Phi_k(i,0);}

    void set_RCAC_data(float, float);
    void buildRegressor(float zkm1, float zkm1_int, float zkm1_diff);
    void filter_data();
    void update_theta();
    float compute_uk(float z, float z_int, float z_diff, float u);

protected:
    const int nf = 2;
    matrix::Matrix<float, 1, 2> filtNu;         // 1st order filter. Gf = filtNu(0) + filtNu(1)/q
                                                // In most cases, filtNu(0) = 0 and filtNu(1) = +-1

    //RCAC internal variables
    matrix::Matrix<float, 3, 3> P;
    matrix::Matrix<float, 3, 1> theta;          // Kp = theta(0,0)
                                                // Ki = theta(1,0)
                                                // Kd = theta(2,0)

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

    int kk = 0;
};
