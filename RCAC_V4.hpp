#ifndef _RCAC_new_HPP_
#define _RCAC_new_HPP_

#include "Eigen/Dense"
#include "Eigen/Sparse"
// #include <iostream>


/**
 * The parent RCAC class. This class handles all the low level computation of RCAC
 * such as the filtering, coefficient updates, and keeping track of the regressors.
 * 
 * Almost all the methods are polymorphic and can be modified by child classes
 * to create RCAC algorithms with more complex filtering.
 */
class RCAC
{
    double P0;
    double lambda;
    double N1;

public:
    RCAC( double, double, double);
    
    /**
         * Returns RCAC's computed value for the control. Must run oneStep at least once.
         */
    //Function getControl: Get the computed control input
    Eigen::VectorXd get_uk()
    {
        return u_k;
    };

    /**
         * Returns a vector of the current RCAC coefficients.
         */
    //Function getCoeff: Get the RCAC coefficients
    Eigen::VectorXd get_thetak()
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

    Eigen::MatrixXd getPhi()
    {
        return Phi_k;
    };

    // void buildRegressor(Eigen::VectorXd &zIn);

    void init_RCAC();
    void buildRegressor(double zk, double zk_int, double zk_diff);
    void set_performance(double zk);
    void update_theta();
    Eigen::VectorXd compute_uk(double zk, double zk_int, double zk_diff);
    void shift_data();
        
protected:
    //RCAC Working variables
    Eigen::MatrixXd P;
    Eigen::VectorXd theta;

    Eigen::VectorXd u_k, u_km1, u_km2, u_filt;
    Eigen::VectorXd z_k, z_km1, z_filt;
    Eigen::MatrixXd Phi_k, Phi_km1, Phi_km2, Phi_filt; 
    Eigen::MatrixXd Gamma;
    Eigen::MatrixXd Idty_lz;
    int kk = 1;
};

#endif