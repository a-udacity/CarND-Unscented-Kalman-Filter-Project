#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "tools.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;

    ///* state covariance matrix
    MatrixXd P_;

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    ///* time when the state is true, in us
    long long time_us_;

    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;

    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;

    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;

    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;

    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_;

    ///* Weights of sigma points
    VectorXd weights_;

    ///* State dimension
    int n_x_;

    ///* Augmented state dimension
    int n_aug_;

    ///* Sigma point spreading parameter
    double lambda_;

    ///* the current NIS for radar
    double NIS_radar_;

    ///* the current NIS for laser
    double NIS_laser_;

    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);

    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);

    void predictedStateCovariance();

    void predictedStateMean();
};

class UKFLaser : public UKF {

public :
    void ProcessMeasurement(MeasurementPackage meas_package) {
        /**
        Complete this function! Make sure you switch between lidar and radar
        measurements.
        */
        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        if (!is_initialized_) {
            // first measurement
            /**
             Initialize state.
             */
            float px = meas_package.raw_measurements_(0);
            float py = meas_package.raw_measurements_(1);
            x_(0) = px;
            x_(1) = py;
            P_(0, 0) = std_laspx_;  // used lidar distance uncertainty px as initial px covariance uncertainty
            P_(1, 1) = std_laspy_;  // used lidar distance uncertainty py as initial py covariance uncertainty

            time_us_ = meas_package.timestamp_;

            // done initializing, no need to predict or update
            is_initialized_ = true;
            return;
        }
        /*****************************************************************************
         *  Prediction
         ****************************************************************************/
        //Time is measured in seconds.
        //compute the time elapsed between the current and previous measurements
        float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    // expressed in seconds
        time_us_ = meas_package.timestamp_;
        Prediction(dt);
        /*****************************************************************************
         *  Predict Sensor Measurements and Update
         ****************************************************************************/
        UpdateLidar(meas_package);
    }
    /**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
public:
    void UpdateLidar(MeasurementPackage meas_package) {
        /**
        Complete this function! Use lidar data to update the belief about the object's
        position. Modify the state vector, x_, and covariance, P_.

        You'll also need to calculate the lidar NIS.
        */
        /*****************************************************************************
         *  Predict Lidar Measurements
         ****************************************************************************/
        //set measurement dimension, lidar can measure px, and py
        int n_z = 2;

        //create matrix for sigma points in measurement space
        MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

        //transform sigma points into measurement space
        for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  //2n+1 simga points

            // extract values for better readibility
            double p_x = Xsig_pred_(0, i);
            double p_y = Xsig_pred_(1, i);

            // measurement model
            Zsig(0, i) = p_x;  //px
            Zsig(1, i) = p_y;  //py
        }

        //mean predicted measurement
        VectorXd z_pred = VectorXd(n_z);
        z_pred.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            z_pred += weights_(i) * Zsig.col(i);
        }

        //measurement covariance matrix S
        MatrixXd S = MatrixXd(n_z, n_z);
        S.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;

            S += weights_(i) * z_diff * z_diff.transpose();
        }

        //add measurement noise covariance matrix
        MatrixXd R = MatrixXd(n_z, n_z);
        R << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;
        S += R;

        /*****************************************************************************
         *  Update Lidar Measurements
         ****************************************************************************/
        //create matrix for cross correlation Tc
        MatrixXd Tc = MatrixXd(n_x_, n_z);
        //calculate cross correlation matrix
        Tc.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;

            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            //angle normalization
            while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
            while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

            Tc += weights_(i) * x_diff * z_diff.transpose();
        }

        //Kalman gain K;
        MatrixXd K = Tc * S.inverse();

        //residual
        VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

        //update state mean and covariance matrix
        x_ += K * z_diff;
        P_ += -K * S * K.transpose();

        /*****************************************************************************
         *  Compute Normalized Innovation Squared (NIS)
         ****************************************************************************/
        NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
    }
};

class UKFRadar : public UKF {
    /**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
public :
    void ProcessMeasurement(MeasurementPackage meas_package) {
        /**
        Complete this function! Make sure you switch between lidar and radar
        measurements.
        */
        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        if (!is_initialized_) {
            // first measurement
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            float ro = meas_package.raw_measurements_(0);
            float phi = meas_package.raw_measurements_(1);
            //float ro_dot = meas_package.raw_measurements_(2);
            x_(0) = ro * cos(phi);
            x_(1) = ro * sin(phi);
            P_(0, 0) = std_radr_;  // used radar distance uncertainty r as initial px covariance uncertainty
            P_(1, 1) = std_radr_;  // used radar distance uncertainty r as initial px covariance uncertainty

            time_us_ = meas_package.timestamp_;

            // done initializing, no need to predict or update
            is_initialized_ = true;
            return;
        }

        /*****************************************************************************
         *  Prediction
         ****************************************************************************/
        //Time is measured in seconds.
        //compute the time elapsed between the current and previous measurements
        float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;    // expressed in seconds
        time_us_ = meas_package.timestamp_;
        Prediction(dt);

        /*****************************************************************************
         *  Predict Sensor Measurements and Update
         ****************************************************************************/
        UpdateRadar(meas_package);
    }
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
public:
    void UpdateRadar(MeasurementPackage meas_package) {
        /**
        Complete this function! Use radar data to update the belief about the object's
        position. Modify the state vector, x_, and covariance, P_.

        You'll also need to calculate the radar NIS.
        */
        /*****************************************************************************
         *  Predict Radar Measurements
         ****************************************************************************/
        //set measurement dimension, radar can measure r, phi, and r_dot
        int n_z = 3;

        //create matrix for sigma points in measurement space
        MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

        //transform sigma points into measurement space
        for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  //2n+1 simga points

            // extract values for better readibility
            double p_x = Xsig_pred_(0, i);
            double p_y = Xsig_pred_(1, i);
            double v = Xsig_pred_(2, i);
            double yaw = Xsig_pred_(3, i);

            double v1 = cos(yaw) * v;
            double v2 = sin(yaw) * v;

            // measurement model
            Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                        //r
            Zsig(1, i) = atan2(p_y, p_x);                                 //phi
            Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
        }

        //mean predicted measurement
        VectorXd z_pred = VectorXd(n_z);
        z_pred.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            z_pred += weights_(i) * Zsig.col(i);
        }

        //measurement covariance matrix S
        MatrixXd S = MatrixXd(n_z, n_z);
        S.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;

            //angle normalization
            while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
            while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

            S += weights_(i) * z_diff * z_diff.transpose();
        }

        //add measurement noise covariance matrix
        MatrixXd R = MatrixXd(n_z, n_z);
        R << std_radr_ * std_radr_, 0, 0,
                0, std_radphi_ * std_radphi_, 0,
                0, 0, std_radrd_ * std_radrd_;
        S += R;

        /*****************************************************************************
         *  Update Radar Measurements
         ****************************************************************************/
        //create matrix for cross correlation Tc
        MatrixXd Tc = MatrixXd(n_x_, n_z);
        //calculate cross correlation matrix
        Tc.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
            //residual
            VectorXd z_diff = Zsig.col(i) - z_pred;
            //angle normalization
            while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
            while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            //angle normalization
            while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
            while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

            Tc += weights_(i) * x_diff * z_diff.transpose();
        }

        //Kalman gain K;
        MatrixXd K = Tc * S.inverse();

        //residual
        VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        //update state mean and covariance matrix
        x_ += K * z_diff;
        P_ += -K * S * K.transpose();

        NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
    }


};

#endif /* UKF_H */
