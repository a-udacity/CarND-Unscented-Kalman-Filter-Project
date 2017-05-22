#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class UKFLaser : public UKF {

    void UKFLaser::ProcessMeasurement(MeasurementPackage meas_package) {
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
};
