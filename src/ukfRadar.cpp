#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

class UKFRadar : public UKF {
    /**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
    void UKFRadar::ProcessMeasurement(MeasurementPackage meas_package) {
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

};
