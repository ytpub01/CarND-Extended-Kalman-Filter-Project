#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  /* TODO DONE: Finish initializing the FusionEKF.
   * TODO DONE: Set the process noise omega ~ N(0,Q)
           and measurement noise nu ~ N(0,R) */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
  ekf_.F_ = MatrixXd(4, 4);
  dt = 0.05;
  noise_ax = 9.0;
  noise_ay = 9.0;
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  // for process noise Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt*dt/4*noise_ax, 0, dt/2*noise_ax, 0,
            0, dt*dt/4*noise_ay, 0, dt/2*noise_ay,
            dt/2*noise_ax, 0, noise_ax, 0,
            0, dt/2*noise_ay, 0, noise_ay;
  ekf_.Q_ *= dt*dt;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
   /* Initialization */
  if (!is_initialized_) {
    /* TODO DONE: Initialize the state ekf_.x_ with the first measurement
     * TODO DONE above: Create the covariance matrix (P).
     * You'll need to convert radar from polar to cartesian coordinates.*/
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO DONE: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << rho*cos(phi), rho*sin(phi), rho_dot*cos(phi), rho_dot*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO DONE: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO DONE: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO DONE: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt=(measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  ekf_.Q_(0,0) = dt*dt/4*noise_ax;
  ekf_.Q_(0,2) = dt/2*noise_ax;
  ekf_.Q_(1,1) = dt*dt/4*noise_ay;
  ekf_.Q_(1,3) = dt/2*noise_ay;
  ekf_.Q_(2,0) = dt/2*noise_ax;
  ekf_.Q_(2,2) = noise_ax;
  ekf_.Q_(3,1) = dt/2*noise_ay;
  ekf_.Q_(3,3) = noise_ay;  
  ekf_.Q_ *= dt*dt;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO DONE:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO DONE: Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // TODO DONE: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
