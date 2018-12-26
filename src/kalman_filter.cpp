#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO DONE: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO DONE: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateState(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO DONE: update the state by using Extended Kalman Filter equations
   */

  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float norm_p = sqrt(pow(px, 2) + pow(py, 2));

  // calculate z_pred = h(x)
  double rho, phi, rho_dot;
  rho = norm_p;
  if (norm_p < 0.0001)
    rho_dot = 0;
  else
    rho_dot = (px*vx+py*vy)/norm_p;
  if (abs(px) < 0.0001)
    if (px*py > 0)
      phi = M_PI/2;
    else if (px*py < 0)
      phi = -M_PI/2;
    else
      phi = 0;
  else // returns principal angle value so between -pi and pi
    phi = atan2(py, px);
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  while (y(1) > M_PI/2 || y(1) <= -M_PI/2)
    if (y(1) > M_PI/2)
      y(1) -= M_PI;
    else
      y(1) += M_PI;
  UpdateState(y);
}
  
void KalmanFilter::UpdateState(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}
