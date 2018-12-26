#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO DONE: Calculate the RMSE here.
   */
  //VectorXd result;
  //return result;
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  unsigned int n = estimations.size();
  if (n == 0 || n != ground_truth.size()) {
    cout << "Zero estimate vector or incorrect vector length."
         << endl;
    return rmse;
  }
  for(unsigned int i=0; i < n; i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    rmse += (residual.array()*residual.array()).matrix();
  }
  rmse /= n;
  rmse = ((rmse.array()).sqrt()).matrix();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO DONE:
   * Calculate a Jacobian here.
   */
  //VectorXd result;
  //return result;
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float norm_p = sqrt(pow(px, 2) + pow(py, 2));
  //check division by zero
  if (norm_p < 0.0001)
      Hj << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  else
  //compute the Jacobian matrix
      Hj << px/norm_p, py/norm_p, 0, 0,
      -py/pow(norm_p, 2), px/pow(norm_p, 2), 0, 0,
      py*(vx*py - vy*px)/pow(norm_p, 3),
      px*(vy*px - vx*py)/pow(norm_p, 3),
      px/norm_p, py/norm_p;
  return Hj;
}
