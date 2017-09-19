#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_+1);

  //set weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; ++i) {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // convert radar data from polar to cartesian coordinates 
      // and initialize state
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      float pos_x = rho * cos(phi);
      float pos_y = rho * sin(phi);
      // float vel = rho_dot;
      // float vx = rho_dot * cos(phi);
      // float vy = rho_dot * sin(phi);
      // float psi = atan2(vy, vx);
      float vel = 0.0;
      float psi = 0.0;
      // we can't know the acceleration of psi angle from the first measurement.
      float psi_dot = 0.0;

      x_ << pos_x, pos_y, vel, psi, psi_dot;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  VectorXd z;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    z = VectorXd(3);
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // fill the process noise covariance matrix
  // P_ << 0.5 * dt*dt * cos(x_(3)) * std_a_*std_a_, 0, 0, 0, 0,
  //       0, 0.5 * dt*dt * sin(x_(3)) * std_a_*std_a_, 0 ,0, 0,
  //       0, 0, dt * std_a_ * std_a_, 0, 0, 
  //       0, 0, 0, 0.5 * dt*dt * std_yawdd_ * std_yawdd_, 0,
  //       0, 0, 0, 0, dt * std_yawdd_ * std_yawdd_;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0 ,0, 0,
        0, 0, 1, 0, 0, 
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;


  // generate sigma points
  
  // MatrixXd sigma_points = MatrixXd(n_x_, 2*n_x_ + 1);
  // GenerateSigmaPoints(&sigma_points);

  // print result
  // std::cout << "sigma_points = " << std::endl << sigma_points << std::endl;


  MatrixXd Xsig_aug = MatrixXd(2*n_aug_, n_aug_);
  AugmentedSigmaPoints(&Xsig_aug);
  
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);
  SigmaPointPrediction(&Xsig_aug, &Xsig_pred, dt);


  PredictMeanAndCovariance(&Xsig_pred);

  MatrixXd Zsig_points;
  VectorXd z_out = VectorXd(3);
  MatrixXd S_out = MatrixXd(3, 3);
  PredictRadarMeasurement(Xsig_pred, &Zsig_points, &z_out, &S_out);

  VectorXd x_out = VectorXd(5);
  MatrixXd P_out = MatrixXd(5, 5);
  ukf.UpdateState(Xsig_pred, &Zsig_points, z_out, S_out, &x_out, &P_out);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //define spreading parameter
  lambda_ = 3 - n_x_;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig
  MatrixXd sig1 = (sqrt(lambda_+n_x_) * A).colwise() + x_;
  MatrixXd sig2 = (sqrt(lambda_+n_x_) * A * -1).colwise() + x_;
  Xsig.col(0) = x_;
  Xsig.block(0, 1, n_x_, n_x_) = sig1;
  Xsig.block(0, 1+n_x_, n_x_, n_x_) = sig2;

  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd sig1 = (sqrt(lambda_ + n_aug_) * A_aug).colwise() + x_aug;
  MatrixXd sig2 = (sqrt(lambda_ + n_aug_) * A_aug * -1).colwise() + x_aug;
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = sig1;
  Xsig_aug.block(0, 1+n_aug_, n_aug_, n_aug_) = sig2;

  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(const MatrixXd* Xsig_in, MatrixXd* Xsig_out, const double& delta_t) {

  //create example sigma point matrix
  MatrixXd Xsig_aug = *Xsig_in;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double px_pred, py_pred;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_pred = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw) );
      py_pred = p_y + v/yawd * ( -cos(yaw + yawd*delta_t) + cos(yaw) );
    }  else {
      px_pred = p_x + v*delta_t*cos(yaw);
      py_pred = p_y + v*delta_t*sin(yaw);
    }

    double v_pred = v;
    double yaw_pred = yaw + yawd*delta_t;
    double yawd_pred = yawd;

    // add noise
    px_pred += 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_pred += 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_pred += nu_a*delta_t;
    yaw_pred += 0.5*nu_yawdd*delta_t*delta_t;
    yawd_pred += nu_yawdd*delta_t;

    //write predicted sigma points into right column
    Xsig_pred(0,i) = px_pred;
    Xsig_pred(1,i) = py_pred;
    Xsig_pred(2,i) = v_pred;
    Xsig_pred(3,i) = yaw_pred;
    Xsig_pred(4,i) = yawd_pred;
  }

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}

void UKF::PredictMeanAndCovariance(const MatrixXd* Xsig_pred) {

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    x += weights_(i) * Xsig_pred->col(i);
  }

  //predict state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    // state difference
    VectorXd x_diff = Xsig_pred->col(i) - x;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2.*M_PI;

    P += weights_(i) * x_diff * x_diff.transpose();
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  x_ = x;
  P_ = P;
}

void UKF::PredictRadarMeasurement(const MatrixXd& Xsig_pred, MatrixXd* Zsig_points, VectorXd* z_out, MatrixXd* S_out) {
  // set measurement dimension, radar can measure rho, phi, and rho_dot
  int n_z = 3;

  // create matrix for simga points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double v_x = cos(yaw)*v;
    double v_y = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                       // rho
    Zsig(1, i) = atan2(p_y, p_x);                               // phi (bearing)
    Zsig(2, i) = (p_x*v_x + p_y*v_y) / sqrt(p_x*p_x + p_y*p_y);   // rho_dot
  }

  // mean of predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    z_pred += weights_(i)*Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;

  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  *Zsig_points = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateState(const MatrixXd& Xsig_pred, const Matrix* Zsig_points, VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;


  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = *Zsig_points;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << 

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2.*M_PI;

    Tc += weights(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2.*M_PI;

  // update state mean and covariance matrix
  x += K * z_diff;
  P -= K * S * K.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}