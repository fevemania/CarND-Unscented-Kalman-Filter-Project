#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// set measurement dimension, lidar can measure px, py
#define N_LIDAR 2
// set measurement dimension, radar can measure rho, phi, and rho_dot
#define N_RADAR 3

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  n_x_ = 5;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3.45;

  // Laser measurement noise standard deviation position_x in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position_y in m
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
  
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);

  // fill the inital state covariance matrix
  P_ << 1, 0, 0, 0, 0,
      0, 1, 0 ,0, 0,
      0, 0, 1, 0, 0, 
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;

  //set weights
  weights_ = VectorXd(2*n_aug_+1);
  
  weights_(0) = lambda_ / (lambda_+n_aug_);
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

      float v = 0.0;
      float psi = 0.0;
      // we can't know the acceleration of psi angle from the first measurement.
      float psi_dot = 0.0;

      x_ << pos_x, pos_y, v, psi, psi_dot;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Prediction Stage
  Prediction(dt);

  // Update Stage -- Radar data
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  // write NIS data into file
  // 
  // if (lidar_NIS_record_.size() == 249) {
  //   ofstream out_file_lidar;
  //   out_file_lidar.open ("lidar_NIS_record.txt");
  //   for (int i = 0; i < lidar_NIS_record_.size(); ++i) {
  //     out_file_lidar << i+1 << " " << lidar_NIS_record_[i] << " " << 5.991 << endl;
  //   }
  //   out_file_lidar.close();
  // }
  // if (radar_NIS_record_.size() == 249) {
  //   ofstream out_file_radar;
  //   out_file_radar.open ("radar_NIS_record.txt");
  //   for (int i = 0; i < radar_NIS_record_.size(); ++i) {
  //     out_file_radar << i+1 << " " << radar_NIS_record_[i] << " " << 7.815 << endl;
  //   }
  //   out_file_radar.close();    
  // }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /*
  Estimate the object's location.
  */
 
  // generate augmented sigma points
  GenerateAugmentSigmaPoints();

  // insert augment sigma points into non-linear process model 
  // to generate predicted sigma points and store into Xsig_pred_.
  SigmaPointPrediction(delta_t);

  // use predicted simga points to calculate the mean and covariance of 
  // the prediction distribution.
  // 
  // Input: Xsig_pred_, Output: predicted x_ and P_
  PredictMeanAndCovariance();
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
 
  ///* radar data dimension. [rho, theta, rho_dot].
  MatrixXd Zsig_points;
  VectorXd z_pred = VectorXd(N_RADAR);
  MatrixXd S = MatrixXd(N_RADAR, N_RADAR);

  
  // Generate sigma points of measurement state from data member Xsig_pred_
  // and predict the the mean vector and covariance matrix of the measurement state 
  PredictRadarMeasurement(&Zsig_points, &z_pred, &S);
  
  // Update the state mean vector and state covariance matrix of object of interest with radar data
  VectorXd z_meas = VectorXd(N_RADAR);
  z_meas << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  UpdateState_Radar(Zsig_points, z_pred, S, z_meas);

  VectorXd z_diff = z_meas - z_pred;
  double radar_NIS = Calc_NIS(z_diff, S);
  radar_NIS_record_.push_back(radar_NIS);
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

  // Predict the the mean vector and covariance matrix of the measurement state 
  MatrixXd Zsig_points;
  VectorXd z_pred = VectorXd(N_LIDAR);
  MatrixXd S = MatrixXd(N_LIDAR, N_LIDAR);
  PredictLidarMeasurement(&Zsig_points, &z_pred, &S);

  // Update the state mean vector and state covariance matrix of object of interest with radar data
  VectorXd z_meas = VectorXd(N_LIDAR);
  z_meas << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];


  UpdateState_Lidar(Zsig_points, z_pred, S, z_meas);

  VectorXd z_diff = z_meas - z_pred;
  double lidar_NIS = Calc_NIS(z_diff, S);
  lidar_NIS_record_.push_back(lidar_NIS);
}

void UKF::GenerateAugmentSigmaPoints() {

  //create augmented state mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance matrix 
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix of augmented state

  // fill mean vector
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0; // longitudinal acceleration noise has mean 0
  x_aug(6) = 0; //          yaw acceleration noise has mean 0

  // fill covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  // insert process noise covariance matrix Q at the bottom right of P_aug_
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // calculate square root matrix of the P_aug
  MatrixXd A_aug = P_aug.llt().matrixL();

  // calculate augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  MatrixXd sig1 = (sqrt(lambda_ + n_aug_) * A_aug).colwise() + x_aug;
  MatrixXd sig2 = (sqrt(lambda_ + n_aug_) * A_aug * -1).colwise() + x_aug;  
  Xsig_aug_.block(0, 1, n_aug_, n_aug_) = sig1;
  Xsig_aug_.block(0, 1+n_aug_, n_aug_, n_aug_) = sig2;

  // print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
}

/*
 * Put sigma points of augmented state into process model to generate corresponding 
 * simga points of predicted state, and store the coorepsonding simga points into   
 * data member Xsig_pred_
 *
 * @param Xsig_in The Xsig_aug calculated from state mean and covariance at k step
 *                by GenerateAugmentSigmaPoints function.
 * @param delta_t The time difference from k step to k+1 step
 */
void UKF::SigmaPointPrediction(const double& delta_t) {

  //predict sigma points
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    double px_pred, py_pred;

    // add deterministic part
    // careful for the situation when yawd is 0 to avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_pred = px + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw) );
      py_pred = py + v/yawd * ( -cos(yaw + yawd*delta_t) + cos(yaw) );
    }  else {
      px_pred = px + v*delta_t * cos(yaw);
      py_pred = py + v*delta_t * sin(yaw);
    }

    double v_pred = v;
    double yaw_pred = yaw + yawd*delta_t;
    double yawd_pred = yawd;

    // add noise part
    double delta_t_2 = delta_t*delta_t;
    
    px_pred += 0.5*nu_a*delta_t_2 * cos(yaw);
    py_pred += 0.5*nu_a*delta_t_2 * sin(yaw);
    v_pred += nu_a*delta_t;
    yaw_pred += 0.5*nu_yawdd*delta_t_2;
    yawd_pred += nu_yawdd*delta_t;

    //write predicted sigma points into right column
    Xsig_pred_(0,i) = px_pred;
    Xsig_pred_(1,i) = py_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;
    Xsig_pred_(4,i) = yawd_pred;
  }

  //print result
  // std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
}

/*
 * use predicted simga points to calculate the mean and covariance of 
 * the prediction distribution, and store the result into data member x_ and P_
 */
void UKF::PredictMeanAndCovariance() {

  //create mean vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for predicted state
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    x += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2.*M_PI;

    P += weights_(i) * x_diff * x_diff.transpose();
  }

  //print result
  // std::cout << "Predicted state" << std::endl;
  // std::cout << x << std::endl;
  // std::cout << "Predicted covariance matrix" << std::endl;
  // std::cout << P << std::endl;

  //write result
  x_ = x;
  P_ = P;
}

void UKF::PredictRadarMeasurement(MatrixXd* Zsig_points, VectorXd* z_out, MatrixXd* S_out) {
  

  // create matrix for simga points in measurement space
  MatrixXd Zsig = MatrixXd(N_RADAR, 2*n_aug_+1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v_x = cos(yaw)*v;
    double v_y = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                         // rho
    Zsig(1, i) = atan2(p_y, p_x);                                 // phi (bearing)
    Zsig(2, i) = (p_x*v_x + p_y*v_y) / sqrt(p_x*p_x + p_y*p_y);   // rho_dot
  }

  // predicted measurement mean vector
  VectorXd z_pred = VectorXd(N_RADAR);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    z_pred += weights_(i)*Zsig.col(i);
  }

  // predicted measurement covariance matrix S  
  MatrixXd S = MatrixXd(N_RADAR, N_RADAR);
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
  MatrixXd R = MatrixXd(N_RADAR, N_RADAR);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  S += R;

  // std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  // std::cout << "S: " << std::endl << S << std::endl;

  *Zsig_points = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::PredictLidarMeasurement(MatrixXd* Zsig_points, VectorXd* z_out, MatrixXd* S_out) {
  
  // create matrix for simga points in measurement space
  MatrixXd Zsig = MatrixXd(N_LIDAR, 2*n_aug_+1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  // mean of predicted measurement
  VectorXd z_pred = VectorXd(N_LIDAR);

  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    z_pred += weights_(i)*Zsig.col(i);
  }

  // mcovariance matrix S of predicted measurement 
  MatrixXd S = MatrixXd(N_LIDAR, N_LIDAR);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(N_LIDAR, N_LIDAR);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  S += R;

  // std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  // std::cout << "S: " << std::endl << S << std::endl;

  *Zsig_points = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateState_Radar(const MatrixXd& Zsig_points, const VectorXd& z_pred, const MatrixXd& S, const VectorXd& z_meas) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, N_RADAR);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {

    // residual
    VectorXd z_diff = Zsig_points.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z_meas - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2.*M_PI;

  // update state mean and state covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  //print result
  // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  // std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
}

void UKF::UpdateState_Lidar(const MatrixXd& Zsig_points, const VectorXd& z_pred, const MatrixXd& S, const VectorXd& z_meas) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, N_LIDAR);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; ++i) {

    // residual
    VectorXd z_diff = Zsig_points.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z_meas - z_pred;

  // update state mean and state covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  //print result
  // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  // std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
}

double UKF::Calc_NIS(const VectorXd& z_diff, const MatrixXd& S) {
  return z_diff.transpose() * S.inverse() * z_diff;
}