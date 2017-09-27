#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
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

  ///* state vector: [px py vel yaw_angle yaw_rate] in SI units and rad
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

  ///* Laser measurement noise standard deviation position_x in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position_y in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  MatrixXd Xsig_aug_;

  ///* Collect NIS data from lidar and radar
  std::vector<double> lidar_NIS_record_;
  std::vector<double> radar_NIS_record_;

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

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateAugmentSigmaPoints();
  void SigmaPointPrediction(const double& delta_t);
  void PredictMeanAndCovariance();

  void PredictRadarMeasurement(MatrixXd* Zsig_points, VectorXd* z_out, MatrixXd* S_out);
  void PredictLidarMeasurement(MatrixXd* Zsig_points, VectorXd* z_out, MatrixXd* S_out);

  void UpdateState_Radar(const MatrixXd& Zsig_points, const VectorXd& z_pred, const MatrixXd& S,
                   const VectorXd& z_meas);
  void UpdateState_Lidar(const MatrixXd& Zsig_points, const VectorXd& z_pred, const MatrixXd& S,
                   const VectorXd& z_meas);

  double Calc_NIS(const VectorXd& z_diff, const MatrixXd& S);
};

#endif /* UKF_H */
