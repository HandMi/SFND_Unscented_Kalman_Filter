#include "ukf.h"
#include "Eigen/Dense"
#include "measurement_package.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
inline void normalize_angle(double &angle) {
  while (angle > M_PI) {
    angle -= 2. * M_PI;
  }
  while (angle < -M_PI) {
    angle += 2. * M_PI;
  }
}
} // namespace
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimensions
  n_x_ = States::NumberOfStates;
  n_aug_ = AugmentedStates::NumberOfStates;
  n_sig_ = 2 * AugmentedStates::NumberOfStates + 1;
  n_lid_ = LidarStates::NumberOfStates;
  n_rad_ = RadarStates::NumberOfStates;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // time when the state is true, in us
  time_us_ = 0.0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values
   */

  // Initialize measurement noise covariance matrices
  R_radar_ = MatrixXd(n_rad_, n_rad_);
  R_radar_ << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0, std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(n_lid_, n_lid_);
  R_lidar_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;

  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  Xsig_pred_.fill(0);
  // Sigma point spreading parameter
  lambda_ = 3. - n_aug_;
  weights_ = VectorXd(n_sig_);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  double weight = 0.5 / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < n_sig_; ++i) {
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_) {
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 10.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 10., 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      double px = rho * std::cos(phi);
      double py = rho * std::sin(phi);
      double vx = rho_dot * std::cos(phi);
      double vy = rho_dot * std::sin(phi);
      double v = std::sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0., 0.;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  } else {
    double dt = static_cast<double>(meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      UpdateLidar(meas_package);
    }
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  // create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(AugmentedStates::NuA) = 0;
  x_aug(AugmentedStates::NuPsiDotDot) = 0;

  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(AugmentedStates::NuA, AugmentedStates::NuA) = std_a_ * std_a_;
  P_aug(AugmentedStates::NuPsiDotDot, AugmentedStates::NuPsiDotDot) = std_yawdd_ * std_yawdd_;

  // Generate Sigma Points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  // calculate square root of P
  MatrixXd L = P_aug.llt().matrixL();

  double lambda_aug = std::sqrt(lambda_ + n_aug_);

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) = x_aug + lambda_aug * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - lambda_aug * L.col(i);
  }

  // Predict sigma points
  for (int i = 0; i < Xsig_aug.cols(); i++) {
    // Extract values for better readability
    double px = Xsig_aug(AugmentedStates::Px, i);
    double py = Xsig_aug(AugmentedStates::Py, i);
    double v = Xsig_aug(AugmentedStates::V, i);
    double psi = Xsig_aug(AugmentedStates::Psi, i);
    double psi_dot = Xsig_aug(AugmentedStates::PsiDot, i);
    double nu_a = Xsig_aug(AugmentedStates::NuA, i);
    double nu_psi_dot_dot = Xsig_aug(AugmentedStates::NuPsiDotDot, i);
    // Predicted state values
    double px_predicted, py_predicted;
    // avoid division by zero
    if (fabs(psi_dot) > 0.00001) {
      px_predicted = px + v / psi_dot * (std::sin(psi + psi_dot * delta_t) - std::sin(psi));
      py_predicted = py + v / psi_dot * (std::cos(psi) - std::cos(psi + psi_dot * delta_t));
    } else {
      px_predicted = px + v * delta_t * std::cos(psi);
      py_predicted = py + v * delta_t * std::sin(psi);
    }

    double v_predicted = v;
    double psi_predicted = psi + psi_dot * delta_t;
    double psi_dot_predicted = psi_dot;

    // add noise
    px_predicted += 0.5 * nu_a * delta_t * delta_t * std::cos(psi);
    py_predicted += 0.5 * nu_a * delta_t * delta_t * std::sin(psi);
    v_predicted += nu_a * delta_t;
    psi_predicted += 0.5 * nu_psi_dot_dot * delta_t * delta_t;
    psi_dot_predicted += nu_psi_dot_dot * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(States::Px, i) = px_predicted;
    Xsig_pred_(States::Py, i) = py_predicted;
    Xsig_pred_(States::V, i) = v_predicted;
    Xsig_pred_(States::Psi, i) = psi_predicted;
    Xsig_pred_(States::PsiDot, i) = psi_dot_predicted;
  }

  // Predicted mean
  x_ = Xsig_pred_ * weights_;

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) { // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    normalize_angle(x_diff(States::Psi));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  // incoming lidar measurement
  VectorXd z = meas_package.raw_measurements_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_lid_, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_lid_);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_lid_, n_lid_);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    // measurement model
    Zsig(LidarStates::Px, i) = Xsig_pred_(LidarStates::Px, i);
    Zsig(LidarStates::Py, i) = Xsig_pred_(LidarStates::Py, i);
  }

  // mean predicted measurement
  z_pred = Zsig * weights_;

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_lidar_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_lid_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  ++lidar_count_;
  if(NIS>5.991){
    ++lidar_nis_count_;
  std::cout << "Lidar exceeds 5% NIS threshold at rate " << static_cast<double>(lidar_nis_count_)/static_cast<double>(lidar_count_) << "\n";
  }
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_rad_, n_sig_);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_rad_);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_rad_, n_rad_);

  // transform sigma points into measurement space
  for (int i = 0; i < n_sig_; ++i) {
    // measurement model
    double px = Xsig_pred_(States::Px,i);
    double py = Xsig_pred_(States::Py,i);
    double v  = Xsig_pred_(States::V,i);
    double psi = Xsig_pred_(States::Psi,i);

    double v1 = cos(psi)*v;
    double v2 = sin(psi)*v;

    // measurement model
    Zsig(RadarStates::R,i) = sqrt(px*px + py*py);
    Zsig(RadarStates::Phi,i) = atan2(py,px);
    Zsig(RadarStates::RDot,i) = (px*v1 + py*v2) / sqrt(px*px + py*py);
  }

  // mean predicted measurement
  z_pred = Zsig * weights_;

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    normalize_angle(z_diff(RadarStates::Phi));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_radar_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_rad_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; ++i) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    normalize_angle(z_diff(RadarStates::Phi));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    normalize_angle(x_diff(States::Psi));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  ++radar_count_;
  if(NIS>7.815){
    ++radar_nis_count_;
  std::cout << "Radar exceeds 5% NIS threshold at rate " << static_cast<double>(radar_nis_count_)/static_cast<double>(radar_count_) << "\n";
  }
}