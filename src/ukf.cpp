#include "ukf.h"
#include "tools.h"
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
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .1;

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

  // Number of State dimenstions
  n_x_ = 5;

  // Number of augmented dimensions
  n_aug_ = 7;

  // Initilize Predicted sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);
  Xsig_pred_.fill(0.0); 

  // Spreading parameter
  lambda_ = 3 - n_aug_;

  // Weight vecotr
  weights_ = VectorXd(2*n_aug_+1);

  weights_[0] = lambda_/(lambda_+n_aug_);

  for (int i=1; i<2*n_aug_+1; i++){
    weights_[i] = 1/(2*(lambda_+n_aug_));
  }

  // Set initialization to false
  is_initialized_ = false;

  // Set initial time to 0
  time_us_ = 0.0;
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
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
  // Initialize states from first measurement

  //cout << "measurement type" << meas_package.sensor_type_ << endl;
  //cout << "measurement values" << meas_package.raw_measurements_ << endl;

  

  if (!is_initialized_){

    P_<<1,0,0,0,0,
      0,1,0,0,0,
      0,0,1,0,0,
      0,0,0,1,0,
      0,0,0,0,1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd meas = meas_package.raw_measurements_;
      x_ << meas[0]*cos(meas[1]), meas[0]*sin(meas[1]),0,0,0;
      time_us_ = meas_package.timestamp_;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;
      time_us_ = meas_package.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;

  }

  // Prediction Step

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000;

  Prediction(delta_t);

  time_us_ = meas_package.timestamp_;

  // Ingnoring all-zero measurements (for sample-data-2)

  if (meas_package.raw_measurements_[0]==0 && meas_package.raw_measurements_[1]==0){
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && meas_package.raw_measurements_[2]==0) {
      time_us_ = meas_package.timestamp_;
      // print the output
      cout<< "All RADAR measurements zero"<< endl;
      //cout << "x_ = " << x_ << endl;
      //cout << "P_ = " << P_ << endl;
      return;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      time_us_ = meas_package.timestamp_;
      cout<< "All LASER measurements zero"<< endl;
      //cout << "x_ = " << x_ << endl;
      //cout << "P_ = " << P_ << endl;
      return;
    }

  }

  // Update step: Check the measurement type and call the correct update

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
  
  //cout << "P_ = " << P_ << endl;
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

  
  /******Generate augmented sigma points******/

  // Augmented x
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug.tail(2) <<0,0;

  // Augmented P
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);

  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  cout<<"P_aug:"<<P_aug<<endl;
  // Squareroot matrix
  MatrixXd A = P_aug.llt().matrixL();

  // Augmented Sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);

  Xsig_aug.col(0) = x_aug;

  for(int i=1; i<=n_aug_; i++){
    VectorXd col = sqrt(lambda_+n_aug_)*A.col(i-1);
    Xsig_aug.col(i) = x_aug + col;
    Xsig_aug.col(i+n_aug_) = x_aug - col;
  }

  /******Predicting New Sigma points******/

  VectorXd transition = VectorXd(n_x_);
  VectorXd noise = VectorXd(n_x_);

  MatrixXd x_k = Xsig_aug.topRows(5);

  double delta_t_2 = delta_t*delta_t;

  for(int i=0; i<2*n_aug_+1; i++){
    // Assign names to values in Xsig_aug
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);


    // Calculation noise vector
    noise[0] = 0.5*delta_t_2*cos(yaw)*nu_a;
    noise[1] = 0.5*delta_t_2*sin(yaw)*nu_a;
    noise[2] = delta_t*nu_a;
    noise[3] = 0.5*delta_t_2*nu_yawdd;
    noise[0] = delta_t*nu_yawdd;

    // Calculating transition vector
    transition.fill(0.0);

    if (yawd==0){
      transition[0] = v*cos(yaw)*delta_t;
      transition[1] = v*sin(yaw)*delta_t;
    }
    else{
      transition[0] = (sin(yaw+yawd*delta_t)-sin(yaw))*v/yawd;
      transition[1] = (-cos(yaw+yawd*delta_t)+cos(yaw))*v/yawd;
      transition[3] = yawd*delta_t;
    }

    Xsig_pred_.col(i) = x_k.col(i) + transition + noise;
  }

  /******Predict Weighted Mean and Covariance of Sigma points******/

  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);

  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);

  for(int i=0; i<2*n_aug_+1; i++){
    x_pred += weights_[i]*Xsig_pred_.col(i);
  }

  for(int i=0; i<2*n_aug_+1; i++){
    MatrixXd sub = Xsig_pred_.col(i) - x_pred;
    P_pred += weights_[i]*sub*sub.transpose();
  }

  x_ = x_pred;
  P_ = P_pred;
  //cout<<"P_"<<P_<<endl;

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
  
  /******Predic LASER Measurement******/

  // Laser measurement dimension
  int n_laser = 2;

  // Matrix for sigma points converted to measurement space
  MatrixXd Zsig = MatrixXd(n_laser, 2*n_aug_+1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_laser);
  z_pred.fill(0.0);

  // Measurement covariance matrix 
  MatrixXd S = MatrixXd(n_laser, n_laser);
  S.fill(0.0);

  // Measurement noise covariance
  MatrixXd R = MatrixXd(n_laser, n_laser);
  R.fill(0.0);
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  // Convert Sigma points to measurements space and calc mean
  Zsig = Xsig_pred_.topRows(n_laser);

  for(int i=0; i<2*n_aug_+1; i++){
    z_pred += weights_[i]*Zsig.col(i);
  }

  // Calculate measurement covariance matrix
  S += R;
  for(int i=0; i<2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_[i]*z_diff*z_diff.transpose();
  }

  /******Laser measurement update******/

  // Current measurement
  VectorXd z = meas_package.raw_measurements_;

  // Cross correlation matrix
  MatrixXd T = MatrixXd(n_x_, n_laser);
  T.fill(0.0);

  // Kalman Gain
  MatrixXd K;

  // Calculate cross correlation
  for (int i=0; i<2*n_aug_+1; i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      while (x_diff(3)>M_PI) x_diff(3) -= 2*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;
      
      T += weights_[i]*x_diff*z_diff.transpose();
  }

  // Calculate Kalman gain
  MatrixXd Si = S.inverse();
  K = T*Si;

  // Update state x_ and covariance P_
  VectorXd meas_diff = z - z_pred;
  x_ += K*meas_diff;
  P_ -= K*S*K.transpose();

  /******Calculate NIS******/
  NIS_laser_ = meas_diff.transpose()*Si*meas_diff;

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

  /******Predic RADAR Measurement******/

  // Radar measurement dimension
  int n_radar = 3;

  // Matrix for sigma points converted to measurement space
  MatrixXd Zsig = MatrixXd(n_radar, 2*n_aug_+1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_radar);
  z_pred.fill(0.0);

  // Measurement covariance matrix 
  MatrixXd S = MatrixXd(n_radar, n_radar);
  S.fill(0.0);

  // Measurement noise covariance
  MatrixXd R = MatrixXd(n_radar, n_radar);
  R.fill(0.0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  // Convert Sigma points to measurements space and calc mean
  for(int i=0; i<2*n_aug_+1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double mag = sqrt(p_x*p_x + p_y*p_y);
    if (mag==0) return;

    Zsig(0,i) = mag;
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw))/mag;

    z_pred += weights_[i]*Zsig.col(i);
  }

  // Calculate measurement covariance matrix
  S += R;
  for(int i=0; i<2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)>M_PI) z_diff(1) -= 2*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2*M_PI;

    S += weights_[i]*z_diff*z_diff.transpose();
  }

  /******Radar measurement update******/

  // Current measurement
  VectorXd z = meas_package.raw_measurements_;

  // Cross correlation matrix
  MatrixXd T = MatrixXd(n_x_, n_radar);
  T.fill(0.0);

  // Kalman Gain
  MatrixXd K;

  // Calculate cross correlation
  for (int i=0; i<2*n_aug_+1; i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      while (x_diff(3)>M_PI) x_diff(3) -= 2*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;
      
      while (z_diff(1)>M_PI) z_diff(1) -= 2*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1) += 2*M_PI;
      
      T += weights_[i]*x_diff*z_diff.transpose();
  }

  // Calculate Kalman gain
  MatrixXd Si = S.inverse();
  K = T*Si;

  // Update state x_ and covariance P_
  VectorXd meas_diff = z - z_pred;
  x_ += K*meas_diff;
  P_ -= K*S*K.transpose();

  /******Calculate NIS******/
  NIS_radar_ = meas_diff.transpose()*Si*meas_diff;


}
