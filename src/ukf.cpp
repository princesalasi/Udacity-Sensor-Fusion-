#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>


using Eigen::MatrixXd;
using Eigen::VectorXd;
#define ESP 0.001 //Just a small number
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
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
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
  
  //Augumented State dimension
  n_aug_ = 7;

  //State dimension
  n_x_ = 5;

  // Compute lambda
  lambda_ = 3 - n_aug_;


  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    float px, py, vx, vy, v, yaw, yaw_dot = 0;
    //1. Initialization
    if (!is_initialized_)
    {
        P_ = MatrixXd::Identity(5, 5);

      if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
      {
        px = meas_package.raw_measurements_[0];
        py = meas_package.raw_measurements_[1];

        //Deal with special cases
        if (fabs(px) < ESP && fabs(py) < ESP)
          px = py = ESP;

        x_ << px, py, 0, 0, 0;
        // P_.diagonal() <<  (std_laspx_ * std_laspx_), 
        //                   (std_laspy_ * std_laspy_),
        //                   1,
        //                   1,
        //                   1;
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
      {
        // Read raw measurements from radar
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double rho_dot = meas_package.raw_measurements_[2];

        // Convert radar data from polar to catesian coordinate
        px = rho * cos(phi);
        py = rho * sin(phi);
        vx = rho_dot * cos(phi);
        vy = rho_dot * sin(phi);
        v = sqrt(vx * vx + vy * vy);
        x_ << px, py, v, 0, 0;
        // P_.diagonal() <<  (std_radr_ * std_radr_), 
        //                   (std_radr_ * std_radr_), 
        //                   (std_radrd_ * std_radrd_), 
        //                   (std_radphi_ * std_radphi_),
        //                   1;
      }

      //Initialise weight
      double weight = 0.5 / (lambda_ + n_aug_);
      weights_.resize(2 * n_aug_ +1);
      weights_.fill(weight);
      weights_(0) = lambda_ / (lambda_ + n_aug_);
      
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }

    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    //2. Prediction
    UKF::Prediction(delta_t);

    //3. Update
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
      UpdateLidar(meas_package);

    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
      UpdateRadar(meas_package);
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  MatrixXd Q = MatrixXd(2, 2);
  Q.fill(0);
  Q.diagonal() << (std_a_ * std_a_), 
                  (std_yawdd_ * std_yawdd_);

  //2.1 Augmented State Vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;

   //2.2 Calculate the augmented covariance matrix
  MatrixXd P_aug = UKF::calculateCovarianceMatrix(P_, Q);

  //2.3 Generate Sigma Points
  MatrixXd Xsig_aug = generateSigPts(x_aug, lambda_, n_aug_, P_aug);

  //2.4 Predicted Sigma points
  UKF::predictSigmaPts(Xsig_aug, delta_t);

  //2.5 Predicted mean and Covariance matrix
  UKF::predictedMean();
  UKF::predictedCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;
  int numSigPts = 2 * n_aug_ + 1;
  MatrixXd Zsig = MatrixXd(n_z, numSigPts);
  VectorXd z_pred = VectorXd(n_z); //Measurement model

  for (int i = 0; i < numSigPts; i++)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig.col(i) << px, py;
    z_pred += weights_(i) * Zsig.col(i);
    //measurement model
  }
  UKF::update(meas_package, Zsig, z_pred, n_z);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3; //State dimension Radar is 3D
  int numSigPts = 2 * n_aug_ + 1;
  VectorXd z_pred = VectorXd(n_z); //Measurement model
  //3.1 Predict mean measurement
  MatrixXd Zsig = MatrixXd(n_z, numSigPts);  // Predicted measurement mean
  z_pred.fill(0);

  //Transform sigma pts to measurement space
  for (int i = 0; i < numSigPts; i++)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double yaw_dot = Xsig_pred_(4, i);

    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rho_dot = ((px * v) + (py * v)) / rho;
    //double rho_dot = ((px * cos(yaw) * v) + (py * sin(yaw) * v)) / rho;

    Zsig.col(i) <<  rho, 
                    phi, 
                    rho_dot;
    z_pred += weights_(i) * Zsig.col(i);  //Measurement prediction
  } 
  UKF::update(meas_package, Zsig, z_pred, n_z);
}

//TOOLS 

//Generation of Sigma points
MatrixXd UKF::generateSigPts(VectorXd X_aug, double lambda, int stateDim, MatrixXd P_aug)
{
	int sigPts = 2 * n_aug_ + 1;
	MatrixXd xSig = MatrixXd(stateDim, sigPts);
	MatrixXd A = P_aug.llt().matrixL();
	float sqrtLambda = sqrt(lambda_ + n_aug_);

	xSig.col(0) = X_aug;

	for (int i = 0; i < n_aug_; i++)
	{
		xSig.col(i + 1) =  X_aug + sqrtLambda * A.col(i);
		xSig.col(i + 1 + n_aug_) = X_aug - sqrtLambda * A.col(i);
	}

	return xSig;
}

//Calculate the augmented Covariance matrix
MatrixXd UKF::calculateCovarianceMatrix(MatrixXd P, MatrixXd Q)
{
	int row = P.rows() + Q.rows();
	int col = P.cols() + Q.cols();
	MatrixXd P_c = MatrixXd(row, col);
	P_c.fill(0);
	P_c.topLeftCorner(P.rows(), P.cols()) = P;
	P_c.bottomRightCorner(Q.rows(), Q.cols()) = Q;

	return P_c;
}

//Prediction of the mean covariance matrix
void UKF::predictedCovariance()
{
  int numSigPts = n_aug_ * 2 + 1;

  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  { 
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    UKF::NormAngle(&x_diff(3)); 
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

//Prediction of Sigma points
void UKF::predictSigmaPts(Eigen::MatrixXd Xsig, float delta_t)
{
  MatrixXd Xsig_pred = MatrixXd(n_x_, Xsig.cols()); // Predicted sigma points
  Xsig_pred_.resize(n_x_, Xsig.cols());
  for (int i = 0; i < 2 * n_aug_ +1; i++)
  {
    double px = Xsig(0, i); // Position X
    double py = Xsig(1, i); // Position Y
    double v = Xsig(2, i);  // Velocity
    double yaw = Xsig(3, i);  // Yaw
    double yaw_dot = Xsig(4, i);  // Yaw rate
    double nu_a = Xsig(5, i);  // Velocity noise
    double nu_yawdd = Xsig(6, i); // Yaw noise
    double px_pred;
    double py_pred;
    double yaw_pred;
    double v_pred;

    VectorXd predictedState = VectorXd(n_x_); //Predicted state
    VectorXd nu_1 = VectorXd(n_x_);  // Noise factor 1
    VectorXd nu_2 = VectorXd(n_x_);  // Noise factor 2

    //Avoid division by zero
    if (fabs(yaw_dot) > 0.001)
    {
      px_pred = px + (v / yaw_dot) * (sin(yaw + yaw_dot * delta_t) - sin(yaw)); 
      py_pred = py + (v / yaw_dot) * (-cos(yaw + yaw_dot * delta_t) + cos(yaw));
    }
    else
    {
      px_pred = px + v * cos(yaw) * delta_t;
      py_pred = py + v * sin(yaw) * delta_t;
    }
    v_pred = v;
    yaw_pred = yaw + yaw_dot * delta_t;
    predictedState << px_pred, py_pred, v_pred, yaw_pred, yaw_dot;
      // Add noise
    nu_2 << 0.5 * (delta_t * delta_t) * cos(yaw) * nu_a,
            0.5 * (delta_t * delta_t) * sin(yaw) * nu_a,
            delta_t * nu_a,
            0.5 * (delta_t * delta_t) * nu_yawdd,
            delta_t * nu_yawdd;
    
    // Add Predicted state, noise1 and noise2 to get predicted sigma point for each column
    Xsig_pred_.col(i) = predictedState + nu_2;
  }
}

void UKF::predictedMean()
{
  x_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
}

void UKF::update(MeasurementPackage meas_package, MatrixXd Zsig, VectorXd z_pred, int n_z)
{
  MatrixXd S = MatrixXd(n_z, n_z);  //Predicted measurement covariance
  MatrixXd R = MatrixXd(n_z, n_z);  // Measurement noise covariance

  S.fill(0);
  R.fill(0);
  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
  {
    R.diagonal() << pow(std_laspx_, 2), 
                    pow(std_laspy_, 2); 
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
  {
    R.diagonal() << pow(std_radr_, 2),
                    pow(std_radphi_, 2),
                    pow(std_radrd_, 2);
  }
  //3.2 Predicted measurement covariance
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd zdiff = Zsig.col(i) - z_pred;
    NormAngle(&(zdiff(1)));
    S += weights_(i) * zdiff * zdiff.transpose();
  }
  //3.3 Add measurement noise
  S += R;


  //3.4 Compute Cross Correlation Matrix Zsig, z_pred,
  MatrixXd T_c = MatrixXd(n_x_, n_z);
  T_c.fill(0);

  for (int i = 0; i < 2 * n_aug_+ 1; i++)
  {
    VectorXd zdiff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
      NormAngle(&(zdiff(1)));
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    NormAngle(&(xdiff(3)));
    T_c = T_c + weights_(i) * xdiff * zdiff.transpose();
  }

  //3.5 Compute Kalman gain
  MatrixXd K = T_c * S.inverse();

  //3.6 Update State mean and covariance matrix
  VectorXd zdiff = meas_package.raw_measurements_ - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    NormAngle(&(zdiff(1)));
   x_ = x_ + K * zdiff;
   P_ = P_ - K * S * K.transpose();

   //Calculate NIS
   auto NIS = zdiff.transpose() * S.inverse() * zdiff; 
}


//Normalize Angle to [-PI PI]
void UKF::NormAngle(double *alpha)
{
  *alpha = fmod(*alpha + M_PI, 2 * M_PI);
  *alpha >= 0 ? (* alpha - M_PI) : (*alpha + M_PI);
  
  //while (*alpha > M_PI) *alpha -= 2. * M_PI;
  //while (*alpha < -M_PI) *alpha += 2. * M_PI;
}