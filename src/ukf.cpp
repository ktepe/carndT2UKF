#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "ket.h"
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
  //std_a_ = 30;
  std_a_ = 3;
	//double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  
  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
	std_yawdd_ = 0.2;
	
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

  /**  TODO:Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...  */
  is_initialized_=false; 
  
  P_ << 1, 0, 0, 0, 0,
  			0, 1, 0, 0, 0,
  			0, 0, 1, 0, 0,
  			0, 0, 0, 1, 0,
  			0, 0, 0, 0, 1;
  n_x_=5;
  n_aug_=7;
  time_us_= 0;
  use_radar_=true;
  use_laser_=false; 
  

}

UKF::~UKF() {}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix. */
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  GenerateSigmaPoints(Xsig);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(Xsig, Xsig_aug);
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  MatrixXd Xsig_pred = MatrixXd (n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0);
  SigmaPointPrediction(Xsig_aug, Xsig_pred, delta_t); 
	PredictMeanAndCovariance(Xsig_pred); 
	Xsig_pred_=Xsig_pred;
#if ket_debug
	cout << " end of prediction" <<endl;
	cout<< Xsig_pred_ <<endl;
#endif
  
}
void UKF::GenerateSigmaPoints(MatrixXd& Xsig) {

	double lambda = 3 - n_x_;
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  //your code goes here 
  Xsig.col(0)=x_;
	double sqrt_lam_nx=sqrt(lambda+n_x_);
  for(int i=1; i<n_x_+1; i++){
      Xsig.col(i)=x_+sqrt_lam_nx*A.col(i-1);
  }
  for(int i=n_x_+1; i<2*n_x_+1; i++){
      Xsig.col(i)=x_-sqrt_lam_nx*A.col(i-(n_x_+1));
  }
#if ket_debug
	cout << "Generate Sigma Points "<<endl;
  cout << Xsig << endl;
#endif	
	return;
}
void UKF::AugmentedSigmaPoints(MatrixXd &Xsig, MatrixXd  &Xsig_aug) {

  //Process noise standard deviation longitudinal acceleration in m/s^2
  
  //define spreading parameter
  double lambda = 3 - n_aug_;
 
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

 
  //create augmented mean state
  x_aug.head(n_x_)=x_;
	x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(n_x_,n_x_)=std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1)=std_yawdd_*std_yawdd_;
 
  //create square root matrix
    MatrixXd L=P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }
 
  //print result
#if ket_debug
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
#endif
 return; 

}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug,  MatrixXd &Xsig_pred, double delta_t) {

#if ket_debug
  cout << "in Sigma point prediction " << Xsig_aug.size()<< endl;
  cout << Xsig_aug <<endl;
#endif
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

#if ket_debug
  //print result
  cout << "Xsig_pred = " << endl << Xsig_pred << endl;
#endif
  return;
}

void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred) {

  //set state dimension
	//later replace with n_x_ and n_aug_
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
 	/* 
  //create vector for predicted state
  VectorXd x = x_

  //create covariance matrix for prediction
  MatrixXd P = P_
  */
  // set weights
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }

#if ket_debug
  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;
#endif
}





/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
  /**  TODO:Complete this function! Make sure you switch between lidar and radar
  measurements.  */
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF initialization: " << endl;
    x_ << 1, 1, 1, 1, 1;
    double px=0.0;
    double py=0.0;
    double v=0.0;
    double yaw=0.0;
    double yaw_rate=0.0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /** Convert radar from polar to cartesian coordinates and initialize state.*/
   		cout << "Radar initialization: " << endl;
    
      double rho=meas_package.raw_measurements_[0];
      double phi=meas_package.raw_measurements_[1];
     
      px=rho*cos(phi);
      py=rho*sin(phi); 

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
     	cout<< "LASER initialization" << endl;
   	
   		px=meas_package.raw_measurements_[0];
   		py=meas_package.raw_measurements_[1];
   		      
    }
    // if px and py are too small, problems with divide with zero etc.
    if (fabs(px) < 0.0001){
    	px=0.1;
    }
    
    if (fabs(py)< 0.0001){
    	py=0.1;
    }
	
		x_ << px, py, v, yaw, yaw_rate;
		time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
#if ket_debug
	cout<<"Initialization ends  :" << x_.transpose() << endl;
#endif   
    return;
  }
  //dt - expressed in seconds
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	
	time_us_ = meas_package.timestamp_;
#if ket_debug
	cout<<"Prediction called"<<endl;

#endif   
  Prediction(delta_t);
#if ket_debug
	cout<<"Prediction ended"<<endl;
#endif    
 
 
 /*****************************************************************************
   *  Update
   ****************************************************************************/

 if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)){
    // Radar updates
#if ket_debug	
		cout << "Radar update" << endl;
#endif
		if(use_laser_){
			return;
		}
  
	  UpdateRadar(meas_package);
  	
  } 

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER)){
	return;
	}
// print the output
  cout << "x_ = " << x_.transpose() << endl;
  cout << "P_ = " << P_ << endl;
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
#if ket_debug
	cout<< "in Update radar " << endl;
#endif 
  int n_z=3;
  
  VectorXd z_pred = VectorXd(n_z);
 
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
 
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  PredictRadarMeasurement(z_pred, S, Zsig);
  //measurements
  VectorXd z = VectorXd(n_z);
  z(0)=meas_package.raw_measurements_[0];
  z(1)=meas_package.raw_measurements_[1];
  z(2)=meas_package.raw_measurements_[2];
  
  
  UpdateState(z, z_pred, S, Zsig);
  
  return;
}

//rest is copied from my lecture programming assingments

void UKF::PredictRadarMeasurement(VectorXd &z_pred, MatrixXd  &S, MatrixXd &Zsig) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //radar measurement noise standard deviation radius in m
//  double std_radr = 0.3;
	double std_radr = std_radr_; //
  //radar measurement noise standard deviation angle in rad
//  double std_radphi = 0.0175;
  double std_radphi = std_radphi_;
  //radar measurement noise standard deviation radius change in m/s
  //double std_radrd = 0.1;
	double std_radrd = std_radrd_;
  //create example matrix with predicted sigma points
 
#if ket_debug
	cout<< "in Predict Radar Measurement; and Xsig_pred_" << endl;
	cout<< Xsig_pred_<<endl;
#endif 

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug+1; i++) {
      double px=Xsig_pred_(0,i);
      double py=Xsig_pred_(1,i);
      double v=Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
      double rho=sqrt(px*px+py*py);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    double rho_dot=(px*v1+py*v2)/rho;
    double phi=atan2(py, px);

     /* 
     double phi=0;
      if (fabs(px)<0.0001){
          std::cout<<"error in atan"<<std::endl;
      } else {
          phi=atan2(py, px);
      }
      
      double rho_dot=0;
     if (rho>0.0001){
         rho_dot=(px*v1+py*v2)/rho;
     }
     */
     // measurement model
    Zsig(0,i) = sqrt(px*px + py*py);                        //r
    Zsig(1,i) = atan2(py,px);                                 //phi
    Zsig(2,i) = (px*v1 + py*v2 ) / sqrt(px*px + py*py);   //r_dot
      
      
  }
#if ket_debug
	cout<< "in Predict Radar Measurement; and Zsig" << endl;
	cout<< Zsig<<endl;
#endif 

  //calculate mean predicted measurement
  for(int i=0; i<2*n_aug+1; i++){
      z_pred+=weights(i)*Zsig.col(i);
  }
  //calculate measurement covariance matrix S
  S.fill(0.0);
#if ket_debug
	cout<< "in Predict Radar Measurement; and z_pred" << endl;
	cout<< z_pred<<endl;
#endif  
  
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
    while (x_diff(1)<-M_PI) x_diff(1)+=2.*M_PI;

    S = S + weights(i) * x_diff * x_diff.transpose() ;
  }
#if ket_debug
	cout<< "in Predict Radar Measurement; and S" << endl;
	cout<< S<<endl;
#endif 

//add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr*std_radr, 0, 0,
          0, std_radphi*std_radphi, 0,
          0, 0,std_radrd*std_radrd;
  S = S + R;

#if ket_debug
  //print result
  cout << "z_pred: " << endl << z_pred << endl;
  cout << "S: " << endl << S << endl;
#endif
 return;
 
}

void UKF::UpdateState(VectorXd &z, VectorXd &z_pred, MatrixXd &S, MatrixXd &Zsig) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
   double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  //create example matrix with predicted sigma points
//  MatrixXd Xsig_pred = Xsig_pred_;
  

  //create example vector for predicted state mean
//  VectorXd x = x_;
  
  //create example matrix for predicted state covariance
//  MatrixXd P = P_;
  
 

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

//std part begin
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
#if ket_debug
  //print result
  cout << "Updated state x: " << endl << x_ << endl;
  cout << "Updated state covariance P: " << endl << P_ << endl;
#endif

  //write result
//  x_ = x;
//  P_ = P;
	return;
	
}
