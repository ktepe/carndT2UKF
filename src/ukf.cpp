#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
//#include <fstream>
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
  //use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  //use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_=1.0 ;
  // Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 1.0;
  
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
  
  P_ << 0.5, 0, 0, 0, 0,
  			0, 0.5, 0, 0, 0,
  			0, 0, 10, 0, 0,
  			0, 0, 0, 10, 0,
  			0, 0, 0, 0, 0.5;
  
  n_x_=5;
  n_aug_=7;
  time_us_= 0;
  use_radar_=true;
  use_laser_=true; 
  //Lidar updates
  Hlidar_=MatrixXd(2,5);
  Rlidar_=MatrixXd(2,2);
  
  Hlidar_ << 1, 0, 0, 0, 0,
  					0, 1, 0, 0, 0;
  
  Rlidar_ << 0.0225, 0,
        			0, 0.0225;
        			
  NISradar_=VectorXd(252);
 	NISradar_counter_=0;
 	NISradar_.fill(0.0);

  NISlidar_=VectorXd(252);
  NISlidar_counter_=0;
  NISlidar_.fill(0.0);
	
	
 	NISradar_file_.open ("NISradar_data.txt");
 	NISlidar_file_.open ("NISlidar_data.txt");
	
	epoch_counter_=0;

}

UKF::~UKF() {}

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
      double rho_dot=meas_package.raw_measurements_[2];
     
      px=rho*cos(phi);
      py=rho*sin(phi); 
      v=rho_dot;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
     	cout<< "LASER initialization" << endl;
   	
   		px=meas_package.raw_measurements_[0];
   		py=meas_package.raw_measurements_[1];
   		      
    }
 	
 		// very small px and py are can cause problems.
		if ( fabs(px) < 0.001 || fabs(py) < 0.001 ) // <--- values are flexible
		{
  	 	px = 0.01;  //  <--- values are flexible
  	 	py = 0.01;  //  <--- choose to satisfy rmse requirements
		}	
	
		x_ << px, py, v, yaw, yaw_rate;
		time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
#if ket_debug
	cout<<"UKF::ProcessMeasurement--Initialization ends  :" << x_.transpose() << endl;
#endif   

    return;
  }
  //dt - expressed in seconds
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	
	time_us_ = meas_package.timestamp_;

#if ket_debug
	cout<<"UKF::ProcessMeasurement--Prediction called"<<endl;
#endif   

  //restart of the simulator needs reinitialization 
  //negative time confuses sigma predictions for reruns.
  if (delta_t < 0){
  	is_initialized_ = false;
 		P_ << 1, 0, 0, 0, 0,
  				0, 1, 0, 0, 0,
  				0, 0, 1, 0, 0,
  				0, 0, 0, 1, 0,
  				0, 0, 0, 0, 1;
  				
	  return;
  }

  //suggested by reviewer to remove some timing instabilities. 
  while (delta_t > 0.2)  // <--- value is flexible
  {
      double step = 0.1;  // <--- value is flexible
      Prediction(step);
      delta_t -= step;
  }
  Prediction(delta_t);
  	
#if ket_debug
	cout<<"UKF::ProcessMeasurement--Prediction ended"<<endl;
#endif    
 
 
 /*****************************************************************************
   *  Update
   ****************************************************************************/

 if ((meas_package.sensor_type_ == MeasurementPackage::RADAR)&&(use_radar_)){
    // Radar updates
#if ket_debug_	
		cout << "UKF::ProcessMeasurement--Radar update" << endl;
#endif
  
	  UpdateRadar(meas_package);
  	
  } 

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_)){
#if ket_debug_	
		cout << "UKF::ProcessMeasurement--Lidar update" << endl;
#endif	
	
		UpdateLidar(meas_package);
		
	}
// print the output
  cout << "UKF::ProcessMeasurement--epoch: " << epoch_counter_++<< endl;
  cout<< " x_ "<< x_.transpose() << endl;
  cout << "UKF::ProcessMeasurement--P_ = " << endl << P_ << endl;
 
	
	if ((NISradar_counter_) == 201){
		cout << "writing NISradar to the file" << endl;
	  NISradar_file_ << NISradar_.head(200).transpose()<< endl; 
  	NISradar_file_.close();
	 }
	 
	if ((NISlidar_counter_) == 201){
		cout << "writing NISlidar to the file" << endl;
	  NISlidar_file_ << NISlidar_.head(200).transpose()<< endl; 
  	NISlidar_file_.close();
	}
}

//rest is taken from my UKF lecture programming assignments

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

  MatrixXd Xsig_pred = MatrixXd (n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0);
  SigmaPointPrediction(Xsig_aug, Xsig_pred, delta_t); 

	PredictMeanAndCovariance(Xsig_pred); 
	Xsig_pred_=Xsig_pred;
#if ket_debug
	cout << "UKF::Prediction" <<endl;
	cout<<"x_ "<< x_.transpose() <<endl;
	cout<<"P_ " << P_ <<endl;
#endif 
}

void UKF::GenerateSigmaPoints(MatrixXd& Xsig) {

	double lambda = 3 - n_x_;
	double sqrt_lam_nx=sqrt(lambda+n_x_);

  //calculate square root of P_
  MatrixXd A = P_.llt().matrixL();
  //your code goes here 
  Xsig.col(0)=x_;

  //set remaining sigma points
  Xsig.fill(0.0);
  
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda+n_x_) * A.col(i);
    Xsig(3,i+1)=PhiNorm(Xsig(3,i+1));
    Xsig(3,i+1+n_x_)=PhiNorm(Xsig(3,i+1+n_x_));
  }

#if ket_debug
	cout << "UKF::GenerateSigmaPoints"<<endl;
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
	x_aug(n_x_)   = 0.0;
  x_aug(n_x_+1) = 0.0;
  
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
  cout << "UKF::AugmentedSigmaPoints "<<endl;
  cout << "Xsig_aug = " << endl << Xsig_aug << endl;
#endif
 return; 

}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug,  MatrixXd &Xsig_pred, double delta_t) {

 //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++)
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
#if ket_debug
cout<< "UKF::SigmaPointPrediction: delta_t, p_x, p_y, v, yaw, yawd " << delta_t<< ", "<< p_x << ", "<< p_y << ", " << v<< ", "<< yaw << ", " << yawd<< endl;
cout<< "UKF::SigmaPointPrediction: px_, py_, v_p, yaw_p, yawd_p " << px_p << ", "<< py_p << ", " << v_p<< ", "<< yaw_p << ", " << yawd_p<< endl;
#endif

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
  cout << "UKF::SigmaPointPrediction"<<endl;
  cout << "Xsig_pred = " << endl << Xsig_pred << endl;
#endif

}


void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred) {

  //set state dimension
  int n_x = n_x_;

  //set augmented dimension
  int n_aug = n_aug_;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  
  //set weights
  double weight_0;
  
  weight_0=1/(2*(lambda+n_aug));
  weights(0)=2*lambda*weight_0;
  for (int i=1; i<2*n_aug+1; i++){
    weights(i)=weight_0;
  }

  //predict state mean
  x_.fill(0.0);
  for(int i=0; i<2*n_aug+1; i++){
      x_+=weights(i)*Xsig_pred.col(i);
 }
 x_(3)=PhiNorm(x_(3));
 
  //predict state covariance matrix
  
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
  	x_diff(3)=PhiNorm(x_diff(3));	
  
    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }
 
#if ket_debug
  //print result
  cout << "UKF::PredictMeanAndCovariance" << endl;
  cout << "prediction: x_ " << x_.transpose() << endl;
  cout << "Predicted covariance matrix" << endl;
  cout << P_ << endl;
#endif
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
	double std_radr = std_radr_; //
  //radar measurement noise standard deviation angle in rad
  double std_radphi = std_radphi_;
  //radar measurement noise standard deviation radius change in m/s
	double std_radrd = std_radrd_;
  //create example matrix with predicted sigma points

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug+1; i++) {
      double px=Xsig_pred_(0,i);
      double py=Xsig_pred_(1,i);
      double v=Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
  		//suggested by reviewer.
  		if (fabs(px) < 0.001 || fabs(py) < 0.001){
	  		cout<< "UKF::PredictRadarMeasurement possible zero division problems px, py " << px << ", "<< py << endl;
  			px=0.01;
  			py=0.01;

  		}
      
      double rho=sqrt(px*px+py*py);

	    double v1 = cos(yaw)*v;
  	  double v2 = sin(yaw)*v;
  	  //double rho_dot=0.001;
  	  //since we took care px, py before we do not expect zero division here
  	  double rho_dot=(px*v1+py*v2)/rho;
  	  	
  	  double phi=atan2(py, px);
#if ket_debug
	cout<< "UKF::PredictRadarMeasurement: px, py, rho, phi, rhodot " << px << ", "<< py << ", " <<rho << ", "<< phi << ", " <<	 	rho_dot<< endl;
#endif

  	  Zsig(0,i) = sqrt(px*px + py*py);                        //r
  	  Zsig(1,i) = atan2(py,px);                                 //phi
  	  Zsig(2,i) = (px*v1 + py*v2 ) / sqrt(px*px + py*py);   //r_dot
      
      
  	}

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i=0; i<2*n_aug+1; i++){
      z_pred+=weights(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  S.fill(0.0);
#if ket_debug
	cout<< "UKF::PredictRadarMeasurement; and z_pred" << endl;
	cout<< z_pred.transpose()<<endl;
#endif  
	z_pred(1)=PhiNorm(z_pred(1));
  
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
  	z_diff(1)=PhiNorm(z_diff(1));	

    S = S + weights(i) * z_diff * z_diff.transpose() ;
  }
#if ket_debug
	cout<< "UKF::PredictRadarMeasurement; and S" << endl;
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
  cout << "UKF::PredictRadarMeasurement: " << endl << z_pred << endl;
  cout << "UKF::PredictRadarMeasurement: " << endl << S << endl;
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

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

	//std part begin
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    z_diff(1)=PhiNorm(z_diff(1));	

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization    
  	x_diff(3)=PhiNorm(x_diff(3));	
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1)=PhiNorm(z_diff(1));	
  //update state mean and covariance matrix
  x_(3)=PhiNorm(x_(3));
  x_ = x_ + K * z_diff;
  x_(3)=PhiNorm(x_(3));
  P_ = P_ - K*S*K.transpose();
#if ket_debug
  //print result
  cout << "UKF::UpdateState, updated state x_: " << endl << x_ << endl;
  cout << "UKF::UpdateState, updated state covariance P_: " << endl << P_ << endl;
#endif

  // Limit the counter so we do not exceed the vector size for 
	//multiple runs
  if (NISradar_counter_ < 220) {
	  NISradar_(NISradar_counter_++)=z_diff.transpose()*S.inverse()*z_diff;
  }
#if ket_debug  
  if ( NISradar_(NISradar_counter_-1) > 50) {
  cout<< S << endl;
  }
#endif
	return;
	
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
  VectorXd z=VectorXd(2);
  z(0)=meas_package.raw_measurements_[0];
  z(1)=meas_package.raw_measurements_[1];
  x_(3)=PhiNorm(x_(3));

	VectorXd z_pred = Hlidar_ * x_;
	
	VectorXd y = z - z_pred;
	
	MatrixXd Ht = Hlidar_.transpose();
	MatrixXd S = Hlidar_ * P_ * Ht + Rlidar_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_(3)=PhiNorm(x_(3));

	x_ = x_ + (K * y);
	
	x_(3)=PhiNorm(x_(3));

	int x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hlidar_) * P_;
	// INS calculation
	// Limit the counter so we do not exceed the vector size for 
	//multiple runs
	if (NISlidar_counter_ < 220) {
		NISlidar_(NISlidar_counter_++)=y.transpose()*Si*y;
	}
}

double UKF::PhiNorm(double phi){

	//while (phi> M_PI) phi-=2.*M_PI;
  //while (phi<-M_PI) phi+=2.*M_PI;
  
  //suggested by reviewer for code efficiency
  phi = atan2( sin(phi), cos(phi) );
  return phi;
}
  




