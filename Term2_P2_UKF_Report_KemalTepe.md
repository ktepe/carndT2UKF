# **Unscented Kalman Filter Project Starter Code**
Self-Driving Car Engineer Nanodegree Program

## Kemal Tepe, ketepe@gmail.com

### Objective: Objective of this project to implement Unscented Kalman Filter (UKF) for sensor fusion application which combines lidar and radar sensor information.

### Summary: 

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. Original project description is provided by [README.md](./README.md). In this file, my explorations would be summarized.

### The goals / steps of this project are the following:

* Implement UKF for Radar and Lidar sensor data fusion.

* Explore parameter dependensies.

* Discuss results and future enhancements.

### 1. Implement UKF for Radar and Lidar sensor data fusion.

The implementation utilizes UKF lectures provided by Udacity. Once UKF lecture is completed, the rest is the integration of these modules (functions) in the ukf.cc and ukf.h files.

Prediction phase can be accomplished by following function.

```c++
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
}
```

We need two update implementation one for Lidar, and one for Radar. Lidar would be acomplished by simply using Linear Kalman Filter (LKF) using the following code:

```c++

void UKF::UpdateLidar(MeasurementPackage meas_package) {
 
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
  NISlidar_(NISlidar_counter_++)=y.transpose()*Si*y;
	
}

```

Update function  for radar needs to be implemented using UKF, and it is done by the following  function:

```c++

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
  for (int i=1; i < 2 * n_aug+1; i++) {  
    //2n+1 weights
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

  //VectorXd z_diff = z-z_pred;
  
  NISradar_(NISradar_counter_++)=z_diff.transpose() * S.inverse() * z_diff;

	return;
	
}
```
### 2. Explore parameter dependensies.

Performance of the UKF was measured using RMSE values of x and y positions as well as velocities in x and y directions. Another parameter which is important to observe Normalized Innovation Squared (NIS) values for Lidar and Radar. These two are provided below.

![alt text](./Doc/NISLidar.jpg) *NIS of Lidar, the values are consistenly below 5.9.*
![alt text](./Doc/NISRadar.jpg) *NIS of Radar, the values are consistenly below 7.8.* (P_diagonal=[0.5 0.5 10 10 0.5].*

Process covariance values also plays an important role in the performance of UKF. I changed these to reduce the error as well as get consistent NIS values.

*Table: for standard variances for acceleration and yaw change during the  different parameters for Radar only UKF*


| Parameters | | | | | | | 
| :------------- |-------------:| -----:|-----:|------:|------:|------:|
| Std_a         | 0.5 	| 1.0 	| 1.5 	| 2.0 	| 1.0 | 1.0 | 
| Std_yawdd     | 1.0 	| 1.0 	| 1.0 	| 1.0 	| 0.5 | 1.5 |
| RMSE x |  0.1568    |   0.1521 | 0.1570 | 0.1629 | 0.1511 | 0.1532 |
| RMSE y |  0.2106     |   0.2160 | 0.2262 | 0.2353 |  0.2057 | 0.2194 |
| RMSE Vx |  0.3673   |   0.3575 |  0.3636 | 0.3761 | 0.3606 | 0.3651 |
| RMSE Vy |  0.3297     |    0.3150 | 0.3264 | 0.3364 | 0.3158 | 0.3304 |

Std_a=1.0, and Std_yawdd=1.0 seem to be a good comprimise among different valuesto use from this table


### 3. Results and Future Enhancements

*. Larger Std_a values yield larger errors in x direction, and larger Std_yawdd values generate larger errors in y direction.

*. Better initialization covariance matrix (P) provides a better error estimates initially. However, the process goes to steady state values even initial P is not that well optimized.

As future enhancement these can be further optimized by investigating large number of samples.

Sensor fusion significantly reduce the estimation errors as evidenced by the following table, where simulation was run only with Lidar, only with Radar and with both.

|RMSE Metric| Only Lidar | Only Radar | Both Lidar and Ridar |
|----|----|-----|-----|
|RMSE x | 0.1022 |0.2175 | 0.0642 |
|RMSE y | 0.0981 | 0.2577 | 0.0827 |
|RMSE Vx | 0.6090 |0.4143 | 0.3262 |
|RMSE Vy | 0.2573 |0.3608 | 0.1865 |


