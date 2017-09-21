#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

int main()
{
  double std_a_ = 30;
  double std_yawdd_ = 30;
	const int n_x_ = 5;
	const int n_aug_ = 7;
	double lambda_ = 3 - n_aug_;

	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  VectorXd x_ = VectorXd(n_x_);
  x_ <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;
 
  MatrixXd P_ = MatrixXd(n_x_,n_x_);
  P_ <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

 	VectorXd x_aug = VectorXd(n_aug_);
 	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	x_aug.head(n_x_) = x_;
  x_aug(5) = 0; // longitudinal acceleration noise has mean 0
  x_aug(6) = 0; //          yaw acceleration noise has mean 0

	P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

	MatrixXd A_aug = P_aug.llt().matrixL();

	int start = clock();


	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; ++i) {
		Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * A_aug.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A_aug.col(i);
	}
	
  // MatrixXd sig1 = (sqrt(lambda_ + n_aug_) * A_aug).colwise() + x_aug;
  // MatrixXd sig2 = (sqrt(lambda_ + n_aug_) * A_aug * -1).colwise() + x_aug;
  // Xsig_aug.col(0) = x_aug;
  // Xsig_aug.block(0, 1, n_aug_, n_aug_) = sig1;
  // Xsig_aug.block(0, 1+n_aug_, n_aug_, n_aug_) = sig2;

  int end = clock();
  std::cout << "it took " << end - start << "ticks, or " << ((float)end - start)/CLOCKS_PER_SEC << "seconds." << std::endl;

  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;  
}