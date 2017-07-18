// Implementation of a PF on CPU in C/C++
// Using PF to estimate the position (state) on a 1D mass-spring-damper system

#include "stdafx.h"
//#include "pf_func.h"
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
using namespace std;

#define n 2 // number of states
#define nm 1 // number of measurements
#define b 0.5 // Ns/m
#define k 3.0 //  N/m
#define m 10.0 // kg
#define tstep 1000 // number of steps for the entire simulation
#define total_time 30.0 // seconds
#define Np 1000 // number of particles in each simulation


// Global variables
double dt = total_time / tstep;

int dynamics(double* arr1, double* arr2);
//int particle_generator(double ***xp);
// int covariance(double w[n], double *v);
// double estimate(double particles[Np][n]);

int main(void) {

	ofstream out_data("plottingdata.txt");

	// Truth Equation
	double **x;
	//double x[tstep][n]; //truth x
	double *y;
	//double y[tstep]; //measure x
	double ***xp;
	//double xp[tstep][Np][n]; //particle x
	double **x_hat;
	// double x_hat[tstep][n]; //x estimate from particles

	// Allocate memory
	x = new double*[tstep];
	y = new double[tstep];
	xp = new double**[tstep];
	x_hat = new double*[tstep];
	for (int i = 0; i < tstep; ++i) {
		x[i] = new double[n];
		xp[i] = new double*[Np];
		x_hat[i] = new double[n];
		for (int p = 0; p < Np; ++p) {
			xp[i][p] = new double[n];
		}
	}
	
	// Initialize truth and noise
	x[0][0] = 0.0; // initial x = 0 m
	x[0][1] = 1.0; // initial dx = 0 m/s
	x_hat[0][0] = 0.0;
	x_hat[0][1] = 1.0;
	y[0] = x[0][0]; // initial measurement
	double Q[n][n] = { { 0.1, 0 },{ 0, 0.05 } };
	double R[nm][nm] = { 0.1 };
	double w[n];
	double v[nm];
	double *v_ptr = v;

	// distribution generator
	std::normal_distribution<double> distribution_w1(0.0, Q[0][0]);
	std::normal_distribution<double> distribution_w2(0.0, Q[1][1]);
	std::normal_distribution<double> distribution_v(0.0, R[0][0]);

	// iterate each time step
	int i = 0; int j = 0;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	// mean of 0, std_dev of 1
	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int p = 0; p < Np; p++) {
		for (int d = 0; d < n; d++) {
			xp[0][p][d] = distribution(generator);
			xp[1][p][d] = distribution(generator);
		}
		//cout << "xp initial: " << xp[0][p][0] << endl;
	}

	//TIME STEP LOOP
	while (i < (tstep - 1)) {

		// grab previous and current state
		double *x_prev = x[i];
		double *x_curr = x[i + 1];


		// update covariance (random gen)
		w[0] = distribution_w1(generator);
		w[1] = distribution_w2(generator);
		v[0] = distribution_v(generator);
		// TRUTH: run dynamics
		dynamics(x_prev, x_curr);

		for (j = 0; j < n; j++) {
			// TRUTH: add noise
			x[i + 1][j + 1] += (w[j]/100);
		}
		// MEASURE: which is dynamics_with_noise
		y[i + 1] = x[i + 1][0] + v[0];

		// PARTICLES LOOP: propogate dynamics and measurement for each
		double yp[Np][nm];
		for (int p = 0; p < (Np); p++) {
			// Particle state estimate: propogate dynamics
			dynamics(xp[i][p], xp[i + 1][p]);
			//cout << "i: " << i << "\np: " << p << "\nposition" << xp[i][p][0] << " velocity "  << xp[i][p][1] << endl;
			// Particle measurement
			yp[p][0] = xp[i + 1][p][0];
		}
		
		// ESTIMATE AND RESAMPLE
		double residual[Np][nm], q_i[Np];
		double est_sum[n] = { 0, 0 };
		double sumq = 0;
		// calculating the residuals for each particle
		for (int p = 0; p < Np; p++) {
			for (int j = 0; j < n; j++) {
				est_sum[j] += xp[i + 1][p][j];
			}
			residual[p][0] = y[i + 1] - yp[p][0];
			q_i[p] = exp(-(residual[p][0] * residual[p][0] * (1 / R[0][0])) / 2);
			//cout << "q: " << q_i[p] << endl;
			sumq += q_i[p];
		}
		x_hat[i + 1][0] = est_sum[0] / Np;
		x_hat[i + 1][1] = est_sum[1] / Np;
		//cout << "x_hat " << x_hat[i + 1][0] <<endl;

		// dividing by the q_total
		for (int p = 0; p < Np; p++) {
			q_i[p] /= sumq;
			//cout << i << " q " << q_i[p] << endl;
		}
		// calculating c vector
		double c[Np];
		c[0] = q_i[0];
		for (int p = 1; p < Np; p++) {
			c[p] = c[p - 1] + q_i[p];
			//cout << p << " c " << c[p] << endl;
		}

		double tmp_xp[Np][n];

		for (int p = 0; p < Np; p++) {
			double r = ((double)rand() / (RAND_MAX));
			//cout << "r " << r << endl;
			int ind = 0;
			while (c[ind] < r && ind < (Np-1)) {
				ind++;
			}
			//cout << "ind " << ind << endl;
			tmp_xp[p][0] = xp[i + 1][ind][0];
			tmp_xp[p][1] = xp[i + 1][ind][1];
		}
		for (int p = 0; p < Np; p++) {
			xp[i + 1][p][0] = tmp_xp[p][0];
			xp[i + 1][p][1] = tmp_xp[p][1];
		}
		i++;
	}

	// print results to file
	

	// Deallocate memory: x, xp, y
	// De-Allocate memory to prevent memory leak
	for (int i = 0; i < tstep; ++i) {
		for (int j = 0; j < n; ++j) {
			delete[] xp[i][j];
		}
		//cout << "x: " << x[i][0] << " x_hat: " << x_hat[i][0] << endl;
		out_data << i << " ";
		out_data << x[i][0] << " ";
		out_data << x_hat[i][0] << " ";
		out_data << x[i][1] << " ";
		out_data << x_hat[i][1] << "\n";
		delete[] xp[i];
		//delete[] x[i];
	}
	delete[] xp;
	delete[] x;
	delete[] y;
	delete[] x_hat;

	return 0;
}

// numerically calculate dynamics
int dynamics(double* arr1, double* arr2) {

	double x_prev = arr1[0];
	double dx_prev = arr1[1];

	// acceleration from previous state to find velocity
	double ddx = -(b / m)* (dx_prev)-(k / m)* (x_prev);


	// velocity (index 1) and position (index 2) calcs
	arr2[1] = dx_prev + (ddx * dt);
	double c = arr2[1];
	arr2[0] = x_prev + (dx_prev * dt);
	double d = arr2[0];
	//std::cout << "print " << arr1;

	return 0;
}

// initial particle generator
/* int particle_generator(double xp) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	// mean of 0, std_dev of 1
	std::normal_distribution<double> distribution(0.0, 1.0);

	for (int p = 0; p < Np; p++) {
		for (int d = 0; d < n; d++) {
			xp[0][p][d] = distribution(generator); 
			xp[1][p][d] = distribution(generator);
		}
	}
	return 0;
} */

/*int random_noise(double w[n], double *v, double *R, double *Q) {
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double cov_v = *R;
double cov_w = *Q;

std::normal_distribution<double> distribution_w(0.0, *Q[0][0]);
std::normal_distribution<double> distribution_w(0.0, c);
std::normal_distribution<double> distribution_v(0.0, cov_v);
for (int i = 0; i < n; i++) {
w[i] = distribution_w(generator);
}
*v = distribution_v(generator);

return 0;
}*/

/* double estimate(double particles[Np][n]) {
// calculate the state estimate for a particular time step
double sum = 0.0;
for (int p = 0; p < Np; p++) {
sum += particles[p][0];
}
double x_hat = (1 / Np) * sum;
return sum;
} */