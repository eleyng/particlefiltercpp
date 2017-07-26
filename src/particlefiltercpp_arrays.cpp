// Implementation of a PF on CPU in C/C++
// Using PF to estimate the position (state) on a 1D mass-spring-damper system

//#include "stdafx.h"
#include <functional>
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
#define Np 2000 // number of particles in each simulation


// Global variables
double dt = total_time / tstep;

// GPU FUNCTIONS
int dynamics(double* arr1, double* arr2);

int main(void) {

	ofstream out_data("plottingdata.dat");

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

	// Initialize truth, estimate, and noise
	x[0][0] = 0.0; // initial x = 0 m
	x[0][1] = 1.0; // initial dx = 0 m/s
	x_hat[0][0] = 1.0;
	x_hat[0][1] = 0.0;
	y[0] = x[0][0]; // initial measurement
	double sigma_w[n] = { 0.01, 0.01 };
	// Q is diagonal matrix of sigma_w; currently hardcoded
	double Q[n][n] = { { pow(sigma_w[0], 2.0), 0 },{ 0, pow(sigma_w[1],2.0) } };
	double sigma_v[nm] = { 0.2 };
	double R[nm] = { pow(sigma_v[0], 2.0) }; //R is diagonal matrix of sigma_v; currently hardcoded
	double w[n];
	double v[nm];
	double *v_ptr = v;
	double sigma_x0[n] = { 0.5, 0.5 };
	// distribution generator
	std::normal_distribution<double> distribution_w1(0.0, Q[0][0]);
	std::normal_distribution<double> distribution_w2(0.0, Q[1][1]);
	std::normal_distribution<double> distribution_v(0.0, R[0]);

	// Initialize particles
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	// mean of 0, std_dev of 1
	std::normal_distribution<double> distribution_x(0.0, 1.0);
	for (int p = 0; p < Np; p++) {
		for (int d = 0; d < n; d++) {
			xp[0][p][d] = distribution_x(generator)*sigma_x0[0] + x_hat[0][0];
			xp[1][p][d] = distribution_x(generator)*sigma_x0[1] + x_hat[0][1];
		}
	}

	// Create random number generators for w1, w2, and v
	auto dice_w1 = std::bind(distribution_w1, generator);
	auto dice_w2 = std::bind(distribution_w2, generator);
	auto dice_v = std::bind(distribution_v, generator);
	auto dice_x = std::bind(distribution_x, generator);

	int i = 0; int j = 0;
	//TIME STEP LOOP
	while (i < (tstep - 1)) {

		// UPDATE COVARIANCE (random gen)
		//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		//std::default_random_engine generator(seed);


		// TRUTH: update dynamics
		double *x_prev = x[i];
		double *x_curr = x[i + 1];
		dynamics(x_prev, x_curr);
		x[i + 1][0] += dice_w1();
		x[i + 1][1] += dice_w2();

		// MEASUREMENT: (heavily simplified) which is dynamics_with_noise
		y[i + 1] = x[i + 1][0] + dice_v();

		// randomize noise variables again

		// PARTICLES LOOP: propogate dynamics and measurement for each
		double yp[Np][nm];
		// for estimation and resampling
		double residual[Np][nm], q_i[Np];
		double est_sum[n] = { 0, 0 };
		double sumq = 0;

		for (int p = 0; p < (Np); p++) {
			// Particle state estimate: propogate dynamics
			dynamics(xp[i][p], xp[i + 1][p]);
			for (int j = 0; j < n; j++) {
				xp[i + 1][p][j] += dice_x();
				// calculate resiuduals
				est_sum[j] += xp[i + 1][p][j];
			}
			// Particle measurement
			yp[p][0] = xp[i + 1][p][0] + dice_v();
			residual[p][0] = y[i + 1] - yp[p][0]; //hardcoded
			q_i[p] = exp(-(residual[p][0] * residual[p][0] * (1 / R[0])) / 2);
			sumq += q_i[p];
		}
		x_hat[i + 1][0] = est_sum[0] / Np;
		x_hat[i + 1][1] = est_sum[1] / Np;
		// obtaining the weights and cumulative sum
		double c[Np];
		for (int p = 0; p < Np; p++) {
			q_i[p] /= sumq;
			if (p != 0) {
				c[p] = c[p - 1] + q_i[p];
			}
			else {
				c[0] = q_i[0];
			}
		}

		// PARTICLE estimate calculation
		double tmp_xp[Np][n];
		for (int p = 0; p < Np; p++) {
			double r = ((double)rand() / (RAND_MAX));
			//cout << "r " << r << endl;
			int ind = 0;
			while (c[ind] < r && ind < (Np - 1)) {
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
	
	
	
	
	
	// CPU CODE: FREE MEMORY
	// Deallocate memory: x, xp, y
	// De-Allocate memory to prevent memory leak
	for (int i = 0; i < tstep; ++i) {
		for (int j = 0; j < Np; ++j) {
			delete[] xp[i][j];
		}
		//cout << "x: " << x[i][0] << " x_hat: " << x_hat[i][0] << endl;
		out_data << i << " ";
		out_data << x[i][0] << " ";
		out_data << x_hat[i][0] << " ";
		out_data << x[i][1] << " ";
		out_data << x_hat[i][1] << "\n";
		delete[] xp[i];
		delete[] x[i];
		delete[] x_hat[i];
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
