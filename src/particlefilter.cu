// Implementation of a PF on CPU in C/C++
// Using PF to estimate the position (state) on a 1D mass-spring-damper system

//#include "stdafx.h"
#include <functional>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda.h>
using namespace std;

#define n 2 // number of states
#define nm 1 // number of measurements
#define b 0.5 // Ns/m
#define k 3.0 //  N/m
#define m 10.0 // kg
#define tstep 1000 // number of steps for the entire simulation
#define total_time 30.0 // seconds
#define Np 1024*1 // number of particles in each simulation = num_SM * threads
#define THREADS 1024


// Global variables
float dt = total_time / tstep;
int BLOCKS = Np / THREADS + (Np % THREADS);

// Device constant memory
__constant__ float sigma_x0[n] = { 0.5, 0.5 };

// Device variables
__device__ float yp[Np][nm];
__device__ float residual[Np][nm];
__device__ float q_i[Np];
__device__ float est_sum[n] = {0.0, 0.0};
__device__ float sumq = 0.0;
__device__ float xp_tmp[2][Np][n];



// Device function: propogate particles
//__global__ void propogate (float ***xp, float **x_hat) {
	// call device function to propogate dynamics
	

// Device function: Initialize particles
__global__ void initializeParticles (float ***xp, float **x_hat) {

	int pid = threadIdx.x + blockIdx.x * blockDim.x;
	curandState state;
	curand_init(1234, pid, 0, &state); //seed, squence, offset, 
	for (int j=0; j<n; j++){
		xp[0][pid][j] = sigma_x0[j] + x_hat[0][j];
	}
}

// Device function: project particles for this particular time step
__global__ void proj_particles ( float ***xp_tmp, float *x_hat) { //passing in x_hat[i+1]; xp_tmp is xp_curr and xp_prev
	
	float **xp_curr = xp_tmp[1];
	float **xp_prev = xp_tmp[0];
	dynamics(xp_prev, xp_curr);
	
	int pid = threadIdx.x + blockIdx.x * blockDim.x;
	for (int j = 0; j < n; j++) {
		xp_curr[pid][j] += dice_x();
		// calculate resiuduals
		est_sum[j] += xp_curr[pid][j];
	}
	// Particle measurement
	yp[pid][0] = xp_curr[pid][0] + dice_v();
	residual[pid][0] = y[i + 1] - yp[pid][0]; //hardcoded
	q_i[pid] = exp(-(residual[pid][0] * residual[pid][0] * (1 / R[0])) / 2);
	sumq += q_i[pid];
	x_hat_curr[0] = est_sum[0] / Np;
	x_hat_curr[1] = est_sum[1] / Np;
}

// Host function
__host__ __device__ int dynamics(float* arr1, float* arr2);

	float x_prev = arr1[0];
	float dx_prev = arr1[1];

	// acceleration from previous state to find velocity
	float ddx = -(b / m)* (dx_prev)-(k / m)* (x_prev);


	// velocity (index 1) and position (index 2) calcs
	arr2[1] = dx_prev + (ddx * dt);
	arr2[0] = x_prev + (dx_prev * dt);
	//std::cout << "print " << arr1;

	return 0;
}



// Host MAIN
int main(void) {

	ofstream out_data("plottingdatagpu.dat");

	// Truth Equation
	float **x;
	//float x[tstep][n]; //truth x
	float *y;
	//float y[tstep]; //measure x
	float ***xp;
	//float xp[tstep][Np][n]; //particle x
	float **x_hat;
	// float x_hat[tstep][n]; //x estimate from particles

	// Allocate memory
	x = new float*[tstep];
	y = new float[tstep];
	xp = new float**[tstep];
	x_hat = new float*[tstep];
	for (int i = 0; i < tstep; ++i) {
		x[i] = new float[n];
		xp[i] = new float*[Np];
		x_hat[i] = new float[n];
		for (int p = 0; p < Np; ++p) {
			xp[i][p] = new float[n];
		}
	}

	// Initialize truth, estimate, and noise
	x[0][0] = 0.0; // initial x = 0 m
	x[0][1] = 1.0; // initial dx = 0 m/s
	x_hat[0][0] = 1.0;
	x_hat[0][1] = 0.0;
	y[0] = x[0][0]; // initial measurement
	float sigma_w[n] = { 0.01, 0.01 };
	// Q is diagonal matrix of sigma_w; currently hardcoded
	float Q[n][n] = { { powf(sigma_w[0], 2.0), 0 },{ 0, powf(sigma_w[1],2.0) } };
	float sigma_v[nm] = { 0.2 };
	float R[nm] = { powf(sigma_v[0], 2.0) }; //R is diagonal matrix of sigma_v; currently hardcoded
	//float w[n];
	//float v[nm];
	//float *v_ptr = v;
	//float sigma_x0[n] = { 0.5, 0.5 };
	// distribution generator
	std::normal_distribution<float> distribution_w1(0.0, Q[0][0]);
	std::normal_distribution<float> distribution_w2(0.0, Q[1][1]);
	std::normal_distribution<float> distribution_v(0.0, R[0]);
	std::normal_distribution<float> distribution_x(0.0, 1.0);
	// Create random number generators for w1, w2, and v
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	auto dice_w1 = std::bind(distribution_w1, generator);
	auto dice_w2 = std::bind(distribution_w2, generator);
	auto dice_v = std::bind(distribution_v, generator);
	auto dice_x = std::bind(distribution_x, generator);

	//CPU: TRUTH AND MEASUREMENT WITH NOISE
	int i = 0; 
	while (i < (tstep - 1)) {

		// TRUTH: update dynamics
		float *x_prev = x[i];
		float *x_curr = x[i + 1];
		dynamics(x_prev, x_curr);
		x[i + 1][0] += dice_w1();
		x[i + 1][1] += dice_w2();

		// MEASUREMENT: (heavily simplified) which is dynamics_with_noise
		y[i + 1] = x[i + 1][0] + dice_v();
		i++;
	}
	float ***d_xp; float **d_x_hat;
	cudaMalloc((void**)&d_xp, tstep * Np * n * sizeof(float));
	cudaMalloc((void**)&d_x_hat, tstep * n *sizeof(float));
	cudaMemcpy(d_xp, xp, tstep * Np * n * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x_hat, x_hat, tstep * n * sizeof(float), cudaMemcpyHostToDevice);
	 
	// GPU: Initialize particles
	initializeParticles <<<BLOCKS, THREADS>>> (d_xp, d_x_hat);
	
	/*for (int p = 0; p < Np; p++) {
		for (int d = 0; d < n; d++) {
			xp[0][p][d] = distribution_x(generator)*sigma_x0[0] + x_hat[0][0];
			xp[1][p][d] = distribution_x(generator)*sigma_x0[1] + x_hat[0][1];
		}
	}*/
	cudaMemcpy(xp , d_xp, tstep * Np * n * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(x_hat, d_x_hat, tstep * n * sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(d_xp); cudaFree(d_x_hat);
	
	// GPU: Particles loop
	i = 0;
	while (i < (tstep - 1)) {
		/* yp, residual, q_i, est_sum, and sumq go on device only 
		// PARTICLES LOOP: propogate dynamics and measurement for each
		float yp[Np][nm];
		// for estimation and resampling
		float residual[Np][nm], q_i[Np];
		float est_sum[n] = { 0, 0 };
		float sumq = 0; */
		
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
		// obtaining the weights and cumulative sum - faster on CPU?
		float c[Np];
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
		float tmp_xp[Np][n];
		for (int p = 0; p < Np; p++) {
			float r = ((float)rand() / (RAND_MAX));
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


