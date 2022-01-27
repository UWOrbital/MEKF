/** @file main.c
*	@brief Multiplicative Extended Kalman Filter (MEKF) source file
*	@date 4-Jan-2022
*	@version 01.00.00
* 
*	This file contains:
*	- Necessary function and calculations for a MEKF
*	.
*	which are relavant for Attitude Determination and Control Systems
*/

/** The University of Waterloo Orbital design team (UW Orbital)
*	The Attitude Determination and Control Systems sub-team (ADCS)
*	At the University of Waterloo
*	Written by: Hugh Dawson
*/

/**	Reference:
*	Name:		main_mekf
*	Author:		Rishav
*	Publisher:	risherlock (GitHub)
*	Date:		2021
*	Language:	C++
*	Source:		https://github.com/risherlock/Kalman/blob/master/mekf_murrell/cpp/main_mekf.cpp
*/

// Required libraries
#include <iostream>

// Required files
#include "Kalman.h"
#include "meschach-1.2/matrix.h"												// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "meschach-1.2/matrix2.h"												// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "quaternion/Quaternion.h"												// Courtesy of Martin Weigel					- Copyright (C) 2019


int main()
{
	// Constants
	float DELTA_TIME = 0.134;											// Timestep
		// Initial Values
	Quaternion initQ = { 0, 0, 0, 1 };									// Initial quaternion
	VEC initB = { 1.324, 2.4532, 3.4532 };								// Initial bias
	MAT initP;															// Initial state noise covariance
	initP = { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 1 } };
		// Noise standard deviation
	float sigmaR[2] = { 1.345, 2.654 };									// Vector measurements noise
	float sigmaG[3] = { 4.453, 2.534, 3.234 };							// Gyro noise
	float sigmaB[3] = { 2.435, 1.2345, 4.55 };							// Gyro-bias noise
		// MEKF Variables
	MAT Q;																// Process noise covariance
	MAT pHat;															// Updated estimate covariance
	Q = { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 1 } };
	pHat = { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 1 } };
	VEC dxHat = { 0, 0, 0, 0, 0, 0 };									// Updated state estimate
	VEC bHat = { 0, 0, 0 };												// Updated estiated covariance
	Quaternion qHat = { 0, 0, 0, 0 };
	float varR[2] = { 0, 0 };											// Vector sensors' noise variance
	

	// Variables
	float magnetometer[3] = { 0, 0, 0 };								// Magnetometer readings
	float accelerometer[3] = { 0, 0, 0 };								// Accelerometer readings
	float gyroscope[3] = { 0, 0, 0 };									// Gyroscope readings
	float magVector[3] = { 0, 0, 0 };									// Magnetometer vector in inertial frame
	float accVector[3] = { 0, 0, 0 };									// Accelerometer vector in inertial frame

	// Example Values
	magnetometer = { 1.32, 2.65, 3.63 };
	accelerometer = { 4.6534, 5.634, 6.634 };
	gyroscope = { 0.223, 2.34, 3.45 };
	magVector = { 7.53, 8.23, 9.97 };
	accVector = { 1.42, 2.432, 3.546 };

	MEKF(initQ, initB, initP, sigmaR, sigmaG, sigmaB, DELTA_TIME, bHat, qHat, pHat, varR, Q);	
	Estimate(gyroscope - bHat, magnetometer, accVector, DELTA_TIME, bHat, qHat, pHat, varR, dxHat);

	return bHat;
}
