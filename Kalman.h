/** @file Kalman.h
*	@brief Multiplicative Extended Kalman Filter (MEKF) source file
*	@date 14-Jan-2022
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
*	Name:		mekf
*	Author:		Rishav
*	Publisher:	risherlock (GitHub)
*	Date:		2021
*	Language:	C++
*	Source:		https://github.com/risherlock/Attitude-Estimation/blob/master/mekf_murrell/cpp/mekf.h
*/

#pragma once

// Required libraries
#include <iostream>

// Required files
#include "Kalman.h"
#include "meschach-1.2/matrix.h"					// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "meschach-1.2/matrix2.h"					// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "quaternion/Quaternion.h"					// Courtesy of Martin Weigel					- Copyright (C) 2019

void ComputeProcessCovariance(const float sigmaG[3], const float sigmaB[3], float dt, Quaternion qHat, MAT& Q);

void PropQuaternion(const Quaternion& q, const VEC& w, const float dt);

void UpQuaternion(const Quaternion& q, const VEC& delAlpha, Quaternion&	qHat);

void PropStateCovariance(const MAT& p, const MAT& q, const VEC& w, float dt, MAT& pHat);

void Estimate(const VEC w, const VEC mb, const VEC mr, float dt);

void MEKF(const Quaternion qInit, const VEC& bInit, const MAT& pInit, const float sigmaR[2], 
		  const float sigmaG[3], const float sigmaB[3], float dt, VEC& bHat, Quaternion& qHat, MAT& pHat, 
		  float& varR, MAT& Q)
{
	qHat = qInit;
	bHat = bInit;
	pHat = pInit;

	ComputeProcessCovariance(sigmaG, sigmaB, dt, qHat, Q);

	for (int i = 0; i < 2; i++)
	{
		varR[i] = sigmaR[i] * sigmaR[i];
	}
}

/*
	Inputs:
		mb[n] = n sensor measuremnt vectors in body fram [3xn]
		mr[n] = n measurement vectors in inertial frame [3xn]
		dt	  = Sample time (sec)

	Outputs:
		Updates quaternion and bias
*/

void Estimate(const VEC w, const VEC mb, const VEC mr, float dt, VEC& bHat, Quaternion& qHat, MAT& pHat,
			  float& varR, VEC& dxHat)
{
	// Propagate state, gyro-bias, quaternion, and estimate covariance
	// Note: Gyro bias propagation is constant
	PropQuaternion(qHat, w, dt);
	PropStateCovariance(pHat, Q, w, dt, pHat);
	// TO-DO: Find/Create an appropriate replacement for this function
	// MAT A = qHat.getDCM;
	MAT A;
	A = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	// TO-DO: Setup input sensors
	const VEC mb = { 0, 0, 0 };
	const VEC mr = { 0, 0, 0 };

	// Murrell's loop
	for (int i = 0; i < 2; i++)
	{
		MAT I3;
		MAT I6;
		I3 = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
		I6 = { { 1, 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 1 } };

		// Sensitivity matrix (H[3x6]) and residual (e[3x1])
		VEC Ar = A * mr[i];
		VEC e = mb[i] - Ar;
		MAT H;
		H = { { 0, -Ar[2], Ar[1] }, { Ar[2], 0, -Ar[0] }, { -Ar[1], Ar[0], 0 } };

		// MEKF update
		MAT K;
		K = pHat * m_transp(H) * m_inverse(H * pHAt * m_transp(H) + varR[i] * I3);
		dxHat = dxHat + K * (e - H * dxHat);
		pHat = (I6 - K * H) * pHat;		
	}

	// Update qHat and bHat
	VEC dQ = { dxHat[0], dxHat[1], dxHat[2] };
	VEC dB = { dxHat[3], dxHat[4], dxHat[5] };
	UpQuaternion(qHat, dQ, qHat);
	bHat = bHat + dB;
}

// Discrete Q using gyro and gyro-bias noise standard deviations
void ComputeProcessCovariance(const float sigmaG[3], const float sigmaB[3], float dt, Quaternion qHat, MAT& Q)
{
	float sqGx = sigmaG[0] * sigmaG[0];
	float sqGy = sigmaG[1] * sigmaG[1];
	float sqGz = sigmaG[2] * sigmaG[2];
	float sqBx = sigmaB[0] * sigmaB[0];
	float sqBy = sigmaB[1] * sigmaB[1];
	float sqBz = sigmaB[2] * sigmaB[2];

	float a = dt * dt * dt / 3.0;
	float b = 0.5 * dt * dt;

	Q = { { sqGx * dt + sqBx * a, 0, 0,		sqBx * b, 0, 0 },
		  { 0, sqGy * dt + sqBy * a, 0,		0, sqBy * b, 0 },
		  { 0, 0, sqGz * dt + sqBz * a,		0, 0, sqBz * b },
		  { sqBx * b, 0, 0,					sqBx * dt, 0, 0 },
		  { 0, sqBy * b, 0,					0, sqBy * dt, 0 },
		  { 0, 0, sqBz * b,					0, 0, sqBz * dt } };

	MAT Y;
	Y = { { -1,  0,  0,  0,  0,  0 },
		  {  0, -1,  0,  0,  0,  0 },
		  {  0,  0, -1,  0,  0,  0 },
		  {  0,  0,  0,  1,  0,  0 },
		  {  0,  0,  0,  0,  1,  0 },
		  {  0,  0,  0,  0,  0,  1 } };

	Q = Y * Q * bd_transp(Y)										// YQY'
}

void PropQuaternion(const Quaternion& q, const VEC& w, const float dt)
{
	float omegaTol = 1e-5;
	float n = v_norm1(w);
	if (n > omegaTol)
	{
		MAT Omega;
		float c = cos(0.5 * n * dt);
		float s = sin(0.5 * n * dt) / n;
		float x = w[0] * s;
		float y = w[1] * s;
		float z = w[2] * s;
		Omega = { { c, z, -y, x }, { -z, c, x, y }, { y, -x, c, z }, { -x, -y, -z, c } };
		qHat = Omega * qHat;
	}
}

void UpQuaternion(const Quaternion& q, const VEC& delAlpha, Quaternion& qHat)
{
	MAT Xi;
	Xi = { {  q[3], -q[2],  q[1] }, 
		   {  q[2],  q[3], -q[0] }, 
		   { -q[1],  q[0],  q[3] }, 
		   { -q[0], -q[1], -q[2] } };

	qHat = q + 0.5 * Xi * del_alpha;
	qHat = qHat / v_norm1(qHat);
}

void PropStateCovariance(const MAT& p, const MAT& q, const VEC& w, float dt, MAT& pHat)
{
	float omegaTol = 1e-5;
	float n = v_norm1(w);
	MAT Phi;

	if (n > omegaTol)
	{
		float s = sin(n * dt);
		float c = cos(n * dt);
		float p = s / n;
		float q = (1 - c) / (n * n);
		float r = (n * dt - s) / (n * n * n);

		MAT wX;
		MAT I;
		MAT wXwX;
		MAT Phi00;
		MAT Phi01;
		wX = { { 0, -w[2], w[1] }, { w[2], 0, -w[0] }, { -w[1], w[0], 0 } };
		I = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
		wXwX = wX * wX;
		Phi00 = I - wX * p + wXwX * q;
		Phi01 = wX * q - I * dt - wXwX * r;

		Phi = { { Phi00[0,0],	Phi00[0,1],		Phi00[0,2],		Phi01[0,0],		Phi01[0,1],		Phi01[0,2]	},
				{ Phi00[1,0],	Phi00[1,1],		Phi00[1,2],		Phi01[1,0],		Phi01[1,1],		Phi01[1,2]	},
				{ Phi00[2,0],	Phi00[2,1],		Phi00[2,2],		Phi01[2,0],		Phi01[2,1],		Phi01[2,2]	},
				{ 0,			0,				0,				1,				0,				0			},
				{ 0,			0,				0,				0,				1,				0			},
				{ 0,			0,				0,				0,				0,				1			} };
	}
	else
	{
		Phi = { {   1,   0,   0, -dt,   0,   0 },
				{   0,   1,   0,   0, -dt,   0 },
				{   0,   0,   1,   0,   0, -dt },
				{   0,   0,   0,   1,   0,   0 },
				{   0,   0,   0,   0,   1,   0 },
				{   0,   0,   0,   0,   0,   1 } };
	}

	pHat = Phi * p * m_transp(Phi) + q;
}
