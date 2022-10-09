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
 *	Written by: Nathan Tang
 */

/**	Reference:
 *	Name:		mekf
 *	Author:		Rishav
 *	Publisher:	risherlock (GitHub)
 *	Date:		2021
 *	Language:	C++
 *	Source:		https://github.com/risherlock/Attitude-Estimation/blob/master/mekf_murrell/cpp/mekf.h
 */

/** Notes:
 * Elements of vectors and matricies (->ve or ->me) are of type double. May need to type cast to float for use of elements in MAT *Q
 * Expressions for dxHat and pHat have complex combinations of adding and multiplying operations. Need to use intermediate steps
 * Addition of quaternions and 4d vectors with require manual expressions
 * Add a getDCM() function. Returns a 3x3 matrix using the qHat quaternion. Inputs quaternion and outputs a 3x3 matrix
 * Note order for matrix and vector mulitplication
 * Might have to convert all quaternions to quaternion pointers
 * mb and mr are 3d vectors of 2d arrays. Use a 3x2 matrix instead and extract the columns (function for it) as the vectors
 */

// Required libraries
#include <math.h>

// Required files
#include "matrix.h"				   // Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "matrix2.h"			   // Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "quaternion/Quaternion.h" // Courtesy of Martin Weigel					- Copyright (C) 2019

/* Struct for MEKF Members */
struct Kalman
{
	MAT *Q;		   // Process noise covariance
	VEC *dxHat;	   // Updated state estimate
	MAT *pHat;	   // Updated estimate covariance
	float varR[2]; // Vector sensors' noisevariance [1xn]
	Quaternion *qHat;
	VEC *bHat;
};

// Function prodotype for getDCM. Uses qHat in *k to modify DCM
void getDCM(struct Kalman *k, MAT *DCM);

void MEKFinit(struct Kalman *k, const Quaternion qInit, const VEC *bInit, const MAT *pInit, const float sigmaR[2],
			  const float sigmaG[3], const float sigmaB[3], float dt);
void ComputeProcessCovariance(struct Kalman *k, const float sigmaG[3], const float sigmaB[3], float dt);
void PropQuaternion(struct Kalman *k, const Quaternion *q, const VEC *w, const float dt);
void UpQuaternion(struct Kalman *k, const Quaternion *q, const VEC *delAlpha);
void PropStateCovariance(struct Kalman *k, const MAT *P, const MAT *YQY_, const VEC *w, float dt);
void Estimate(struct Kalman *k, const VEC *w, const MAT *mb, const MAT *mr, float dt);

// Constructor function for Kalman Struct
void MEKFinit(struct Kalman *k, const Quaternion *qInit, const VEC *bInit, const MAT *pInit, const float sigmaR[2],
			  const float sigmaG[3], const float sigmaB[3], float dt)
{
	k->Q = m_get(6, 6);
	k->qHat = (Quaternion *) qInit;
	v_copy(bInit, k->bHat);
	m_copy(pInit, k->pHat);
	ComputeProcessCovariance(k, sigmaG, sigmaB, dt);

	for (uint8_t i = 0; i < 2; i++)
	{
		k->varR[i] = sigmaR[i] * sigmaR[i];
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

void Estimate(struct Kalman *k, const VEC *w, const MAT *mb, const MAT *mr, const float dt)
{
	// Propagate state, gyro-bias, quaternion and estimate covariance
	// NoteL Gyro bias propagation is constant. Eq.(7.42b) (i.e. b=b)
	PropQuaternion(k, w, dt);					  // Eq.(7.39)
	PropStateCovariance(k, k->pHat, k->Q, w, dt); // Eq.(7.43)
	v_zero(k->dxHat);							  // initialize dxHat to a vector of zeros

	MAT *A;
	A = m_get(3, 3); // Declare Attitude Matrix
	getDCM(k, A);	 // Set Attitude Matrix

	// Murrell's loop
	for (uint8_t i = 0; i < 2; i++)
	{
		MAT *I3;
		I3 = m_get(3, 3);
		I3 = m_ident(MNULL); // Set a 3x3 identity matrix
		MAT *I6;
		I6 = m_get(6, 6);
		I6 = m_ident(MNULL); // Set a 3x3 identity matrix

		// Sensitivity matrix (H[3x6]) and residal (e[3x1])
		VEC *Ar;
		Ar = v_get(3); // Vector Ar = A * mr[i]
		Ar = vm_mlt(A, get_col(mr, i, VNULL), VNULL);

		MAT *H;
		H = m_get(3, 6);
		H = m_zero(MNULL);

		H->me[0][1] = -Ar->ve[2];
		H->me[0][2] = Ar->ve[1];
		H->me[1][0] = Ar->ve[2];
		H->me[1][2] = -Ar->ve[0];
		H->me[2][0] = -Ar->ve[1];
		H->me[2][1] = Ar->ve[0];

		VEC *e;
		e = v_get(3);
		e = v_sub(get_col(mb, i, VNULL), Ar, VNULL);

		// MEKF update
		MAT *K;
		K = m_get(6, 3);

		MAT *A1;
		A1 = m_get(3, 3);
		A1 = sm_mlt(k->varR[i], I3, MNULL); // holds varR[i] * I3
		MAT *A2;
		A2 = m_get(3, 3);
		A2 = m_mlt(H, k->pHat, MNULL); // holds H * pHat
		MAT *A3;
		A3 = m_get(3, 3);
		A3 = m_mlt(A2, m_transp(H, MNULL), MNULL); // holds H * pHat * trans(H)
		K = m_mlt(m_mlt(k->pHat, m_transp(H, MNULL), MNULL), m_inverse(m_add(A3, A1, MNULL), MNULL), MNULL);

		// Eq.(7.30)
		// dx_hat = dx_hat + K * (e - H * dx_hat);
		k->dxHat = v_add(k->dxHat, vm_mlt(K, v_sub(e, vm_mlt(H, k->dxHat, VNULL), VNULL), VNULL), VNULL);
		k->pHat = m_mlt(m_sub(I6, m_mlt(K, H, MNULL), MNULL), k->pHat, MNULL);
	}
	// Update q_hat and b_hat
	VEC *dq;
	dq = v_get(3);
	dq->ve[0] = k->dxHat->ve[0];
	dq->ve[1] = k->dxHat->ve[1];
	dq->ve[2] = k->dxHat->ve[2];

	VEC *db;
	db = v_get(3);
	db->ve[0] = k->dxHat->ve[3];
	db->ve[1] = k->dxHat->ve[4];
	db->ve[2] = k->dxHat->ve[5];

	UpQuaternion(k, k->qHat, dq); // Eq.(7.34)
	k->bHat = v_add(k->bHat, db, VNULL);  // Eq.(7.32)
}

// Discrete Q using gyro and gyro-bias noise standard deviations
void ComputeProcessCovariance(struct Kalman *k, const float sigmaG[3], const float sigmaB[3], float dt)
{
	// Eq.(7.46)
	float sq_gx = sigmaG[0] * sigmaG[0];
	float sq_gy = sigmaG[1] * sigmaG[1];
	float sq_gz = sigmaG[2] * sigmaG[2];
	float sq_bx = sigmaB[0] * sigmaB[0];
	float sq_by = sigmaB[1] * sigmaB[1];
	float sq_bz = sigmaB[2] * sigmaB[2];

	float a = (dt * dt * dt) / 3.0;
	float b = 0.5 * dt * dt;

	// First Quadrant
	k->Q->me[0][0] = sq_gx + dt + sq_bx * a;
	k->Q->me[1][1] = sq_gy + dt + sq_by * a;
	k->Q->me[2][2] = sq_gz + dt + sq_bz * a;

	// Second Quadrant
	k->Q->me[0][3] = sq_bx * b;
	k->Q->me[1][4] = sq_by * b;
	k->Q->me[2][5] = sq_bz * b;

	// Third Quadrant
	k->Q->me[3][0] = k->Q->me[0][3];
	k->Q->me[4][1] = k->Q->me[1][4];
	k->Q->me[5][2] = k->Q->me[2][5];

	// Fourth Quadrant
	k->Q->me[3][3] = sq_bx * dt;
	k->Q->me[4][4] = sq_by * dt;
	k->Q->me[5][5] = sq_bz * dt;

	MAT *Y;
	Y = m_get(6, 6);
	Y = m_zero(MNULL);
	Y->me[0][0] = -1;
	Y->me[1][1] = -1;
	Y->me[2][2] = -1;
	Y->me[3][3] = 1;
	Y->me[4][4] = 1;
	Y->me[5][5] = 1;

	// YQY' is used in Eq.(7.43)
	MAT *A1;
	A1 = m_get(6, 6); // holds trans(Y)
	MAT *A2;
	A2 = m_get(6, 6); // holds Y * Q
	A1 = m_transp(Y, MNULL);
	A2 = m_mlt(Y, k->Q, MNULL);
	k->Q = m_mlt(A2, A1, MNULL); // Not sure if the OUT paramater can be the same as one of the input paramaters
}

// Discrete-time quaternion propagation
void PropQuaternion(struct Kalman *k, const VEC *w, const float dt)
{
	float omega_tol = 1e-5;
	float n = (float)v_norm1(w); // v_norm1() returns a double -> cast to a float
	if (n > omega_tol)
	{
		// Eq.(7.40)
		MAT *Omega;
		Omega = m_get(4, 4); // Declare a 4x4 matrix
		float c = cos(0.5 * n * dt);
		float s = sin(0.5 * n * dt) / n;
		float x = w->ve[0] * s;
		float y = w->ve[1] * s;
		float z = w->ve[2] * s;

		// Omega = {{c, z, -y, x}, {-z, c, x, y}, {y, -x, c, z}, {-x, -y, -z, c}};
		Omega->me[0][0] = c;
		Omega->me[0][1] = z;
		Omega->me[0][2] = -y;
		Omega->me[0][3] = x;

		Omega->me[1][0] = -z;
		Omega->me[1][1] = c;
		Omega->me[1][2] = x;
		Omega->me[1][3] = y;

		Omega->me[2][0] = y;
		Omega->me[2][1] = -x;
		Omega->me[2][2] = c;
		Omega->me[2][3] = z;

		Omega->me[3][0] = -x;
		Omega->me[3][1] = -y;
		Omega->me[3][2] = -z;
		Omega->me[3][3] = c;

		// Eq.(7.39)
		// qHat = Omega * qHat;
		MAT *M1;
		M1 = m_get(4, 4);
		m_copy(k->Q, M1); // holds a copy of k->Q to simplify code expression (less syntax)

		// Matrix-vector multiplication with matrix and quaternion
		k->qHat->w = M1->me[0][0] * k->qHat->w + M1->me[1][1] * k->qHat->v[0] + M1->me[2][2] * k->qHat->v[1] + M1->me[3][3] * k->qHat->v[2];
		for (uint8_t i = 0; i < 3; i++)
		{
			k->qHat->v[i] = M1->me[i + 1][0] * k->qHat->w + M1->me[i + 1][1] * k->qHat->v[0] + M1->me[i + 1][2] * k->qHat->v[1] + M1->me[i + 1][3] * k->qHat->v[2];
		}
	}
}

// Update quaternion using error-quaternion
void UpQuaternion(struct Kalman *k, const Quaternion *q, const VEC *delAlpha)
{
	MAT *Xi;
	Xi = m_get(4, 3);

	Xi->me[0][0] = q->v[2];
	Xi->me[0][1] = -q->v[1];
	Xi->me[0][2] = q->v[0];

	Xi->me[1][0] = q->v[1];
	Xi->me[1][1] = q->v[2];
	Xi->me[1][2] = q->w;

	Xi->me[2][0] = -q->v[0];
	Xi->me[2][1] = q->w;
	Xi->me[2][2] = q->v[2];

	Xi->me[3][0] = -q->w;
	Xi->me[3][1] = -q->v[0];
	Xi->me[3][2] = -q->v[1];

	// Eq.(7.34)
	VEC *V1;
	V1 = v_get(4); // holds Xi * delAlpha
	VEC *V2;
	V2 = v_get(4); // holds 0.5 * Xi * delAlpha
	mv_mlt(Xi, delAlpha, V1);
	sv_mlt(0.5, V1, V2);
	// Quaternion adding
	k->qHat->w = q->w + V2->ve[0];
	for (uint8_t i = 1; i < 3; i++)
	{
		k->qHat->v[i] = q->v[i] + V2->ve[i + 1];
	}

	Quaternion_normalize(k->qHat, k->qHat); // Quaternion normalizations
}

// Discrete propagation of the state error covariance
void PropStateCovariance(struct Kalman *k, const MAT *P, const MAT *YQY_, const VEC *w, float dt)
{
	float OmegaTol = 1e-5;
	float n = (float)v_norm1(w); // v_norm1() returns a double -> cast to a float
	MAT *Phi;
	Phi = m_get(6, 6);

	if (n > OmegaTol)
	{
		// Eq.(7.45)
		float s = sin(n * dt);
		float c = cos(n * dt);
		float p = s / n;
		float q = (1 - c) / (n * n);
		float r = (n * dt - s) / (n * n * n);

		// Declare 3x3 matricies
		MAT *wX;
		wX = m_get(3, 3);
		MAT *I;
		I = m_get(3, 3);
		m_ident(I);
		MAT *wXwX;
		wXwX = m_get(3, 3);

		// setting wX
		wX->me[0][0] = 0;
		wX->me[0][1] = -w->ve[2];
		wX->me[0][2] = w->ve[1];

		wX->me[1][0] = w->ve[2];
		wX->me[1][1] = 0;
		wX->me[1][2] = -w->ve[0];

		wX->me[2][0] = -w->ve[1];
		wX->me[2][1] = w->ve[0];
		wX->me[2][2] = 0;

		m_mlt(wX, wX, wXwX);

		// Declare 3x3 matricies
		MAT *Phi_00;
		Phi_00 = m_get(3, 3);
		MAT *Phi_01;
		Phi_01 = m_get(3, 3);

		MAT *A1;
		A1 = m_get(3, 3);
		A1 = sm_mlt(p, wX, MNULL);
		MAT *A2;
		A2 = m_get(3, 3);
		sm_mlt(q, wXwX, MNULL);
		MAT *A3;
		A3 = m_get(3, 3);
		A3 = m_add(A1, A2, MNULL);
		Phi_00 = m_sub(I, A3, MNULL);

		A1 = sm_mlt(dt, I, MNULL);
		A2 = sm_mlt(r, wXwX, MNULL);
		A3 = m_sub(A1, A2, MNULL);
		MAT *A4;
		A4 = m_get(3, 3);
		A4 = sm_mlt(q, wX, MNULL);
		Phi_01 = m_sub(A4, A3, MNULL);

		Phi->me[0][0] = Phi_00->me[0][0];
		Phi->me[0][3] = Phi_01->me[0][0];
		Phi->me[0][1] = Phi_00->me[0][1];
		Phi->me[0][4] = Phi_01->me[0][1];
		Phi->me[0][2] = Phi_00->me[0][2];
		Phi->me[0][5] = Phi_01->me[0][2];
		Phi->me[1][0] = Phi_00->me[1][0];
		Phi->me[1][3] = Phi_01->me[1][0];
		Phi->me[1][1] = Phi_00->me[1][1];
		Phi->me[1][4] = Phi_01->me[1][1];
		Phi->me[1][2] = Phi_00->me[1][2];
		Phi->me[1][5] = Phi_01->me[1][2];
		Phi->me[2][0] = Phi_00->me[2][0];
		Phi->me[2][3] = Phi_01->me[2][0];
		Phi->me[2][1] = Phi_00->me[2][1];
		Phi->me[2][4] = Phi_01->me[2][1];
		Phi->me[2][2] = Phi_00->me[2][2];
		Phi->me[2][5] = Phi_01->me[2][2];
	}
	else
	{
		// Steady-state Phi: Eq.(7.45) and Eq.(7.51b)
		Phi->me[0][0] = 1;
		Phi->me[0][3] = -dt;
		Phi->me[1][1] = 1;
		Phi->me[1][4] = -dt;
		Phi->me[2][2] = 1;
		Phi->me[2][5] = -dt;
	}

	// Phi_11
	Phi->me[3][3] = 1;
	Phi->me[4][4] = 1;
	Phi->me[5][5] = 1;

	MAT *A1;
	A1 = m_get(3, 3);
	A1 = m_transp(Phi, MNULL); // holds trans(Phi)
	MAT *A2;
	A2 = m_get(3, 3);
	A2 = m_mlt(Phi, P, MNULL); // holds Phi * P
	MAT *A3;
	A3 = m_get(3, 3);
	A3 = m_mlt(A2, A1, MNULL); // holds Phi * P * trans(Phi)

	k->pHat = m_add(A3, YQY_, MNULL);
}

void getDCM(struct Kalman *k, MAT *DCM)
{
	float q0 = k->qHat->v[2];
	float q1 = k->qHat->w;
	float q2 = k->qHat->v[0];
	float q3 = k->qHat->v[1];

	DCM->me[0][0] = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
	DCM->me[0][1] = 2.0 * (q1 * q2 + q0 * q3);
	DCM->me[0][2] = 2.0 * (q1 * q3 - q0 * q2);
	DCM->me[1][0] = 2.0 * (q1 * q2 - q0 * q3);
	DCM->me[1][1] = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
	DCM->me[1][2] = 2.0 * (q2 * q3 + q0 * q1);
	DCM->me[2][0] = 2.0 * (q1 * q2 + q0 * q2);
	DCM->me[2][1] = 2.0 * (q2 * q3 - q0 * q1);
	DCM->me[2][2] = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
}
