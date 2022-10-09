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
*	Written by: Nathan Tang
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
//#include <iostream>

// Required files
#include "Kalman.h"
#include "matrix.h"												// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "matrix2.h"												// Courtesy of David E. Steward & Zbigniew Leyk - Copyright (C) 1993
#include "quaternion/Quaternion.h"												// Courtesy of Martin Weigel					- Copyright (C) 2019


int main(void)
{
    Quaternion *qInit; Quaternion_set(0, 0, 0, 1, &qInit); // Initial quaternion
    VEC *bInit;    bInit = v_get(3);
    bInit->ve[0] = 1.324;   bInit->ve[1] = 2.4532;  bInit->ve[2] = 3.4532; // inital bias
    MAT *pInit;     pInit = m_get(6,6);     pInit = m_ident(MNULL);
    float dt = 0.134; // Time step
    
    // Noise standard deviation
    float sigmaR[2] = {1.345, 2.654}; // Vector measurements noise
    float sigmaG[3] = {4.453, 2.534, 3.234}; // Gyro bias
    float sigmaB[3] = {2.435, 1.2345, 4.55}; // Gyro-bias noise

    // Vector measurements
    MAT *mb;    mb = m_get(3,2);
    // Eg. Magnetometer     // Eg. Accelerometer
    mb->me[0][0] = 1.32;    mb->me[0][1] = 4.6534;
    mb->me[1][0] = 2.65;    mb->me[1][1] = 5.634;
    mb->me[2][0] = 3.63;    mb->me[2][1] = 6.634;

    // Reference vectors
    MAT *mr;    mr = m_get(3,2);
    // Magnetometer vector  // Acc. due to gravity
    // in inertial frame    // in inertial frame
    mr->me[0][0] = 7.53;    mr->me[0][1] = 1.42;
    mr->me[1][0] = 8.23;    mr->me[1][1] = 2.432;
    mr->me[2][0] = 9.97;    mr->me[2][1] = 3.546;

    // Gyro reading
    VEC *w;     w = v_get(3);
    w->ve[0] = 0.223;   w->ve[1] = 2.34;    w->ve[2] = 3.45;

    struct Kalman *k; // Declare struct to hold MEKF members

    MEKFinit(k, qInit, bInit, pInit, sigmaR, sigmaG, sigmaB, dt); // Initialize MEKF

    Estimate(k, w, mb, mr, dt); // Run MEKF

    // Display
    printf("Estimated quaternion: %d %d %d %d\n", k->qHat->w, k->qHat->v[0], k->qHat->v[1], k->qHat->v[2]);
    printf("Estimated bias: %d %d %d %d\n", k->bHat->ve[0], k->bHat->ve[1], k->bHat->ve[2], k->bHat->ve[3]);
	
    return 0;
}