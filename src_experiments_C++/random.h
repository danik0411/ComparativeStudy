/*
 * random.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include "standardHeaders.h"
#include "utilities.h"
#include "matrixOperations.h"

/* rnd_int: generates a random integer
   between 0 and N-1 (inclusive) */
// For ran1
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)

// For ran2
#define IM1 2147483563
#define IM2 2147483399
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV1 (1+IMM1/NTAB)

// For both ran1 & ran2
#define NTAB 32
#define AM (1.0/IM)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define rnd_int(N) \
  rand()%(N)
void UniformArray(double* a, int n, long* seed);
void GaussianArray(double* a, int n, double mu, double sigma, long* seed);
float ran1(long* idum);
float ran2(long* idum);
float randn(long* idum);
void shuffleArray(int N, long* seed, int* array);

double ncdf(double x);

#endif /* RANDOM_H_ */
