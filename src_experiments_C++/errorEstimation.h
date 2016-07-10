/*
 * errorEstimation.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef ERRORESTIMATION_H_
#define ERRORESTIMATION_H_

#include "standardHeaders.h"
#include "utilities.h"
#include "featureSelection.h"
#include "classifiers.h"
#include "random.h"
#include "main.h"

double ldaBolster(double** X, int* y, int N, int d, int* ind, LDA lda);
double qdaBolster(double** X, int* y, int N, int d, int* ind, QDA qda);
double lsvmBolster(double** X, int *y, int N, int d, int* ind, int N_bolster, long* seed);
double ksvmBolster(double** X, int *y, int N, int d, int* ind, int N_bolster, long* seed);

double ldaLOO(double** X, int* y, int N, int D, int d, double prior);
double qdaLOO(double** X, int* y, int N, int D, int d, double prior);
double lsvmLOO(double** X, int* y, int N, int D, int d);
double ksvmLOO(double** X, int* y, int N, int D, int d);

double ldaCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double qdaCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double lsvmCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double ksvmCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);

double ldaCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double qdaCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double lsvmCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);
double ksvmCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed);

double ldaBoot632(double** X, int* y, int N, int D, int d, int B, double resub, double prior, long* seed);
double qdaBoot632(double** X, int* y, int N, int D, int d, int B, double resub, double prior, long* seed);
double lsvmBoot632(double** X, int* y, int N, int D, int d, int B, double resub, long* seed);
double ksvmBoot632(double** X, int* y, int N, int D, int d, int B, double resub, long* seed);
#endif /*ERRORESTIMATION_H_*/
