#ifndef FEATURESELECTION_H_
#define FEATURESELECTION_H_

#include "standardHeaders.h"
#include "matrixOperations.h"

int featureSelection(double** X, int* y, int N, int D, int d, int* best_features);
double t_test(double* feature_class0, double* feature_class1, int n_class0, int n_class1);

#endif /*FEATURESELECTION_H_*/
