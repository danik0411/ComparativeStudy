/*
 * dataGeneration.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef DATAGENERATION_H_
#define DATAGENERATION_H_

#include "standardHeaders.h"
#include "utilities.h"
#include "matrixOperations.h"
#include "random.h"

struct SimulationData {
  double** data;
  int* labels; /* labels for samples */
  int N;
  int D;
};

void dataGeneration(double** X,
		int* y,
		int N_trn,
		int N_tst,
		int D,
		long* seed,
		SimulationData* data_trn,
		SimulationData* data_tst);

void dataLogfy(double** X, int N, int D, double** Xlog);
void dataStandardize(double** X, int N, int D, double** Xnew);

#endif /* DATAGENERATION_H_ */
