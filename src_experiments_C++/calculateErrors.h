/*
 * calculateErrors.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef CALCULATEERRORS_H_
#define CALCULATEERRORS_H_

#include "standardHeaders.h"
#include "print.h"
#include "utilities.h"
#include "matrixOperations.h"
#include "dataGeneration.h"
#include "featureSelection.h"
#include "errorEstimation.h"
#include "classifiers.h"


struct AllErrors {
  double lda_true_error;
  double edc_true_error;
  double dlda_true_error;
  double g13_true_error;
  double serd_true_error;
  double serd_true_error_cov;
  double zar_true_error;
  double rlda_true_error;
  //double qda_true_error;
  double lsvm_true_error;
  double ksvm_true_error;

  //double lda_resub_error;
  //double qda_resub_error;
  //double lsvm_resub_error;
  //double ksvm_resub_error;
  //
  //double lda_bolster_error;
  //double qda_bolster_error;
  //double lsvm_bolster_error;
  //double ksvm_bolster_error;
  //
  //double lda_loo_error;
  //double qda_loo_error;
  //double lsvm_loo_error;
  //double ksvm_loo_error;
  //
  //double lda_cvkfold_error;
  //double qda_cvkfold_error;
  //double lsvm_cvkfold_error;
  //double ksvm_cvkfold_error;
  //	
  //double lda_cvkfold_Mixture_error;
  //double qda_cvkfold_Mixture_error;
  //double lsvm_cvkfold_Mixture_error;
  //double ksvm_cvkfold_Mixture_error;
  //
  //double lda_boot632_error;
  //double qda_boot632_error;
  //double lsvm_boot632_error;
  //double ksvm_boot632_error;
};

struct AllTimes {

  double lda_time;
  double edc_time;
  double dlda_time;

  double g13_time;
  double serd_time;
  double serd_time_cov;
  double zar_time;
  double rlda_time;

  double lsvm_time;
  double ksvm_time;

};

struct AllAUC {

  double lda_auc;
  double edc_auc;
  double dlda_auc;

  double g13_auc;
  double serd_auc;
  double serd_auc_cov;
  double zar_auc;
  double rlda_auc;

  double lsvm_auc;
  double ksvm_auc;

};

#ifdef REC_TIME
double calculateTime ( double t2, double t1, double REF );
#endif

void calculateErrors(
		SimulationData data_trn,
		SimulationData data_tst,
		int N_trn,
		int N_tst,
		int D,
		int d,
		double prior,
		long* seed,
		AllErrors* all_errors,
		AllTimes* all_trn_errors,
		AllTimes* all_tst_errors,
		AllAUC* all_auc,
		int* best_features,
		int serd_gourp,
		int serd_eling,
		int serd_cov_gourp,
		int serd_cov_eling);
		//double* lda_true_error,
		//double* qda_true_error,
		//double* lsvm_true_error,
		//double* ksvm_true_error,
		//double* dlda_true_error,
		//double* edc_true_error,
		//double* g13_true_error,
		//double* serd_true_error,
		//double* zar_true_error,
		//double* rlda_true_error
		//);
		


#endif /* CALCULATEERRORS_H_ */
