/*
 * calculateErrors.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "calculateErrors.h"
#include "print.h"

#ifdef CALC_ROC
#include "ROC.h"
#endif


#ifdef REC_TIME
#define MY_CLCKS ( ( (double) CLOCKS_PER_SEC )/1000 )

double calculateTime( double t1, double t2, double REF) {

  return ( (t1-t2)/REF );

}

#endif

void calculateErrors (
		SimulationData data_trn,
		SimulationData data_tst,
		int N_trn,
		int N_tst,
		int D,
		int d,
		double prior,
		long* seed,
		AllErrors* all_errors,
		AllTimes* all_trn_times,
		AllTimes* all_tst_times,
		AllAUC* all_auc,
		int* best_features,
		int serd_group,
		int serd_elingr,
		int serd_cov_group,
		int serd_cov_elingr
		      )
		//double* lda_true,
		//double* qda_true,
		//double* lsvm_true,
		//double* ksvm_true,
		//double* dlda_true,
		//double* edc_true,
		//double* g13_true,
		//double* serd_true,
		//double* zar_true,
		//double* rlda_true )
//		
{

  double dinit = 0.00;
  
  int i;
  
  //ModelParams params;
  	
  int N_bolster = 25; // number of bolster sample point
  int K = 5; // k-fold CV
  int R = 1; // number of iterations in CV
  int B = 100; // number of bootstrap replicates
  double priorT=0.5;
  //best_features = new int [d];
  	
  
  DLDA model_DLDA;
  G13 model_G13;
  EDC model_EDC;
  LDA model_LDA;
  //QDA model_QDA;	
  SERD model_SERD_cov;
  SERD model_SERD;
  ZAR model_ZAR;
  RLDA model_RLDA;

  regparam gammas;

  double gammaBase = 1.318256739; // 1000^{1/25}
  gammas.amount = 51;
  gammas.gamma = new double [gammas.amount];
  for (i = 0; i < gammas.amount; i++){
    gammas.gamma[i] = pow(gammaBase,(double)(i-gammas.amount/2));
    //printf("i is %d\tgammas.amount/2 is %f\tgamma is %f\n", i, (double)(i-gammas.amount/2), gammas.gamma[i]);
  }

  
  //printf("gammas.amount %d\n", gammas.amount);

  model_SERD_cov.a = new double [serd_cov_group];
  model_SERD.a = new double [serd_group];
  model_DLDA.a = new double [d];
  model_G13.a = new double [d];
  model_EDC.a = new double [d];
  model_LDA.a = new double [d];
  model_ZAR.a = new double [d];
  model_RLDA.a = new double [d];
	//model_QDA.a = make_2D_matrix(d, d, dinit);
	//model_QDA.b = new double [d];

  	
  svm_model *model_LSVM;	//svm training model
  svm_node *subdata_LSVM;	//svm training data
  svm_problem subcl_LSVM;	//svm training data structure
  	
  svm_model *model_KSVM;	//svm training model
  svm_node *subdata_KSVM;	//svm training data
  svm_problem subcl_KSVM;	//svm training data structure

  clock_t tinit;

  int N_0 = 0, N_1 = 0;

  for ( int i = 0 ; i < data_tst.N ; i++ )
    if ( data_tst.labels[i] == 0 )
      N_0++;
    else
      N_1++;

  //printf("N_0 is %d\tN_1 is %d\n", N_0, N_1 );
  	
  	
  //featureSelection(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, best_features);

#ifndef REC_TIME
  rldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &gammas, &model_RLDA);
  	
  //printf("\n\nlda training\n");
  ldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_LDA);

  //printf("\n\ndlda training\n");
  dldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_DLDA);
   
  //printf("\n\nedc training\n");
  edcTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_EDC);
   
  //printf("\n\ngirko training\n");
  g13Trn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_G13);

  //printf("starting serd training\n\n");

  // without log
  //printf("serd_groups is %d\tserd_elingr is %d\n",serd_group,serd_elingr);
  serdTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_SERD, serd_group, serd_elingr);

  //printf("serd without cov has been finished\n\n");
   
  //printf("\n\nserd training\n");
  // with log and covariance matrix
  serdCovTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_SERD_cov, serd_cov_group, serd_cov_elingr);

  //printf("serd with cov has been finished\n\n");
  
   

  //printf("\n\nzar training\n");
  zarTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_ZAR);
  
  //printf("\n\nrlda training\n");

  //qdaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_QDA);
  //printf("qdaTrn\n");
  	
  //printf("\n\nlsvmTrn\n\n");
  model_LSVM = svmTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, 0, &subdata_LSVM, &subcl_LSVM);
  	
  //printf("\n\nksvmTrn\n\n");
  model_KSVM = svmTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, 2, &subdata_KSVM, &subcl_KSVM);

  
  (*all_errors).lda_true_error = ldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_LDA);
  //(*all_errors).qda_true_error = qdaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_QDA);
  (*all_errors).lsvm_true_error = svmTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_LSVM);
  (*all_errors).ksvm_true_error = svmTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_KSVM);

  (*all_errors).dlda_true_error = dldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_DLDA);
  (*all_errors).edc_true_error = edcTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_EDC);
  (*all_errors).g13_true_error = g13Tst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_G13);
  (*all_errors).serd_true_error = serdTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_SERD, serd_group, serd_elingr);
  (*all_errors).serd_true_error_cov = serdTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_SERD_cov, serd_cov_group, serd_cov_elingr);
  (*all_errors).zar_true_error = zarTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_ZAR);
  (*all_errors).rlda_true_error = rldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_RLDA);


#endif



#ifdef REC_TIME
  tinit = clock();
  rldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &gammas, &model_RLDA);
  (*all_trn_times).rlda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);
  	
  //printf("\n\nlda training\n");
  tinit = clock();
  ldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_LDA);
  (*all_trn_times).lda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  //printf("\n\ndlda training\n");
  tinit = clock();
  dldaTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_DLDA);
  (*all_trn_times).dlda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);
   
  //printf("\n\nedc training\n");
  tinit = clock();
  edcTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_EDC);
  (*all_trn_times).edc_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);
  
   
  //printf("\n\ngirko training\n");
  tinit = clock();
  g13Trn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_G13);
  (*all_trn_times).g13_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  //printf("starting serd training\n\n");

  // without log
  //printf("serd_groups is %d\tserd_elingr is %d\n",serd_group,serd_elingr);
  tinit = clock();
  serdTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_SERD, serd_group, serd_elingr);
  (*all_trn_times).serd_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  //printf("serd without cov has been finished\n\n");
   
  //printf("\n\nserd training\n");
  // with log and covariance matrix
  tinit = clock();
  serdCovTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_SERD_cov, serd_cov_group, serd_cov_elingr);
  (*all_trn_times).serd_time_cov = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  //printf("serd with cov has been finished\n\n");

  //printf("\n\nzar training\n");
  tinit = clock();
  zarTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, prior, &model_ZAR);
  (*all_trn_times).zar_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);
  
  //printf("\n\nlsvmTrn\n\n");
  tinit = clock();
  model_LSVM = svmTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, 0, &subdata_LSVM, &subcl_LSVM);
  (*all_trn_times).lsvm_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);
  	
  //printf("\n\nksvmTrn\n\n");
  tinit = clock();
  model_KSVM = svmTrn(data_trn.data, data_trn.labels, data_trn.N, d, best_features, 2, &subdata_KSVM, &subcl_KSVM);
  (*all_trn_times).ksvm_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);


  // testing

  
  tinit = clock();
  (*all_errors).lda_true_error = ldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_LDA);
  (*all_tst_times).lda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).lsvm_true_error = svmTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_LSVM);
  (*all_tst_times).lsvm_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).ksvm_true_error = svmTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_KSVM);
  (*all_tst_times).ksvm_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).dlda_true_error = dldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_DLDA);
  (*all_tst_times).dlda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).edc_true_error = edcTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_EDC);
  (*all_tst_times).edc_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).g13_true_error = g13Tst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_G13);
  (*all_tst_times).g13_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).serd_true_error = serdTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_SERD, serd_group, serd_elingr);
  (*all_tst_times).serd_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).serd_true_error_cov = serdTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_SERD_cov, serd_cov_group, serd_cov_elingr);
  (*all_tst_times).serd_time_cov = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).zar_true_error = zarTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, prior, model_ZAR);
  (*all_tst_times).zar_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

  tinit = clock();
  (*all_errors).rlda_true_error = rldaTst(data_tst.data, data_tst.labels, data_tst.N, d, best_features, priorT, model_RLDA);
  (*all_tst_times).rlda_time = calculateTime( (double) clock(), (double) tinit, MY_CLCKS);

#endif


  //printf("\n\ntrying calc roc\n\n");
  
// calculating the true false
#ifdef CALC_ROC

  //printf("\n\ncalc roc\n\n");

  linModel linear_LDA      = convertLinModel<LDA>  ( model_LDA , d );
  linModel linear_DLDA     = convertLinModel<DLDA> ( model_DLDA, d );
  linModel linear_RLDA     = convertLinModel<RLDA> ( model_RLDA, d );

  linModel linear_EDC      = convertLinModel<EDC>  ( model_EDC , d );
  linModel linear_ZAR      = convertLinModel<ZAR>  ( model_ZAR , d );

  linModel linear_G13      = convertLinModel<G13>  ( model_G13 , d );

  linModel linear_SERD     = convertLinModel<SERD> ( model_SERD     , serd_group     );
  linModel linear_SERD_cov = convertLinModel<SERD> ( model_SERD_cov , serd_cov_group );


  //printf("\n\nconversion has been finished\n\n");


  int *elab_LDA      = new int[data_tst.N]();
  int *elab_DLDA     = new int[data_tst.N]();
  int *elab_RLDA     = new int[data_tst.N]();

  int *elab_EDC      = new int[data_tst.N]();
  int *elab_ZAR      = new int[data_tst.N]();
  int *elab_G13      = new int[data_tst.N]();
  
  int *elab_SERD     = new int[data_tst.N]();
  int *elab_SERD_cov = new int[data_tst.N]();

  int *elab_lsvm     = new int[data_tst.N]();
  int *elab_ksvm     = new int[data_tst.N]();

  
  double *eval_LDA      = new double[data_tst.N]();
  double *eval_DLDA     = new double[data_tst.N]();
  double *eval_RLDA     = new double[data_tst.N]();

  double *eval_EDC      = new double[data_tst.N]();
  double *eval_ZAR      = new double[data_tst.N]();
  double *eval_G13      = new double[data_tst.N]();

  double *eval_SERD     = new double[data_tst.N]();
  double *eval_SERD_cov = new double[data_tst.N]();

  double *eval_lsvm     = new double[data_tst.N]();
  double *eval_ksvm     = new double[data_tst.N]();

  int *temp_ind      = new int[data_tst.N]();

  int v0 = data_tst.labels[0];
  int v1 = data_tst.labels[data_tst.N-1];

  //int pos = data_tst.labels[data_tst.N-1];
  //int neg = data_tst.labels[0];


  for ( int ii = 0 ; ii < data_tst.N ; ii++ )
    temp_ind[ii] = ii;

  //printf("\n\ndeclaration has been finished\n\n");

  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_LDA,  eval_LDA  , elab_LDA  );
  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_DLDA, eval_DLDA , elab_DLDA );
  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_RLDA, eval_RLDA , elab_RLDA );

  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_EDC,  eval_EDC  , elab_EDC  );
  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_G13,  eval_G13  , elab_G13  );
  linClassTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_ZAR,  eval_ZAR  , elab_ZAR  );

  serdTrueFalse (data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_SERD     , serd_group     , serd_elingr     , eval_SERD     , elab_SERD     );
  serdTrueFalse (data_tst.data, data_tst.labels, data_tst.N, d, best_features, linear_SERD_cov , serd_cov_group , serd_cov_elingr , eval_SERD_cov , elab_SERD_cov );

  svmTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, model_LSVM, eval_lsvm, elab_lsvm );
  svmTrueFalse( data_tst.data, data_tst.labels, data_tst.N, d, best_features, model_KSVM, eval_ksvm, elab_ksvm );


  //                        calc_auc( est_value     , true_label     , est_label     , N         , N0 , N1 , v0, v1
  /*
  (*all_auc).lda_auc      = calc_auc( eval_LDA      , data_tst.labels, elab_LDA      , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).dlda_auc     = calc_auc( eval_DLDA     , data_tst.labels, elab_DLDA     , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).rlda_auc     = calc_auc( eval_RLDA     , data_tst.labels, elab_RLDA     , data_tst.N, N_0, N_1, v0, v1 );

  (*all_auc).g13_auc      = calc_auc( eval_G13      , data_tst.labels, elab_G13      , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).zar_auc      = calc_auc( eval_ZAR      , data_tst.labels, elab_ZAR      , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).edc_auc      = calc_auc( eval_EDC      , data_tst.labels, elab_EDC      , data_tst.N, N_0, N_1, v0, v1 );

  (*all_auc).serd_auc     = calc_auc( eval_SERD     , data_tst.labels, elab_SERD     , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).serd_auc_cov = calc_auc( eval_SERD_cov , data_tst.labels, elab_SERD_cov , data_tst.N, N_0, N_1, v0, v1 );
  													     
  (*all_auc).lsvm_auc     = calc_auc( eval_lsvm     , data_tst.labels, elab_lsvm     , data_tst.N, N_0, N_1, v0, v1 );
  (*all_auc).ksvm_auc     = calc_auc( eval_ksvm     , data_tst.labels, elab_ksvm     , data_tst.N, N_0, N_1, v0, v1 );
  */


  (*all_auc).lda_auc      = aucroc2(eval_LDA, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).dlda_auc     = aucroc2(eval_DLDA, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).rlda_auc     = aucroc2(eval_RLDA, data_tst.labels, data_tst.N, N_0, N_1);
   
  (*all_auc).g13_auc      = aucroc2(eval_G13, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).zar_auc      = aucroc2(eval_ZAR, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).edc_auc      = aucroc2(eval_EDC, data_tst.labels, data_tst.N, N_0, N_1);
   
  (*all_auc).serd_auc     = aucroc2(eval_SERD, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).serd_auc_cov = aucroc2(eval_SERD_cov, data_tst.labels, data_tst.N, N_0, N_1);
   
  (*all_auc).lsvm_auc     = aucroc2(eval_lsvm, data_tst.labels, data_tst.N, N_0, N_1);
  (*all_auc).ksvm_auc     = aucroc2(eval_ksvm, data_tst.labels, data_tst.N, N_0, N_1);

  //printf("ksvm 0 neg %d\t pos %d\n", countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 1, -1 ));
  //printf("ksvm 1 neg %d\t pos %d\n", countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 1, 1 ));

  /*
  printf("edc auc\n");
  (*all_auc).edc_auc      = aucroc2(eval_EDC, data_tst.labels, data_tst.N, N_0, N_1);
  printf("edc 0 neg %d\t pos %d\n", countLabels2(eval_EDC, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_EDC, data_tst.labels, data_tst.N, 1, -1 ));
  printf("edc 1 neg %d\t pos %d\n", countLabels2(eval_EDC, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_EDC, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));

  printf("\ndlda auc\n");
  printf("dlda old auc %.3f\tnew auc %.3f\n",     (*all_auc).dlda_auc, aucroc2( eval_DLDA, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("DLDA 0 neg %d\t pos %d\n", countLabels2(eval_DLDA, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_DLDA, data_tst.labels, data_tst.N, 1, -1 ));
  printf("DLDA 1 neg %d\t pos %d\n", countLabels2(eval_DLDA, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_DLDA, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_DLDA, 1, data_tst.N), countLabels(elab_DLDA, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_DLDA, -1, data_tst.N), countValues(eval_DLDA, 1, data_tst.N));
  //
  printf("\ng13 auc\n");
  printf("g13 old auc %.3f\tnew auc %.3f\n",      (*all_auc).g13_auc, aucroc2( eval_G13, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("G1G13 neg %d\t pos %d\n", countLabels2(eval_G13, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_G13, data_tst.labels, data_tst.N, 1, -1 ));
  printf("G13 1 neg %d\t pos %d\n", countLabels2(eval_G13, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_G13, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_G13, 1, data_tst.N), countLabels(elab_G13, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_G13, -1, data_tst.N), countValues(eval_G13, 1, data_tst.N));
  // 
  printf("\nzar auc\n");
  printf("zar old auc %.3f\tnew auc %.3f\n",      (*all_auc).zar_auc, aucroc2( eval_ZAR, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("ZAR 0 neg %d\t pos %d\n", countLabels2(eval_ZAR, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_ZAR, data_tst.labels, data_tst.N, 1, -1 ));
  printf("ZAR 1 neg %d\t pos %d\n", countLabels2(eval_ZAR, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_ZAR, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_ZAR, 1, data_tst.N), countLabels(elab_ZAR, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_ZAR, -1, data_tst.N), countValues(eval_ZAR, 1, data_tst.N));
  // 

  printf("\nedc auc\n");
  printf("edc old auc %.3f\tnew auc %.3f\n",      (*all_auc).edc_auc, aucroc2( eval_EDC, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("edc 0 neg %d\t pos %d\n", countLabels2(eval_EDC, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_EDC, data_tst.labels, data_tst.N, 1, -1 ));
  printf("edc 1 neg %d\t pos %d\n", countLabels2(eval_EDC, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_EDC, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_EDC, 1, data_tst.N), countLabels(elab_EDC, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_EDC, -1, data_tst.N), countValues(eval_EDC, 1, data_tst.N));
  // 
  printf("\nserd auc\n");
  printf("serd old auc %.3f\tnew auc %.3f\n",     (*all_auc).serd_auc, aucroc2( eval_SERD, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("SERD 0 neg %d\t pos %d\n", countLabels2(eval_SERD, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_SERD, data_tst.labels, data_tst.N, 1, -1 ));
  printf("SERD 1 neg %d\t pos %d\n", countLabels2(eval_SERD, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_SERD, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_SERD, 1, data_tst.N), countLabels(elab_SERD, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_SERD, -1, data_tst.N), countValues(eval_SERD, 1, data_tst.N));
  // 
  printf("\nserd_cov auc\n");
  printf("serd_cov old auc %.3f\tnew auc %.3f\n", (*all_auc).serd_auc_cov, aucroc2( eval_SERD_cov, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("SERD_cov 0 neg %d\t pos %d\n", countLabels2(eval_SERD_cov, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_SERD_cov, data_tst.labels, data_tst.N, 1, -1 ));
  printf("SERD_cov 1 neg %d\t pos %d\n", countLabels2(eval_SERD_cov, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_SERD_cov, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_SERD_cov, 1, data_tst.N), countLabels(elab_SERD_cov, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_SERD_cov, -1, data_tst.N), countValues(eval_SERD_cov, 1, data_tst.N));
  // 
  
  printf("\nrlda auc\n");
  printf("rlda old auc %.3f\tnew auc %.3f\n",     (*all_auc).rlda_auc, aucroc2( eval_RLDA, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("RLDA 0 neg %d\t pos %d\n", countLabels2(eval_RLDA, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_RLDA, data_tst.labels, data_tst.N, 1, -1 ));
  printf("RLDA 1 neg %d\t pos %d\n", countLabels2(eval_RLDA, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_RLDA, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_RLDA, 1, data_tst.N), countLabels(elab_RLDA, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_RLDA, -1, data_tst.N), countValues(eval_RLDA, 1, data_tst.N));
  // 
  printf("\nlsvm auc\n");
  printf("lsvm old auc %.3f\tnew auc %.3f\n",     (*all_auc).lsvm_auc, aucroc2( eval_lsvm, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("lsvm 0 neg %d\t pos %d\n", countLabels2(eval_lsvm, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_lsvm, data_tst.labels, data_tst.N, 1, -1 ));
  printf("lsvm 1 neg %d\t pos %d\n", countLabels2(eval_lsvm, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_lsvm, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_lsvm, 1, data_tst.N), countLabels(elab_lsvm, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_lsvm, -1, data_tst.N), countValues(eval_lsvm, 1, data_tst.N));
  // 
  
  printf("\nksvm auc\n");
  printf("ksvm old auc %.3f\tnew auc %.3f\n",     (*all_auc).ksvm_auc, aucroc2( eval_ksvm, data_tst.labels, data_tst.N, N_0, N_1 ));
  printf("ksvm 0 neg %d\t pos %d\n", countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 0, -1 ), countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 1, -1 ));
  printf("ksvm 1 neg %d\t pos %d\n", countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 0, 1 ), countLabels2(eval_ksvm, data_tst.labels, data_tst.N, 1, 1 ));
  printf("true zeros is %d\t ones is %d\n", countLabels(data_tst.labels, 0, data_tst.N), countLabels(data_tst.labels, 1, data_tst.N));
  printf("est negs = %d\tposs = %d\n", countLabels(elab_ksvm, 1, data_tst.N), countLabels(elab_ksvm, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_ksvm, -1, data_tst.N), countValues(eval_ksvm, 1, data_tst.N));
  
  //printf("lda old auc %.3f\tnew auc %.3f\n",      (*all_auc).lda_auc, aucroc2( eval_LDA, data_tst.labels, data_tst.N, N_0, N_1 ));
  //printf("est negs = %d\tposs = %d\n", countLabels(elab_ksvm, 1, data_tst.N), countLabels(elab_ksvm, 0, data_tst.N));
  //printf("est negv = %d\tposv = %d\n", countValues(eval_ksvm, -1, data_tst.N), countValues(eval_ksvm, 1, data_tst.N));
  */
  

  /*
  printf("ksvm_val is\n");
  a_double_toScreen(eval_ksvm,data_tst.N);
  printf("ksvm_lab is\n");
  a_int_toScreen(elab_ksvm,data_tst.N);
  printf("true_lab is\n");
  a_int_toScreen(data_tst.labels,data_tst.N);
  printf("\n\n");
  */


  //printf("ksvm\t\tedc\t\t\n");
  //for ( int ii = 0 ; ii < data_tst.N ; ii++ ){
  //printf("%.3f\t%d\t%.3f\t%d\n", eval_ksvm[ii], elab_ksvm[ii], eval_EDC[ii], elab_EDC[ii]);
  //}



  delete[] linear_LDA.a;
  delete[] linear_DLDA.a;
  delete[] linear_RLDA.a;

  delete[] linear_EDC.a;
  delete[] linear_G13.a;
  delete[] linear_ZAR.a;

  delete[] linear_SERD.a;
  delete[] linear_SERD_cov.a;

  delete[] temp_ind;

  delete[] elab_LDA      ;
  delete[] elab_DLDA     ;
  delete[] elab_RLDA     ;

  delete[] elab_EDC      ;
  delete[] elab_ZAR      ;
  delete[] elab_G13      ;

  delete[] elab_SERD     ;
  delete[] elab_SERD_cov ;

  delete[] elab_lsvm     ;
  delete[] elab_ksvm     ;
  

  delete[] eval_LDA      ;
  delete[] eval_DLDA     ;
  delete[] eval_RLDA     ;

  delete[] eval_EDC      ;
  delete[] eval_ZAR      ;
  delete[] eval_G13      ;
  
  delete[] eval_SERD     ;
  delete[] eval_SERD_cov ;

  delete[] eval_lsvm     ;
  delete[] eval_ksvm     ;

  

#endif
  

	// //printf("==========");
	//(*all_errors).lda_resub_error = ldaTst(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_LDA);
	//(*all_errors).qda_resub_error = qdaTst(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_QDA);
	//(*all_errors).lsvm_resub_error = svmTst(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_LSVM);
	//(*all_errors).ksvm_resub_error = svmTst(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_KSVM);
  	
	//(*all_errors).lda_bolster_error = ldaBolster(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_LDA);
	//(*all_errors).qda_bolster_error = qdaBolster(data_trn.data, data_trn.labels, data_trn.N, d, best_features, model_QDA);
	//(*all_errors).lsvm_bolster_error = lsvmBolster(data_trn.data, data_trn.labels, data_trn.N, d, best_features, N_bolster, seed);
	//(*all_errors).ksvm_bolster_error = ksvmBolster(data_trn.data, data_trn.labels, data_trn.N, d, best_features, N_bolster, seed);
  	
	//(*all_errors).lda_loo_error = ldaLOO(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d);
	//(*all_errors).qda_loo_error = qdaLOO(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d);
	//(*all_errors).lsvm_loo_error = lsvmLOO(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d);
	//(*all_errors).ksvm_loo_error = ksvmLOO(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d);
  	
	//(*all_errors).lda_cvkfold_error = ldaCVkFold(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("ldaCVkFold\n");
	//(*all_errors).qda_cvkfold_error = qdaCVkFold(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("qdaCVkFold\n");
	//(*all_errors).lsvm_cvkfold_error = lsvmCVkFold(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("svmCVkFold\n");
	//(*all_errors).ksvm_cvkfold_error = ksvmCVkFold(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("ksvmCVkFold\n");
	//(*all_errors).lsvm_cvkfold_error = 0.0;
	//(*all_errors).ksvm_cvkfold_error = ksvmCVkFold(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, seed);
  	
	//(*all_errors).lda_cvkfold_Mixture_error = ldaCVkFold_Mixture(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("ldaCVkFold_Mixture\n");
	//(*all_errors).qda_cvkfold_Mixture_error = qdaCVkFold_Mixture(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("qdaCVkFold_Mixture\n");
	//(*all_errors).lsvm_cvkfold_Mixture_error = lsvmCVkFold_Mixture(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("lsvmCVkFold_Mixture\n");
	//(*all_errors).ksvm_cvkfold_Mixture_error = ksvmCVkFold_Mixture(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, prior, seed);
	//printf("ksvmCVkFold_Mixture\n");
	//(*all_errors).lsvm_cvkfold_Mixture_error = 0.0;
	//(*all_errors).ksvm_cvkfold_Mixture_error = ksvmCVkFold_Mixture(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, K, R, seed);
  	
	//(*all_errors).lda_boot632_error = ldaBoot632(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, B, (*all_errors).lda_resub_error, seed);
	//(*all_errors).qda_boot632_error = qdaBoot632(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, B, (*all_errors).qda_resub_error, seed);
	//(*all_errors).lsvm_boot632_error = lsvmBoot632(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, B, (*all_errors).lsvm_resub_error, seed);
	//(*all_errors).ksvm_boot632_error = ksvmBoot632(data_trn.data, data_trn.labels, data_trn.N, data_trn.D, d, B, (*all_errors).ksvm_resub_error, seed);
  	
  //printf("best_features is\n");
  //for ( int iii = 0 ; iii < d ; iii++ )
  //printf("%d\t", best_features[iii]);
  //printf("\n\n");
  
  
  	
  delete [] model_LDA.a;
  delete [] model_DLDA.a;
  delete [] model_G13.a;
  delete [] model_EDC.a;
  delete [] model_ZAR.a;
  delete [] model_SERD.a;
  delete [] model_SERD_cov.a;
  delete [] model_RLDA.a;
  delete [] gammas.gamma;
  //delete_2D_matrix(d, d, model_QDA.a);
  svmDestroy(model_LSVM, subdata_LSVM, &subcl_LSVM);
  svmDestroy(model_KSVM, subdata_KSVM, &subcl_KSVM);
  //delete [] best_features;	


  return;

}

