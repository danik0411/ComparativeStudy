/*
 * main.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "main.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main (int argc, char *argv[])
{
	//////////////////////// MPI ///////////////////////////////////
	int MYRANK=0, NP=1;   //DEFAULT assume total one node, current is node 0

	printf("NP before USE_MPI is %d\n", NP);

#ifdef USE_MPI
	MPI_Status status;  // return status for receive
#endif

#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NP);
	MPI_Comm_rank(MPI_COMM_WORLD, &MYRANK);
	printf("NP is %d\n", NP);
#endif


#ifdef REC_TIME
	printf("i will record time readings during this run\n");
#endif


	double dinit = 0.00;
	int i, j, flag, rep;

	int N_0; // samples in class 0
	int N_1; // samples in class 1
	int N; // total number of sample points
	int N_trn;
	int N_tst;
	int ijk;
	int ihj;
	int ijkl;
	int D; // total # of features
	int d = 1; // after feature selection
	int max_rep = 1;
	int rlda_optGamma_rep = 0;
	double rlda_optGamma_sr = 0.00;
	double prior=0.5;
	double priorT=0.5;
	//for (i=0; i<3; i++)
	   //printf("n01r_vec[ijk]=%.2f\n",n01r_vec[i]);
	//double n01r;
	long int seed = 0L;
	double temp;

	double** X;
	double** Xlog;
	double** Xstan;
	
	int* y;

	AllErrors all_errors;
	AllTimes all_trn_times;
	AllTimes all_tst_times;
	AllAUC all_auc;
	
	int* best_features;
	int* best_features2;

	int serd_group;
	int serd_elingr;

	int serd_cov_group;
	int serd_cov_elingr;

	all_errors.lda_true_error = 0.00;
	all_errors.edc_true_error = 0.00;
	all_errors.g13_true_error = 0.00;
	all_errors.serd_true_error = 0.00;

	//double* lda_true_error;
	//double* qda_true_error;
	//double* lsvm_true_error;
	//double* ksvm_true_error;
	//
	//double* edc_true_error;
	//double* dlda_true_error;
	//double* g13_true_error;
	//double* serd_true_error;
	//double* zar_true_error;
	//double* rlda_true_error;

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

	char *dataFilename;
	char *outputFilename;
	

	FILE *fpIn;


	printf("parsing arguments\n");

	//===========================================//
	// Arguments are passed through command-line //
	//===========================================//
	switch (parse_common_arguments(
			&dataFilename,
			&outputFilename,
			&N_trn,
			&d,
			&max_rep,
			argc, argv))
	{
	case -1:
	printf("Error in function parse_common_arguments()!\n");
	return(-1);
	case -2:
	printf("Usage printed\n"); // should add some usage printing function here
	return(0);
	default:
		printf("\n");
	}
	
	// make sure that both values are integers
	serd_group = 1;                serd_elingr = d;
	serd_cov_group = 1;        serd_cov_elingr = d;
	

	fpIn = fopen(dataFilename, "r");
	if(fpIn==NULL)
	{
		printf("cannot open data file in node %d!\n", MYRANK);
		return(-1);
	}

	//obtain the size information of two classes
	fscanf(fpIn,"%d, %d\r\n", &N_0, &N_1);

	// total number of features
	fscanf(fpIn,"%d\r\n", &D);

	//size_class1 = size_class0;
	N = N_0 + N_1;
	N_tst = N - N_trn;

	// initializing matrix of original data set
	X = make_2D_matrix(N, D, dinit);
	Xlog = make_2D_matrix(N, D, dinit);
	Xstan = make_2D_matrix(N, D, dinit);
	y = new int [N];

	// obtain the gene data from data file (rows: fetaures, columns: patients)
	// but we transpose it to: (rows: patients, columns: features)
	for (i=0; i<N; i++)
	{
		if (i<N_0)
			y[i] = 0;
		else
			y[i] = 1;
	}

	for (i=0; i<D; i++)
	{
		for (j=0; j<N; j++)
		{
			fscanf(fpIn, "%lf\t", &temp);
			X[j][i] = temp;
		}
		fscanf(fpIn, "\r\n");
	}

	fclose(fpIn);

	srand((unsigned)time(NULL)+MYRANK);
	seed = (long int)(-abs((long)floor(rand()*time(NULL)+MYRANK+1)));
	//seed = -12314001L;
	printf("seed is %ld\n\n", seed);
		

	char outputFilenameFinal[90];
	FILE* fp_lda_true_error_c01;
	FILE* fp_dlda_true_error_c01;
	FILE* fp_rlda_true_error_c01;
	
	FILE* fp_edc_true_error_c01;
	FILE* fp_g13_true_error_c01;
	FILE* fp_serd_true_error_c01;
	FILE* fp_serd_true_error_c01_cov;
	FILE* fp_zar_true_error_c01;
	
	FILE* fp_lsvm_true_error_c01;
	FILE* fp_ksvm_true_error_c01;

	
#ifdef REC_TIME
	FILE* fp_trn_lda_time;
	FILE* fp_trn_dlda_time;
	FILE* fp_trn_rlda_time;

	FILE* fp_trn_edc_time;
	FILE* fp_trn_g13_time;
	FILE* fp_trn_serd_time;
	FILE* fp_trn_serd_time_cov;
	FILE* fp_trn_zar_time;

	FILE* fp_trn_lsvm_time;
	FILE* fp_trn_ksvm_time;

	/////////////////////////
	
	FILE* fp_tst_lda_time;
	FILE* fp_tst_dlda_time;
	FILE* fp_tst_rlda_time;

	FILE* fp_tst_edc_time;
	FILE* fp_tst_g13_time;
	FILE* fp_tst_serd_time;
	FILE* fp_tst_serd_time_cov;
	FILE* fp_tst_zar_time;

	FILE* fp_tst_lsvm_time;
	FILE* fp_tst_ksvm_time;
#endif

#ifdef CALC_ROC

	FILE* roc_lda;
	FILE* roc_dlda;
	FILE* roc_rlda;

	FILE* roc_edc;
	FILE* roc_g13;
	FILE* roc_zar;

	FILE* roc_serd;
	FILE* roc_serd_cov;

	FILE* roc_lsvm;
	FILE* roc_ksvm;

#endif
	
	
	if (MYRANK==0)
	  {
	    
	    // true error filenames
	    printf("outputFileName is %s\n", outputFilename);
	    sprintf(outputFilenameFinal, "%s_lda_true_error.txt", outputFilename);
	    fp_lda_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_lda_true_error_c01==NULL)
	      {
		printf("cannot open %s_lda_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_dlda_true_error.txt", outputFilename);
	    fp_dlda_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_dlda_true_error_c01==NULL)
	      {
		printf("cannot open %s_dlda_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_edc_true_error.txt", outputFilename);
	    fp_edc_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_edc_true_error_c01==NULL)
	      {
		printf("cannot open %s_edc_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_zar_true_error.txt", outputFilename);
	    fp_zar_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_zar_true_error_c01==NULL)
	      {
		printf("cannot open %s_zar_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_rlda_true_error.txt", outputFilename);
	    fp_rlda_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_rlda_true_error_c01==NULL)
	      {
		printf("cannot open %s_rlda_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_g13_true_error.txt", outputFilename);
	    fp_g13_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_g13_true_error_c01==NULL)
	      {
		printf("cannot open %s_g13_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_serd_true_error.txt", outputFilename);
	    fp_serd_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_serd_true_error_c01==NULL)
	      {
		printf("cannot open %s_serd_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_serd_true_error_cov.txt", outputFilename);
	    fp_serd_true_error_c01_cov = fopen(outputFilenameFinal, "w");
	    if(fp_serd_true_error_c01_cov==NULL)
	      {
		printf("cannot open %s_serd_true_error_cov.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_lsvm_true_error.txt", outputFilename);
	    fp_lsvm_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_lsvm_true_error_c01==NULL)
	      {
		printf("cannot open %s_lsvm_true_error.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_ksvm_true_error.txt", outputFilename);
	    fp_ksvm_true_error_c01 = fopen(outputFilenameFinal, "w");
	    if(fp_ksvm_true_error_c01==NULL) {
	      printf("cannot open %s_ksvm_true_error.txt", outputFilename);
	      return(-1);
	    }


	    // training time of each classifier
#ifdef REC_TIME	    

	    sprintf(outputFilenameFinal, "%s_trn_lda_time.txt", outputFilename);
	    fp_trn_lda_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_lda_time==NULL)
	      {
		printf("cannot open %s_trn_lda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_dlda_time.txt", outputFilename);
	    fp_trn_dlda_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_dlda_time==NULL)
	      {
		printf("cannot open %s_trn_dlda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_edc_time.txt", outputFilename);
	    fp_trn_edc_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_edc_time==NULL)
	      {
		printf("cannot open %s_trn_edc_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_zar_time.txt", outputFilename);
	    fp_trn_zar_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_zar_time==NULL)
	      {
		printf("cannot open %s_trn_zar_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_rlda_time.txt", outputFilename);
	    fp_trn_rlda_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_rlda_time==NULL)
	      {
		printf("cannot open %s_trn_rlda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_g13_time.txt", outputFilename);
	    fp_trn_g13_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_g13_time==NULL)
	      {
		printf("cannot open %s_trn_g13_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_serd_time.txt", outputFilename);
	    fp_trn_serd_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_serd_time==NULL)
	      {
		printf("cannot open %s_trn_serd_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_serd_time_cov.txt", outputFilename);
	    fp_trn_serd_time_cov = fopen(outputFilenameFinal, "w");
	    if(fp_trn_serd_time_cov==NULL)
	      {
		printf("cannot open %s_trn_serd_time_cov.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_lsvm_time.txt", outputFilename);
	    fp_trn_lsvm_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_lsvm_time==NULL)
	      {
		printf("cannot open %s_trn_lsvm_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_trn_ksvm_time.txt", outputFilename);
	    fp_trn_ksvm_time = fopen(outputFilenameFinal, "w");
	    if(fp_trn_ksvm_time==NULL) {
	      printf("cannot open %s_trn_ksvm_time.txt", outputFilename);
	      return(-1);
	    }

	    ////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////
	    ////////////////////////////////////////////////////////////////
	    
	    sprintf(outputFilenameFinal, "%s_tst_lda_time.txt", outputFilename);
	    fp_tst_lda_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_lda_time==NULL)
	      {
		printf("cannot open %s_tst_lda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_dlda_time.txt", outputFilename);
	    fp_tst_dlda_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_dlda_time==NULL)
	      {
		printf("cannot open %s_tst_dlda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_edc_time.txt", outputFilename);
	    fp_tst_edc_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_edc_time==NULL)
	      {
		printf("cannot open %s_tst_edc_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_zar_time.txt", outputFilename);
	    fp_tst_zar_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_zar_time==NULL)
	      {
		printf("cannot open %s_tst_zar_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_rlda_time.txt", outputFilename);
	    fp_tst_rlda_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_rlda_time==NULL)
	      {
		printf("cannot open %s_tst_rlda_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_g13_time.txt", outputFilename);
	    fp_tst_g13_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_g13_time==NULL)
	      {
		printf("cannot open %s_tst_g13_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_serd_time.txt", outputFilename);
	    fp_tst_serd_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_serd_time==NULL)
	      {
		printf("cannot open %s_tst_serd_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_serd_time_cov.txt", outputFilename);
	    fp_tst_serd_time_cov = fopen(outputFilenameFinal, "w");
	    if(fp_tst_serd_time_cov==NULL)
	      {
		printf("cannot open %s_tst_serd_time_cov.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_lsvm_time.txt", outputFilename);
	    fp_tst_lsvm_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_lsvm_time==NULL)
	      {
		printf("cannot open %s_tst_lsvm_time.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_tst_ksvm_time.txt", outputFilename);
	    fp_tst_ksvm_time = fopen(outputFilenameFinal, "w");
	    if(fp_tst_ksvm_time==NULL) {
	      printf("cannot open %s_tst_ksvm_time.txt", outputFilename);
	      return(-1);
	    }

#endif	

#ifdef CALC_ROC

	    sprintf(outputFilenameFinal, "%s_roc_lda.txt", outputFilename);
	    roc_lda = fopen(outputFilenameFinal, "w");
	    if(roc_lda==NULL)
	      {
		printf("cannot open %s_roc_lda.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_dlda.txt", outputFilename);
	    roc_dlda = fopen(outputFilenameFinal, "w");
	    if(roc_dlda==NULL)
	      {
		printf("cannot open %s_roc_dlda.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_rlda.txt", outputFilename);
	    roc_rlda = fopen(outputFilenameFinal, "w");
	    if(roc_rlda==NULL)
	      {
		printf("cannot open %s_roc_rlda.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_edc.txt", outputFilename);
	    roc_edc = fopen(outputFilenameFinal, "w");
	    if(roc_edc==NULL)
	      {
		printf("cannot open %s_roc_edc.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_g13.txt", outputFilename);
	    roc_g13 = fopen(outputFilenameFinal, "w");
	    if(roc_g13==NULL)
	      {
		printf("cannot open %s_roc_g13.txt", outputFilename);
		return(-1);
	      }
	    
	    sprintf(outputFilenameFinal, "%s_roc_zar.txt", outputFilename);
	    roc_zar = fopen(outputFilenameFinal, "w");
	    if(roc_zar==NULL)
	      {
		printf("cannot open %s_roc_zar.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_serd.txt", outputFilename);
	    roc_serd = fopen(outputFilenameFinal, "w");
	    if(roc_serd==NULL)
	      {
		printf("cannot open %s_roc_serd.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_serd_cov.txt", outputFilename);
	    roc_serd_cov = fopen(outputFilenameFinal, "w");
	    if(roc_serd_cov==NULL)
	      {
		printf("cannot open %s_roc_serd_cov.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_lsvm.txt", outputFilename);
	    roc_lsvm = fopen(outputFilenameFinal, "w");
	    if(roc_lsvm==NULL)
	      {
		printf("cannot open %s_roc_lsvm.txt", outputFilename);
		return(-1);
	      }

	    sprintf(outputFilenameFinal, "%s_roc_ksvm.txt", outputFilename);
	    roc_ksvm = fopen(outputFilenameFinal, "w");
	    if(roc_ksvm==NULL)
	      {
		printf("cannot open %s_roc_ksvm.txt", outputFilename);
		return(-1);
	      }

#endif
	    

	    
	  }	//if (MYRANK==0)
	
	best_features = new int[d] ();
	best_features2 = new int[d] ();

	dataLogfy(X, N, D, Xlog);
	dataStandardize(Xlog, N, D, Xstan);
	  
	featureSelection(Xstan, y, N, D, d, best_features2);

	//PrintInt_a(best_features2, d, (char*)("best_features.txt"));
	printf("total feature size=%d\n",D);
	printf("code feature size=%d\n",d);

	printf("N_trn is %d\n", N_trn);
	printf("N_tst is %d\n", N_tst);
	
	for (rep=MYRANK; rep<max_rep; rep+=NP) {
	  
	  
	  printf("rep=%d\n",rep);
	  //printf("prior=%.2f\n",prior);
	  //printf("\n\n\n\n");
	  
	  SimulationData data_trn;
	  SimulationData data_tst;		
	  double dinit = 0.00;		
	  
	  data_trn.data = make_2D_matrix(N_trn, D, dinit);
	  data_trn.labels = new int [N_trn];
	  data_tst.data = make_2D_matrix(N_tst, D, dinit);
	  data_tst.labels = new int [N_tst];
	  
	  //printf("data_trn is %d\n", N_trn);
	  //printf("data_tst is %d\n", N_tst);
	  //printf("data_trn is %d\n", N_trn);
	  //printf("data_trn is %d\n", N_trn);

	  dataGeneration(Xstan, y, N_trn, N_tst, D, &seed, &data_trn, &data_tst);
	  
	  ijkl = 0;
	  
	  prior = 2;
	  
	  //printf("calculateErrors start\n");
	  
	  calculateErrors(data_trn,
			  data_tst, 
			  N_trn,
			  N_tst,
			  D,
			  d,
			  prior,
			  &seed,
			  &all_errors,
			  &all_trn_times,
			  &all_tst_times,
			  &all_auc,
			  best_features2,
			  serd_group,
			  serd_elingr,
			  serd_cov_group,
			  serd_cov_elingr);
	  //lda_true_error,
	  //qda_true_error,
	  //lsvm_true_error,
	  //ksvm_true_error,
	  //dlda_true_error,
	  //edc_true_error,
	  //g13_true_error,
	  //serd_true_error,
	  //zar_true_error,
	  //rlda_true_error);
	  
	  //lda_true_error = all_errors.lda_true_error;
	  //lsvm_true_error = all_errors.lsvm_true_error;
	  //ksvm_true_error = all_errors.ksvm_true_error;
	  //dlda_true_error = all_errors.dlda_true_error;
	  //edc_true_error = all_errors.edc_true_error;
	  //g13_true_error = all_errors.g13_true_error;
	  //serd_true_error = all_errors.serd_true_error;
	  //zar_true_error = all_errors.zar_true_error;
	  //rlda_true_error = all_errors.rlda_true_error;
	  //serd_true_error = all_errors.lda_true_error;
	  
	  //printf("2 point\n");
	  
	  //lda_cvkfold_error = all_errors.lda_cvkfold_error;
	  //qda_cvkfold_error = all_errors.qda_cvkfold_error;
	  //lsvm_cvkfold_error = all_errors.lsvm_cvkfold_error;
	  //ksvm_cvkfold_error = all_errors.ksvm_cvkfold_error;
	  //
	  ////	printf("3 point\n");
	  //	
	  //lda_cvkfold_Mixture_error = all_errors.lda_cvkfold_Mixture_error;
	  //qda_cvkfold_Mixture_error = all_errors.qda_cvkfold_Mixture_error;
	  //lsvm_cvkfold_Mixture_error = all_errors.lsvm_cvkfold_Mixture_error;
	  //ksvm_cvkfold_Mixture_error = all_errors.ksvm_cvkfold_Mixture_error;
	  
	  //	printf("4 point\n");
	  
	  //lda_true_error = all_errors.lda_true_error;
	  //qda_true_error = all_errors.qda_true_error;
	  //lsvm_true_error = all_errors.lsvm_true_error;
	  
	  
	  /*
	    lda_cvkfold_error = 0.0;
	    qda_cvkfold_error = 0.0;
	    lsvm_cvkfold_error = 0.0;
	    
	    lda_cvkfold_Mixture_error = 0.0;
	    qda_cvkfold_Mixture_error = 0.0;
	    lsvm_cvkfold_Mixture_error = 0.0;
	  */
	  
	  
	  
	  //	printf("%d=",MYRANK);	
	  
	  
	  //printf("\n\nstart writing to the files\n");
	  //printf("all.errors.lda_true_error is %.5f\n", all_errors.lda_true_error);
	  //printf("all.errors.dlda_true_error is %.5f\n", all_errors.dlda_true_error);
	  //printf("all.errors.edc_true_error is %.5f\n", all_errors.edc_true_error);
	  //printf("all.errors.g13_true_error is %.5f\n", all_errors.g13_true_error);
	  //printf("all.errors.serd_true_error is %.5f\n", all_errors.serd_true_error);
	  //printf("all.errors.zar_true_error is %.5f\n", all_errors.zar_true_error);
	  //printf("all.errors.rlda_true_error is %.5f\n", all_errors.rlda_true_error);
	  //printf("all.errors.lsvm_true_error is %.5f\n", all_errors.lsvm_true_error);
	  //printf("all.errors.ksvm_true_error is %.5f\n", all_errors.ksvm_true_error);
			
		if (MYRANK==0)	{
		  switch(ijkl) {
		  case 0:
		    //printf("I am here at ijkl = 0\n\n");

		    fprintf(fp_lda_true_error_c01, "%f\t", all_errors.lda_true_error);
		    fprintf(fp_dlda_true_error_c01, "%f\t", all_errors.dlda_true_error);
		    fprintf(fp_g13_true_error_c01, "%f\t", all_errors.g13_true_error);
		    fprintf(fp_lsvm_true_error_c01, "%f\t", all_errors.lsvm_true_error);
		    fprintf(fp_serd_true_error_c01, "%f\t", all_errors.serd_true_error);
		    fprintf(fp_serd_true_error_c01_cov, "%f\t", all_errors.serd_true_error_cov);
		    fprintf(fp_rlda_true_error_c01, "%f\t", all_errors.rlda_true_error);
		    fprintf(fp_ksvm_true_error_c01, "%f\t", all_errors.ksvm_true_error);
		    fprintf(fp_edc_true_error_c01, "%f\t", all_errors.edc_true_error);
		    fprintf(fp_zar_true_error_c01, "%f\t", all_errors.zar_true_error);

#ifdef REC_TIME

		    fprintf(fp_trn_lda_time, "%f\n", all_trn_times.lda_time);
		    fprintf(fp_trn_dlda_time, "%f\n", all_trn_times.dlda_time);
		    fprintf(fp_trn_g13_time, "%f\n", all_trn_times.g13_time);
		    fprintf(fp_trn_lsvm_time, "%f\n", all_trn_times.lsvm_time);
		    fprintf(fp_trn_serd_time, "%f\n", all_trn_times.serd_time);
		    fprintf(fp_trn_serd_time_cov, "%f\n", all_trn_times.serd_time_cov);
		    fprintf(fp_trn_rlda_time, "%f\n", all_trn_times.rlda_time);
		    fprintf(fp_trn_ksvm_time, "%f\n", all_trn_times.ksvm_time);
		    fprintf(fp_trn_edc_time, "%f\n", all_trn_times.edc_time);
		    fprintf(fp_trn_zar_time, "%f\n", all_trn_times.zar_time);

		    ////////////////////////////////////////////////////////////////
		    ////////////////////////////////////////////////////////////////
		    ////////////////////////////////////////////////////////////////
		    
		    fprintf(fp_tst_lda_time, "%f\n", all_tst_times.lda_time);
		    fprintf(fp_tst_dlda_time, "%f\n", all_tst_times.dlda_time);
		    fprintf(fp_tst_g13_time, "%f\n", all_tst_times.g13_time);
		    fprintf(fp_tst_lsvm_time, "%f\n", all_tst_times.lsvm_time);
		    fprintf(fp_tst_serd_time, "%f\n", all_tst_times.serd_time);
		    fprintf(fp_tst_serd_time_cov, "%f\n", all_tst_times.serd_time_cov);
		    fprintf(fp_tst_rlda_time, "%f\n", all_tst_times.rlda_time);
		    fprintf(fp_tst_ksvm_time, "%f\n", all_tst_times.ksvm_time);
		    fprintf(fp_tst_edc_time, "%f\n", all_tst_times.edc_time);
		    fprintf(fp_tst_zar_time, "%f\n", all_tst_times.zar_time);
		    
#endif


#ifdef CALC_ROC

		    fprintf(roc_lda      , "%f\n", all_auc.lda_auc      );
		    fprintf(roc_edc      , "%f\n", all_auc.edc_auc      );
		    fprintf(roc_dlda     , "%f\n", all_auc.dlda_auc     );
		    fprintf(roc_g13      , "%f\n", all_auc.g13_auc      );
		    fprintf(roc_serd     , "%f\n", all_auc.serd_auc     );
		    fprintf(roc_serd_cov , "%f\n", all_auc.serd_auc_cov );
		    fprintf(roc_zar      , "%f\n", all_auc.zar_auc      );
		    fprintf(roc_rlda     , "%f\n", all_auc.rlda_auc     );
		    fprintf(roc_lsvm     , "%f\n", all_auc.lsvm_auc     );
		    fprintf(roc_ksvm     , "%f\n", all_auc.ksvm_auc     );

#endif

		    
		    //for(i=0; i<5; i++) {
		      //	printf("5 point\n");
		    //fprintf(fp_lda_true_error_c01, "%f\t", lda_true_error[i]);
		    //fprintf(fp_edc_true_error_c01, "%f\t", edc_true_error[i]);
		    //fprintf(fp_dlda_true_error_c01, "%f\t", dlda_true_error[i]);
		    //fprintf(fp_g13_true_error_c01, "%f\t", g13_true_error[i]);
		    //fprintf(fp_serd_true_error_c01, "%f\t", serd_true_error[i]);
		    //fprintf(fp_zar_true_error_c01, "%f\t", zar_true_error[i]);
		    //fprintf(fp_rlda_true_error_c01, "%f\t", rlda_true_error[i]);
		    ////printf("6 point\n");
		    ////printf("%f\t", lda_true_error[i]);
		    ////printf("7 point\n");
		    ////fprintf(fp_qda_true_error_c01, "%f\t", qda_true_error[i]);
		    //fprintf(fp_lsvm_true_error_c01, "%f\t", lsvm_true_error[i]);
		    //fprintf(fp_ksvm_true_error_c01, "%f\t", ksvm_true_error[i]);
		      //}
		    //printf("8 point\n");
		    fprintf(fp_lda_true_error_c01, "\n");
		    fprintf(fp_edc_true_error_c01, "\n");
		    fprintf(fp_dlda_true_error_c01, "\n");
		    fprintf(fp_g13_true_error_c01, "\n");
		    fprintf(fp_serd_true_error_c01, "\n");
		    fprintf(fp_serd_true_error_c01_cov, "\n");
		    fprintf(fp_zar_true_error_c01, "\n");
		    fprintf(fp_rlda_true_error_c01, "\n");
		    //fprintf(fp_qda_true_error_c01, "\n");
		    fprintf(fp_lsvm_true_error_c01, "\n");
		    fprintf(fp_ksvm_true_error_c01, "\n");
		    
		    //fprintf(fp_lda_cvkfold_error_c01, "%f\n", lda_cvkfold_error);
		    //fprintf(fp_qda_cvkfold_error_c01, "%f\n", qda_cvkfold_error);
		    //fprintf(fp_lsvm_cvkfold_error_c01, "%f\n", lsvm_cvkfold_error);
                    //fprintf(fp_ksvm_cvkfold_error_c01, "%f\n", ksvm_cvkfold_error);
		    //	
		    //fprintf(fp_lda_cvkfold_Mixture_error_c01, "%f\n", lda_cvkfold_Mixture_error);
		    //fprintf(fp_qda_cvkfold_Mixture_error_c01, "%f\n", qda_cvkfold_Mixture_error);
		    //fprintf(fp_lsvm_cvkfold_Mixture_error_c01, "%f\n", lsvm_cvkfold_Mixture_error);
                    //fprintf(fp_ksvm_cvkfold_Mixture_error_c01, "%f\n", ksvm_cvkfold_Mixture_error);
		    break;
				default:
					printf("Unexpected prior: Error");
					break;
			}

#ifdef USE_MPI
			flag = ((max_rep-rep)>(max_rep)%NP ? NP : max_rep%NP);
			for (j=1; j<flag; j++)	{

			  MPI_Recv(&(all_errors.lda_true_error), 1, MPI_DOUBLE, j, 6, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.dlda_true_error), 1, MPI_DOUBLE, j, 7, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.edc_true_error), 1, MPI_DOUBLE, j, 8, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.serd_true_error), 1, MPI_DOUBLE, j, 9, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.zar_true_error), 1, MPI_DOUBLE, j, 10, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.g13_true_error), 1, MPI_DOUBLE, j, 11, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.rlda_true_error), 1, MPI_DOUBLE, j, 12, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.lsvm_true_error), 1, MPI_DOUBLE, j, 13, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.ksvm_true_error), 1, MPI_DOUBLE, j, 14, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_errors.serd_true_error_cov), 1, MPI_DOUBLE, j, 15, MPI_COMM_WORLD,&status);

#ifdef CALC_ROC

			  MPI_Recv(&(all_auc.lda_auc      ), 1, MPI_DOUBLE, j, 6, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.dlda_auc     ), 1, MPI_DOUBLE, j, 7, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.edc_auc      ), 1, MPI_DOUBLE, j, 8, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.serd_auc     ), 1, MPI_DOUBLE, j, 9, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.zar_auc      ), 1, MPI_DOUBLE, j, 10, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.g13_auc      ), 1, MPI_DOUBLE, j, 11, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.rlda_auc     ), 1, MPI_DOUBLE, j, 12, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.lsvm_auc     ), 1, MPI_DOUBLE, j, 13, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.ksvm_auc     ), 1, MPI_DOUBLE, j, 14, MPI_COMM_WORLD,&status);
			  MPI_Recv(&(all_auc.serd_auc_cov ), 1, MPI_DOUBLE, j, 15, MPI_COMM_WORLD,&status);


			  fprintf(roc_lda      , "%f\n", all_auc.lda_auc      );
			  fprintf(roc_edc      , "%f\n", all_auc.edc_auc      );
			  fprintf(roc_dlda     , "%f\n", all_auc.dlda_auc     );
			  fprintf(roc_g13      , "%f\n", all_auc.g13_auc      );
			  fprintf(roc_serd     , "%f\n", all_auc.serd_auc     );
			  fprintf(roc_serd_cov , "%f\n", all_auc.serd_auc_cov );
			  fprintf(roc_zar      , "%f\n", all_auc.zar_auc      );
			  fprintf(roc_rlda     , "%f\n", all_auc.rlda_auc     );
			  fprintf(roc_lsvm     , "%f\n", all_auc.lsvm_auc     );
			  fprintf(roc_ksvm     , "%f\n", all_auc.ksvm_auc     );

#endif


			  
			  //MPI_Recv(lsvm_true_error, 5, MPI_DOUBLE, j, 8, MPI_COMM_WORLD,&status);
			  //MPI_Recv(ksvm_true_error, 5, MPI_DOUBLE, j, 9, MPI_COMM_WORLD,&status);
			  //	
			  //MPI_Recv(&lda_cvkfold_error, 1, MPI_DOUBLE, j, 10, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&qda_cvkfold_error, 1, MPI_DOUBLE, j, 11, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&lsvm_cvkfold_error, 1, MPI_DOUBLE, j, 12, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&ksvm_cvkfold_error, 1, MPI_DOUBLE, j, 13, MPI_COMM_WORLD,&status);
			  //	
			  //MPI_Recv(&lda_cvkfold_Mixture_error, 1, MPI_DOUBLE, j, 14, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&qda_cvkfold_Mixture_error, 1, MPI_DOUBLE, j, 15, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&lsvm_cvkfold_Mixture_error, 1, MPI_DOUBLE, j, 16, MPI_COMM_WORLD,&status);
			  //MPI_Recv(&ksvm_cvkfold_Mixture_error, 1, MPI_DOUBLE, j, 17, MPI_COMM_WORLD,&status);
			  
			  fprintf(fp_lda_true_error_c01, "%f\t", all_errors.lda_true_error);
			  fprintf(fp_edc_true_error_c01, "%f\t", all_errors.edc_true_error);
			  fprintf(fp_dlda_true_error_c01, "%f\t", all_errors.dlda_true_error);
			  fprintf(fp_g13_true_error_c01, "%f\t", all_errors.g13_true_error);
			  fprintf(fp_serd_true_error_c01, "%f\t", all_errors.serd_true_error);
			  fprintf(fp_serd_true_error_c01_cov, "%f\t", all_errors.serd_true_error_cov);
			  fprintf(fp_zar_true_error_c01, "%f\t", all_errors.zar_true_error);
			  fprintf(fp_rlda_true_error_c01, "%f\t", all_errors.rlda_true_error);
			  fprintf(fp_lsvm_true_error_c01, "%f\t", all_errors.lsvm_true_error);
			  fprintf(fp_ksvm_true_error_c01, "%f\t", all_errors.ksvm_true_error);
			  //printf("6 point\n");
			  //printf("%f\t", lda_true_error[i]);
			  //printf("7 point\n");
			  //fprintf(fp_qda_true_error_c01, "%f\t", qda_true_error[i]);
			  //printf("8 point\n");
			  fprintf(fp_lda_true_error_c01, "\n");
			  fprintf(fp_edc_true_error_c01, "\n");
			  fprintf(fp_dlda_true_error_c01, "\n");
			  fprintf(fp_g13_true_error_c01, "\n");
			  fprintf(fp_serd_true_error_c01, "\n");
			  fprintf(fp_serd_true_error_c01_cov, "\n");
			  fprintf(fp_zar_true_error_c01, "\n");
			  fprintf(fp_rlda_true_error_c01, "\n");
			  fprintf(fp_lsvm_true_error_c01, "\n");
			  fprintf(fp_ksvm_true_error_c01, "\n");
			  
				}
#endif
				
			} //if (MYRANK==0)
				
		 
#ifdef USE_MPI
		else {
		  MPI_Send(&(all_errors.lda_true_error), 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.dlda_true_error), 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.edc_true_error), 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.serd_true_error), 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.zar_true_error), 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.g13_true_error), 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.rlda_true_error), 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.lsvm_true_error), 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.ksvm_true_error), 1, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD);
		  MPI_Send(&(all_errors.serd_true_error_cov), 1, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);

#ifdef CALC_ROC

		  printf("\n\ni am in calcroc that is inside use_mpi\n\n");

		  MPI_Send(&(all_auc.lda_auc      ), 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.dlda_auc     ), 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.edc_auc      ), 1, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.serd_auc     ), 1, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.zar_auc      ), 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.g13_auc      ), 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.rlda_auc     ), 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.lsvm_auc     ), 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.ksvm_auc     ), 1, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD);
		  MPI_Send(&(all_auc.serd_auc_cov ), 1, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);

#endif

		  //MPI_Send(qda_true_error, 5, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
		  //MPI_Send(&lda_cvkfold_error, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		  //MPI_Send(&qda_cvkfold_error, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
		  //MPI_Send(&lsvm_cvkfold_error, 1, MPI_DOUBLE, 0, 12, MPI_COMM_WORLD);
		  //MPI_Send(&ksvm_cvkfold_error, 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
		  //	
		  //MPI_Send(&lda_cvkfold_Mixture_error, 1, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD);
		  //MPI_Send(&qda_cvkfold_Mixture_error, 1, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);
		  //MPI_Send(&lsvm_cvkfold_Mixture_error, 1, MPI_DOUBLE, 0, 16, MPI_COMM_WORLD);
		  //MPI_Send(&ksvm_cvkfold_Mixture_error, 1, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD);

		}
#endif
			
		//} //for (ijkl =
	
	delete_2D_matrix(N_trn, D, data_trn.data);
	delete [] data_trn.labels;
	delete_2D_matrix(N_tst, D, data_tst.data);
	delete [] data_tst.labels;
		
	} //for (rep=
		if (MYRANK==0) 		{
		  
		  fclose(fp_lda_true_error_c01);
		  fclose(fp_edc_true_error_c01);
		  fclose(fp_dlda_true_error_c01);
		  fclose(fp_g13_true_error_c01);
		  fclose(fp_serd_true_error_c01);
		  fclose(fp_serd_true_error_c01_cov);
		  fclose(fp_zar_true_error_c01);
		  fclose(fp_rlda_true_error_c01);
		  //fclose(fp_qda_true_error_c01);
		  fclose(fp_lsvm_true_error_c01);
		  fclose(fp_ksvm_true_error_c01);

#ifdef REC_TIME


		  fclose(fp_trn_lda_time);
		  fclose(fp_trn_edc_time);
		  fclose(fp_trn_dlda_time);
		  fclose(fp_trn_g13_time);
		  fclose(fp_trn_serd_time);
		  fclose(fp_trn_serd_time_cov);
		  fclose(fp_trn_zar_time);
		  fclose(fp_trn_rlda_time);
		  fclose(fp_trn_lsvm_time);
		  fclose(fp_trn_ksvm_time);

		  fclose(fp_tst_lda_time);
		  fclose(fp_tst_edc_time);
		  fclose(fp_tst_dlda_time);
		  fclose(fp_tst_g13_time);
		  fclose(fp_tst_serd_time);
		  fclose(fp_tst_serd_time_cov);
		  fclose(fp_tst_zar_time);
		  fclose(fp_tst_rlda_time);
		  fclose(fp_tst_lsvm_time);
		  fclose(fp_tst_ksvm_time);

#endif		  

#ifdef CALC_ROC

		  fclose(roc_lda);
		  fclose(roc_dlda);
		  fclose(roc_rlda);

		  fclose(roc_g13);
		  fclose(roc_edc);
		  fclose(roc_zar);

		  fclose(roc_serd);
		  fclose(roc_serd_cov);

		  fclose(roc_lsvm);
		  fclose(roc_ksvm);

#endif

		}	


	delete_2D_matrix(N, D, X);
	delete_2D_matrix(N, D, Xlog);
	delete_2D_matrix(N, D, Xstan);
	delete [] y;
	delete [] best_features2;
	delete [] best_features;

#ifdef USE_MPI
	MPI_Finalize();
#endif

	printf("the whole program is finished\n\n");

	return 0;
}

