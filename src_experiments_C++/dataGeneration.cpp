/*
 * dataGeneration.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "dataGeneration.h"
#include "print.h"

void dataGeneration(double** X,
		int* y,
		int N_trn,
		int N_tst,
		int D,
		long* seed,
		SimulationData* data_trn,
		SimulationData* data_tst)
{
	double dinit = 0.00;
	int i, j, k;
	int N_trn_0; // training samples in class 0
	int N_trn_1; // training samples in class 1
	int N = N_trn+N_tst;
	int N_0 = 0;
	int N_1 = 0;

	int* I0;
	int* I1;

	(*data_trn).N = N_trn;;
	(*data_trn).D = D;

	(*data_tst).N = N_tst;;
	(*data_tst).D = D;

	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	N_trn_0 = (int)floor((float)(*data_trn).N*(float)N_0/(float)N + 0.5);
	N_trn_1 = (*data_trn).N - N_trn_0;

	//N_trn_0 = (int)(floor((double)(*data_trn).N * n01r)); 
	//N_trn_1 = (*data_trn).N - N_trn_0;

	//printf("N_trn is %d\n", N_trn);
	//printf("N_tst is %d\n", N_tst);
	//
	//printf("N_trn_0 is %d\n", N_trn_0);
	//printf("N_trn_1 is %d\n", N_trn_1);
	//
	//printf("seed is %ld\n", *seed);

	
	//printf("N_trn_0=%d",N_trn_0);

	I0 = new int [N_0];
	I1 = new int [N_1];

	for (i=0; i<N_0; i++)
		I0[i] = i;

	for (i=0; i<N_1; i++)
		I1[i] = i;

	shuffleArray(N_0, seed, I0);
	shuffleArray(N_1, seed, I1);

	//PrintInt_a(I0, N_0, (char*)("I0.txt"));
	//PrintInt_a(I1, N_1, (char*)("I1.txt"));

	for(i=0; i<(*data_trn).N; i++)
	{
		if (i<N_trn_0)
		{
			for(j=0; j<(*data_trn).D; j++)
				(*data_trn).data[i][j] = X[I0[i]][j];
			(*data_trn).labels[i] = y[I0[i]];
		}
		else
		{
			for(j=0; j<(*data_trn).D; j++)
				(*data_trn).data[i][j] = X[I1[i-N_trn_0]+N_0][j];
			(*data_trn).labels[i] = y[I1[i-N_trn_0]+N_0];
		}
	}

	for(i=0; i<(*data_tst).N; i++)
	{
		if (i<N_0-N_trn_0)
		{
			for(j=0; j<(*data_trn).D; j++)
				(*data_tst).data[i][j] = X[I0[i+N_trn_0]][j];
			(*data_tst).labels[i] = y[I0[i+N_trn_0]];
		}
		else
		{
			for(j=0; j<(*data_trn).D; j++)
				(*data_tst).data[i][j] = X[I1[i-(N_0-N_trn_0)+N_trn_1]+N_0][j];
			(*data_tst).labels[i] = y[I1[i-(N_0-N_trn_0)+N_trn_1]+N_0];
		}
	}

	delete [] I0;
	delete [] I1;

	return;
}



void dataLogfy(double** X, int N, int D, double** Xlog){

  bool anyZero = false;

  for ( int i = 0; i < N; i++ ){
    for( int j = 0; j < D; j++ ){
      anyZero = X[i][j] <= 0 ? true : false;
      if (anyZero == true)
	break;

    }
    if (anyZero == true)
      break;
  }

  if ( anyZero == false ){
    for ( int i =  0 ; i < N ; i++ ){
      for ( int j = 0 ; j < D ; j++ ){
	Xlog[i][j] = log(X[i][j]);
      }
    }
  } else if (anyZero == true ){
    printf("non-positive value detected in X, returnring X\n");

    for ( int i = 0; i < N; i++ )
      for( int j = 0; j < D; j++ )
	Xlog[i][j] = X[i][j];

  }

  return;
 
}

/*
void dataStandardize(double** X, int* y, int N, int D, double** Xnew){

  double* fvec = new double[N]();
  double* fvec0, fvec1;
  double fmean0 = 0.00, fmean1 = 0.00;
  double var0 = 0.00, var1 = 0.00;
  int i_class0 = 0, i_class1 = 0;
  int N_class0 = 0, N_class1 = 0;

  for (int i=0; i<N; i++){
    if (y[i]==0) {
      N_class0++;
    } else if(y[i] == 1) {
      N_class1++;
    }
  }

  fvec0 = new double[N_class0]();
  fvec1 = new double[N_class1]();

  for ( int j = 0 ; j < D ; j++ ){
    i_class0 = 0;
    i_class1 = 0;
    for ( int i = 0 ; i < N ; i++ ){

      if (y[i]==0)
	{
	  vfec0[i_class0] = X[i][j];
	  i_class0++;
	}
      else
	{
	  vfec1[i_class1] = X[i][j];
	  i_class1++;
	}
      
    }
    mean( fvec0, 0, N_class0, &fmean0);
    mean( fvec1, 0, N_class1, &fmean1);
  }

  for (int i = 0 ; i < N_class0 ; i++ )
    var0 += (fvec0[i]-fmean0)*(fvec0[i]-fmean);
  
  for (int i = 0 ; i < N_class1 ; i++ )
    var1 += (fvec1[i]-fmean1)*(fvec1[i]-fmean);

  if (N_class0<2)
    {
      var1/=(N_class1-1);
      var0=var1;
    }
  else if (N_class1<2)
    {
      var0/=(N_class0-1);
      var1=var0;
    }
  else
    {
      var1/=(N_class1-1);
      var0/=(N_class0-1);
    }

}*/

void dataStandardize(double** X, int N, int D, double** Xnew){

  double fmean, fstd;

  for (int i = 0; i < D ; i++ ){

    double* fvec = new double[N]();
    fmean = 0; fstd = 0;

    // retrive fvec
    for ( int j = 0; j < N ; j++ ){
      fvec[j] = X[j][i];
    }

    // calculate mean and std
    mean(fvec, N, N, &fmean);

    for ( int j = 0 ; j < N ; j++ )
      fstd += (fvec[j] - fmean)*(fvec[j] - fmean);
    
    fstd /= (N-1);
    fstd = sqrt(fstd);
    //printf("gene is %d\tmean is %.5f\tstd is %.5f\n",i,fmean,fstd);

    // standardize and write
    for ( int j = 0 ; j < N ; j++ )
      Xnew[j][i] = (X[j][i] - fmean)/fstd;

    delete [] fvec;

  }

}
