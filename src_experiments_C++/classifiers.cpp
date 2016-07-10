/*
 * classifiers.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "classifiers.h"
#include "print.h"
#include "random.h"

int findIndOfMinElement ( double* inputArray, int lastElement ) {

  int i;
  
  int tempInd = 0;
  double tempValue = inputArray[0];
  
  for ( i = 0 ; i < lastElement; i++ ) {
    //printf("tempInd is %d\ttempValue is %f\tinputArray[i] is %f\n", tempInd, tempValue, inputArray[i]);
    if ( tempValue > inputArray[i] ) {
      tempValue = inputArray[i];
      tempInd = i;
    }
  }
  //printf("tempInd is %d\n", tempInd);
  return tempInd;
  
}

void estimateOptimumGamma(double* x0, double* x1, double** pooled_cov, regparam* gammas, int d, int n0, int n1 ) {

  double dinit = 0.00;
  double d_hat = dinit;
  double myGamma = dinit;
  double trH = dinit;
  double HCH = dinit;
  double myG = dinit;
  double mynom = dinit;
  double dasym = dinit;
  double alpha0 = dinit, alpha1 = dinit;

  double* Hdiff = new double [d] ();
  double* CHdiff = new double [d] ();
  double* x_diff = new double [d] ();
  double* error0 = new double [gammas->amount] ();
  double* error1 = new double [gammas->amount] ();
  double* dasym_error = new double [gammas->amount] ();
  double* all_dasym0 = new double [gammas->amount] ();
  double* all_dasym1 = new double [gammas->amount] ();
  double** iden = make_2D_identity_matrix(d, d, dinit, 1.00);

  double** tempH = make_2D_matrix(d, d, dinit);
  double** H = make_2D_matrix( d, d, dinit);

  int i, j;
  int x_class;
  int optGammaInd = 0;

  //PrintDouble_A( pooled_cov, d, d, (char*)"pooled_cov.txt");
  //PrintDouble_a( x0, d, (char*)"mu0.txt");
  //PrintDouble_a( x1, d, (char*)"mu1.txt");

  sum_two_vector ( x0 , x1 , d, 1.00, -1.00, 1.00 , x_diff );

  for ( i = 0 ; i < gammas->amount; i++ ){

    myGamma = gammas->gamma[i];

    sum_two_matrix(iden, pooled_cov, d, d, 1.00, myGamma, 1.00, tempH);
    invert_matrix(tempH, d, H);

    // calculate gamma*delta denoted as gd_hat
    trH = dinit;
    for ( j = 0 ; j < d ; j++ ){
      trH += H[j][j];
    }
    d_hat = ( (double) d - trH ) / (myGamma * ( (double)( n0 + n1 - 2 - d ) + trH ) );

    //printf("gamma is %.5f\n", myGamma);
    //printf("trH is %.5f\n", trH);
    //printf("d_hat is %.5f\n", d_hat);

    // calculate function D denoted as HCH
    multiply_A_by_b (H, x_diff, d, d, Hdiff);
    //multiply_a_by_B (x_diff, H, d, d, diffH); I think it is the same as Hdiff since t(H) = H
    multiply_A_by_b ( pooled_cov, Hdiff, d, d, CHdiff );
    multiply_a_by_b ( Hdiff , CHdiff, d, &HCH );

    // calculate G denoted as sHd
    //multiply_a_by_b ( x_diff, Hdiff , d, &myG );
    //myG = myG / 2;

    //printf("D is %.5f\tG is %.5f\n", HCH, myG);

    // calculate denominator and nominator
    HCH = sqrt( ( (double) 1.00 + myGamma*d_hat)*( (double) 1.00+myGamma*d_hat)*HCH );

    for ( x_class = 0 ; x_class < 2 ; x_class++ ){

      multiply_a_by_b ( x_diff, Hdiff, d, &myG ); // can be moved out from x_class loop, but myG would require a slight change
      
      if ( x_class == 0 ) {
	myG = myG/2;
	//mynom = -myG + ((double)(n0+n1-2))*d_hat/((double)(n0));
	mynom = -myG + ( (double)(n0+n1-2) ) *d_hat/( (double)(n0) ) + log( ((double)(n1)) / ((double)(n0)) );
	//mynom = -myG + (n0+n1-2)*d_hat/n0;
      }
      else if ( x_class == 1 ) {
	myG = -myG/2;
	//mynom = myG + ((double)(n0+n1-2))*d_hat/((double)(n1));
	mynom = myG + ((double)(n0+n1-2))*d_hat/((double)(n1)) - log(((double)(n1))/((double)(n0)));
								      //mynom = myG + (n0+n1-2)*d_hat/(n1) - log(n1/n0);
	//mynom = myG + (n0+n1-2)*d_hat/n1;
      }
      else
	printf("error is in estimateOptimumGamma, find nom\n\n");


      //printf("n0 is %d\tn1 is %d\t", n0, n1);
      //printf("D is %.5f\tG is %.5f\n", HCH, myG);
      //printf("x_class is %d\tdenom is %.5f\tnom is %.5f\tnom/denom is %.10f\n", x_class, HCH, mynom, mynom/HCH);
      // 
      dasym = mynom / HCH;
      //printf("dasym is %.5f\n", dasym);
      
      if ( x_class == 0 ) {
	all_dasym0[ i ] = dasym;
	error0[i] = ncdf( dasym );
      }
      else if ( x_class == 1 ) {
	all_dasym1[ i ] = dasym;
	error1[i] = ncdf( dasym );
      }
      else
	printf("error is in estimateOptimumGamma, find error0 or error1 \n\n");

    }
    //printf("\n\n");


  }

  //printf("i reached this position\n\n\n");

  // ask question from professor regarding alpha0 and alpha1, also about the prior bias or (-1)*c/k
  //alpha0 = ( double) n0 / (double)(n0+n1) );    alpha1 = (double) n1 / (double)(n0+n1) );
  sum_two_vector( error0, error1, gammas->amount, (double) n0, (double) n1, (double) (n0+n1), dasym_error );
  //sum_two_matrix( &error0, &error1, 1.00, d, 1.00, 1.00, 1.00, &dasym_error );

  //printf("error0 is \n");
  //for ( int iii = 0 ; iii < gammas->amount ; iii++ )
  //printf("%.5f\t", error0[iii]);
  //printf("\n");
  // 
  //printf("error1 is \n");
  //for ( int iii = 0 ; iii < gammas->amount ; iii++ )
  //printf("%.5f\t", error1[iii]);
  //printf("\n");
  // 
  //printf("dasym error is \n");
  //for ( int iii = 0 ; iii < gammas->amount ; iii++ )
  //printf("%.5f\t", dasym_error[iii]);
  //printf("\n");
  //
  //
  //printf("all_dasym0 error is \n");
  //for ( int iii = 0 ; iii < gammas->amount ; iii++ )
  //printf("%.5f\t", all_dasym0[iii]);
  //printf("\n");
  // 
  //printf("all_dasym1 error is \n");
  //for ( int iii = 0 ; iii < gammas->amount ; iii++ )
  //printf("%.5f\t", all_dasym1[iii]);
  //printf("\n");

  //PrintDouble_a( dasym_error, gammas->amount, (char*)"gamma_error.txt");

  optGammaInd = findIndOfMinElement ( dasym_error, gammas->amount );
  gammas->optGamma = gammas->gamma[optGammaInd];

  //printf("index of optimum gamma is %d\toptimum gamma is %.5f\n\n", optGammaInd, gammas->optGamma);

  delete [] all_dasym0;
  delete [] all_dasym1;


  delete [] Hdiff;
  delete [] CHdiff;
  delete [] x_diff;
  delete [] error0;
  delete [] error1;
  delete [] dasym_error;

  delete_2D_matrix(d, d, iden);
  delete_2D_matrix (d, d, tempH);
  delete_2D_matrix (d, d, H);

  return;
  
}

void rldaTrn(double** X, int* y, int N, int d, int* ind, double prior, regparam* gammas, RLDA* rlda) {
  
  int i, j;
  int N_0 = 0;
  int N_1 = 0;
  int i_0 = 0;
  int i_1 = 0;
  double dinit = 0.00;
  double temp_b = 0.00;
  double** X_0;
  double** X_1;
  double* mu_hat_0 = new double [d];
  double* mu_hat_1 = new double [d];
  double* temp = new double [d];
  double* temp2 = new double [d];
  double** cov_hat_0 = make_2D_matrix(d, d, dinit);
  double** cov_hat_1 = make_2D_matrix(d, d, dinit);
  double** pooled_cov = make_2D_matrix(d, d, dinit);
  double** tempH = make_2D_matrix(d, d, dinit);
  double** inv_pooled_cov = make_2D_matrix(d, d, dinit);

  //printf("tst_N is %d\trlda_optGamma_trnN is %d\n", tst_N, rlda_optGamma_trnN);

  for (i=0; i<N; i++)
    if (y[i]==0)
      N_0++;
    else
      N_1++;
  
  X_0 = make_2D_matrix(N_0, d, dinit);
  X_1 = make_2D_matrix(N_1, d, dinit);
  
  for (i=0; i<N; i++)
    {
      if (y[i]==0)
	{
	  for (j=0; j<d; j++)
	    X_0[i_0][j] = X[i][ind[j]];
	  i_0++;
	}
      else
	{
	  for (j=0; j<d; j++)
	    X_1[i_1][j] = X[i][ind[j]];
	  i_1++;
	}
    }

  mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
  mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector
  cov(X_0, N_0, d, cov_hat_0);
  cov(X_1, N_1, d, cov_hat_1);
  
  sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
  sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp2); // hat{X1} + hat{X0}
  sum_two_matrix(cov_hat_0, cov_hat_1, d, d, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);

  estimateOptimumGamma( mu_hat_0, mu_hat_1, pooled_cov, gammas, d, N_0, N_1 );

  inv_pooled_cov = make_2D_identity_matrix( d, d, 0.00, 1.00 );

  //printf("optimum gamma is %.5f\n", gammas->optGamma);
  sum_two_matrix( pooled_cov, inv_pooled_cov, d, d, gammas->optGamma, 1.00, 1.00, tempH );
  invert_matrix( tempH, d, inv_pooled_cov);

  multiply_A_by_b (inv_pooled_cov, temp, d, d, (*rlda).a); // multiply inv{cov} by // hat{X1} - hat{X0}
  
  sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp); // hat{X1} + hat{X0}
  multiply_a_by_b ((*rlda).a, temp2, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b
  
  if (prior==2.00){
    //printf("We are here 1");
    (*rlda).b = -0.50*temp_b + log((double)N_1/(double)N_0);
  } else {
    //printf("We are here 2");
    (*rlda).b = -0.50*temp_b + log(prior/(1.00-prior));
  }
  
  delete_2D_matrix(N_0, d, X_0);
  delete_2D_matrix(N_1, d, X_1);
  delete [] mu_hat_0;
  delete [] mu_hat_1;
  delete [] temp;
  delete_2D_matrix(d, d, cov_hat_0);
  delete_2D_matrix(d, d, cov_hat_1);
  delete_2D_matrix(d, d, pooled_cov);
  delete_2D_matrix(d, d, inv_pooled_cov);
  delete_2D_matrix(d, d, tempH);

  return;
  
}

double rldaTst(double** X, int* y, int N, int d, int* ind, double prior, RLDA rlda)
{

	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, rlda.a, d, &temp);
		y_hat = ((temp+rlda.b) <= 0) ? 0 : 1; 
		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}
    
	if (prior==2.00)
	{
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
		//printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("rlda error is %.3f\n", error);

	delete temp_X;

	return error;
}

void zarTrn( double** X, int* y, int N, int d, int* ind, double prior,
	     ZAR* zar ){
  
  int i, j;
  int N_0 = 0;
  int N_1 = 0;
  int i_0 = 0;
  int i_1 = 0;
  double dinit = 0.00;
  double temp_b = 0.00;
  double** X_0;
  double** X_1;
  double* mu_hat_0 = new double [d]();
  double* mu_hat_1 = new double [d]();
  double* temp = new double [d];
  double** cov_hat_0 = make_2D_matrix(d, d, dinit);
  double** cov_hat_1 = make_2D_matrix(d, d, dinit);
  double** pooled_cov = make_2D_matrix(d, d, dinit);

  double temp_cov = 0.00;

  double** corMat = make_2D_matrix( d, d, 0);
  double** ZarCovMat2 = make_2D_matrix( d, d, 0);

  double** ZarCovMat = make_2D_matrix( d, d, 0);
  double** ZarCovMatT = make_2D_matrix( d, d, 0);    // transpose of ZarCovMat
  
  int temp_ind = 0;   // to identify the argmax of a row
  int* zar_ind = new int[d-1]();   // write argmax here
  double temp_value = dinit;    // needed to identify argmax

  // calculate amount of class0 and class1
  for (i=0; i<N; i++)
    if (y[i]==0)
      N_0++;
    else
      N_1++;
  
  // create X0 and X1 and write values to them from X
  X_0 = make_2D_matrix(N_0, d, dinit);
  X_1 = make_2D_matrix(N_1, d, dinit);
  
  for (i=0; i<N; i++)
    {
      if (y[i]==0)
	{
	  for (j=0; j<d; j++)
	    X_0[i_0][j] = X[i][ind[j]];
	  i_0++;
	}
      else
	{
	  for (j=0; j<d; j++)
	    X_1[i_1][j] = X[i][ind[j]];
	  i_1++;
	}
    }
  
  mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
  mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector
  cov(X_0, N_0, d, cov_hat_0);
  cov(X_1, N_1, d, cov_hat_1);

  sum_two_matrix(cov_hat_0, cov_hat_1, d, d, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);


  for ( i = 0 ; i < d; i++ ){
    for ( j = 0 ; j < d ; j++ ){
      if ( i == j ){  // if diagonal calculate
	corMat[i][j] = 1;
      } else if ( j > i ) { // if higher triangle new calculation
	corMat[i][j] = pooled_cov[i][j] / sqrt( pooled_cov[i][i]*pooled_cov[j][j] );
      } else if ( j < i ) { // if lower triangle copy
	corMat[i][j] = corMat[j][i];
      }
    }
  }

  // running the total Kruskal Algorithm, which is equievalent to selecting the maximum values of the lower triangle in corMat
  // the easiest way to justify this is by drawing a graph, selecting only one element that is maximum for each branch will never make a loop, and maximum values for each branch can be selected
  Edge *result = new Edge[d-1];
  //maxSpanningTreeMatrix(d,corMat,covMatKruskal); // aboslute value is taken within the maxSpanningTreeMatrix method
  maxSpanningTreeMatrix(d,corMat,result); // aboslute value is taken within the maxSpanningTreeMatrix method

  //printf("Edge is\n");
  //for (int ii = 0 ; ii < d-1; ii++ ){
  //printf( "%d\t%d\n", result[ii].row, result[ii].col );
  //}
  //printf("\n");
  
  /*
    // next lines can replace maxSpanningTreeMatrix generation by taking the index of the absolute maximum element within a row
  for ( i = 1; i < d ; i++ ){
    temp_value = 0;
    for ( j = 0; j < i; j++ ){
      if ( fabs(corMat[i][j]) > temp_value ){
	temp_value = fabs(corMat[i][j]);
	temp_ind = j;
      }
    }
    zar_ind[i-1] = temp_ind;
  }
  */

  // comment these lines if you want to use the simple method via finding the max element within the row
  for ( int ii = 0 ; ii < d - 1; ii++ )
    zar_ind[result[ii].row-1] = result[ii].col;
  delete[] result;

  ZarCovMat2[0][0] = 1;

  for ( i = 0; i < d-1; i++ ){
    temp_cov = pooled_cov[i+1][zar_ind[i]] / sqrt( pooled_cov[i+1][i+1] * pooled_cov[zar_ind[i]][zar_ind[i]] );
    ZarCovMat2[i+1][zar_ind[i]] = - temp_cov / sqrt( pooled_cov[zar_ind[i]][zar_ind[i]]*(1-temp_cov*temp_cov) );
    ZarCovMat2[i+1][i+1] = 1 / sqrt( pooled_cov[i+1][i+1] * ( 1 - temp_cov * temp_cov ) );
  }

  transpose_2D_matrix(ZarCovMat2,d,ZarCovMatT);
  multiply_A_by_B(ZarCovMatT, ZarCovMat2, d, d, d, ZarCovMat );  // final CovMat
  
  //PrintDouble_A(ZarCovMat2,d,d, (char*)"cormat.txt");     // correlation matrix
  //PrintDouble_A(ZarCovMat2,d,d, (char*)"zar_lower.txt"); // lower triange
  //PrintDouble_A(ZarCovMat,d,d, (char*)"zar_total.txt");  // total matrix


  sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
  multiply_A_by_b (ZarCovMat, temp, d, d, (*zar).a); 
  sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp); // hat{X1} + hat{X0}
  multiply_a_by_b ((*zar).a, temp, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b

  if (prior==2.00) {
    //printf("We are here 1");
    (*zar).b = -0.50*temp_b + log((double)N_1/(double)N_0);
  } else {
    //printf("We are here 2");
    (*zar).b = -0.50*temp_b + log(prior/(1.00-prior));
  }
  
  //printf("zar.b is %.3f\n", zar->b);
  // 
  //printf("zar.a is\n");
  //for ( i = 0 ; i < d ; i++ )
  //printf("%.3f\t", zar->a[i]);
  //printf("\n");
  
  delete [] mu_hat_0;
  delete [] mu_hat_1;
  delete [] temp;
  delete [] zar_ind;


  delete_2D_matrix( d, d, cov_hat_0 );
  delete_2D_matrix( d, d, cov_hat_1 );
  delete_2D_matrix( d, d, pooled_cov );

  delete_2D_matrix( d, d, corMat );
  delete_2D_matrix( d, d, ZarCovMat2 );
  delete_2D_matrix( d, d, ZarCovMat );
  delete_2D_matrix( d, d, ZarCovMatT );

  delete_2D_matrix( N_0, d, X_0 );
  delete_2D_matrix( N_1, d, X_1 );

  return;

}

double zarTst(double** X, int* y, int N, int d, int* ind, double prior, ZAR zar)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;

	// uncomment next line if you want to compare with dlda.R
	//int* values = new int[N];
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, zar.a, d, &temp);
		y_hat = ((temp+zar.b) <= 0) ? 0 : 1; 

		// uncomment next line if you want to comapre with dlda.R
		//values[i] = y_hat;

		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}
    
	if (prior==2.00)
	{
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
		//printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("zar error is %.3f\n", error);

	//uncomment this section you want to compare with zar.R
	//PrintDouble( X, N, d, (char*)"terminal/ZAR_tstX.txt");
	//PrintInt( y, N, (char*)"terminal/ZAR_tsty.txt");
	//PrintInt( values, N , (char*)"terminal/ZAR_expy.txt");
	//printf("error is %f\n", error);
	//delete values;

	delete temp_X;

	return error;

}

void serdTrn(double** X, int* y, int N, int d, int* ind, double prior, SERD* serd, int k, int m)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d]();
	double* mu_hat_1 = new double [d]();
	double* temp = new double [k]();

	//int k = 5;
	//const int m = floor( (double)(d/k));

	double* X_hat_0 = new double [k]();
	double* X_hat_1 = new double [k]();

	double* X_hat = new double [k]();

	// calculate amount of class0 and class1
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	// create X0 and X1 and write values to them from X
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	  {
	    if (y[i]==0)
	      {
		for (j=0; j<d; j++)
		  X_0[i_0][j] = X[i][ind[j]];
		i_0++;
	      }
	    else
	      {
		for (j=0; j<d; j++)
		  X_1[i_1][j] = X[i][ind[j]];
		i_1++;
	      }
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_0 is a d dimension vector

	int kk = 0, mm = 0;

	//printf("mu0 is\n");
	//for (int kkk = 0 ; kkk < d; kkk++){
	//printf("%.5f\n",mu_hat_0[kkk]);
	//}

	divide_a_to_k_groups ( mu_hat_0 , d, k, m, X_hat_0 );
	divide_a_to_k_groups ( mu_hat_1 , d, k, m, X_hat_1 );

	//printf("X_hat_0 is\n");
	//for (int kkk = 0 ; kkk < k; kkk++){
	//printf("%.5f\n",X_hat_0[kkk]);
	//}
	//printf("\n");
	//
	//printf("X_hat_1 is\n");
	//for (int kkk = 0 ; kkk < k; kkk++){
	//printf("%.5f\n",X_hat_1[kkk]);
	//}
	//printf("\n");

	sum_two_vector( X_hat_0 , X_hat_1 , k, 1.00, 1.00 , 1.00 , X_hat );
	sum_two_vector( X_hat_0 , X_hat_1 , k, 1.00, -1.00 , 1.00 , temp );

	//printf("sum is is\n");
	//for (int kkk = 0 ; kkk < k; kkk++){
	//printf("%.5f\n",X_hat[kkk]);
	//}
	//printf("\n");
	//
	//printf("diff is is\n");
	//for (int kkk = 0 ; kkk < k; kkk++){
	//printf("%.5f\n",temp[kkk]);
	//}
	//printf("\n");
	

	multiply_a_by_scalar( temp , (double) m, k, (*serd).a);

	multiply_a_by_b( (*serd).a , X_hat, k, &temp_b);

	//printf("serd.a is is\n");
	//for (int kkk = 0 ; kkk < k; kkk++){
	//printf("%.5f\n",(*serd).a[kkk]);
	//}
	//printf("\n");
	//
	//printf("temp_b is %.5f\n", temp_b);

	if (prior==2.00){
	  (*serd).b = -0.50*temp_b;
	} else {
	  (*serd).b = -0.50*temp_b;
	}

	//printf("serd.b is %.5f\n", (*serd).b);
	//printf("\n\n");

	//if (prior==2.00){
	//(*serd).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	//} else {
	//(*serd).b = -0.50*temp_b + log(prior/(1.00-prior));
	//}

	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete [] mu_hat_0;
	delete [] mu_hat_1;
	delete [] temp;
	delete [] X_hat_0;
	delete [] X_hat_1;
	delete [] X_hat;

	return;
}


void serdCovTrn(double** X, int* y, int N, int d, int* ind, double prior, SERD* serd, int k, int m)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double** X_0;
	double** X_1;
	double** cov_hat_0;
	double** cov_hat_1;
	double** pooled_cov;
	double** inv_pooled_cov;
	double** serd_X_0;
	double** serd_X_1;
	double* mu_hat_0 = new double [d]();
	double* mu_hat_1 = new double [d]();
	double* temp = new double [k]();
	double* temp2 = new double [k]();

	//int k = 5;
	//const int m = floor( (double)(d/k));

	double* X_hat_0 = new double [k]();
	double* X_hat_1 = new double [k]();

	double* X_hat = new double [k]();

	// calculate amount of class0 and class1
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	// create X0 and X1 and write values to them from X
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	serd_X_0 = make_2D_matrix(N_0, k, 1.00);
	serd_X_1 = make_2D_matrix(N_1, k, 1.00);

	cov_hat_0 = make_2D_matrix(k,k, dinit);
	cov_hat_1 = make_2D_matrix(k,k, dinit);
	pooled_cov =  make_2D_matrix(k,k, dinit);
	inv_pooled_cov =  make_2D_matrix(k,k, dinit);

	for (i=0; i<N; i++)
	  {
	    if (y[i]==0)
	      {
		for (j=0; j<d; j++)
		  X_0[i_0][j] = X[i][ind[j]];
		i_0++;
	      }
	    else
	      {
		for (j=0; j<d; j++)
		  X_1[i_1][j] = X[i][ind[j]];
		i_1++;
	      }
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector

	divide_a_to_k_groups ( mu_hat_0 , d, k, m, X_hat_0 );
	divide_a_to_k_groups ( mu_hat_1 , d, k, m, X_hat_1 );

	divide_A_to_k_groups ( X_0 , N_0, d, k, m, serd_X_0 );  // X_0 is N_0 by d, serd_X_0 is N_0 by k
	divide_A_to_k_groups ( X_1 , N_1, d, k, m, serd_X_1 );

	cov(serd_X_0, N_0, k, cov_hat_0);
	cov(serd_X_1, N_1, k, cov_hat_1);

	//sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, k, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
	sum_two_matrix(cov_hat_0, cov_hat_1, k, k, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);
	diag_matrix(pooled_cov, k, cov_hat_0);
	invert_matrix(cov_hat_0, k, inv_pooled_cov);

	sum_two_vector( X_hat_0 , X_hat_1 , k, 1.00, 1.00 , 1.00 , X_hat );
	sum_two_vector( X_hat_0 , X_hat_1 , k, 1.00, -1.00 , 1.00 , temp );

	multiply_A_by_b( inv_pooled_cov , temp, k, k, temp2);

	multiply_a_by_scalar( temp2 , (double) m, k, (*serd).a);

	multiply_a_by_b( (*serd).a , X_hat, k, &temp_b);

	if (prior==2.00){
	  (*serd).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	} else {
	  (*serd).b = -0.50*temp_b + log(prior/(1.00-prior));
	}

	delete_2D_matrix(N_0,k,serd_X_0);
	delete_2D_matrix(N_1,k,serd_X_1);

	delete_2D_matrix(k,k,cov_hat_0);
	delete_2D_matrix(k,k,cov_hat_1);
	delete_2D_matrix(k,k,pooled_cov);
	delete_2D_matrix(k,k,inv_pooled_cov);

	delete [] mu_hat_0;
	delete [] mu_hat_1;
	delete [] X_hat_0;
	delete [] X_hat_1;
	delete [] X_hat;

	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);

	delete [] temp;
	delete [] temp2;


	return;
}

double serdTst(double** X, int* y, int N, int d, int* ind, double prior, SERD serd, int k, int m)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat, kk, mm;
	int N_0, N_1;

	//double*  = new double[m]();
	double* X_k = new double[k]();
	
	// uncomment next lines if you want to compare with serd.R
	//int* expValues = new int[N];
	//double* calValues = new double[N];

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		divide_a_to_k_groups ( temp_X , d, k , m , X_k );

		multiply_a_by_b (X_k, serd.a, k, &temp);
		y_hat = ((temp+serd.b) >= 0) ? 0 : 1; 

		// uncomment next lines if you want to comapre with serd.R
		//expValues[i] = y_hat;
		//calValues[i] = temp+serd.b;

		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}

	if (prior==2.00)
	  {
	    error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
	    //printf("ohoh prior = %.2f",prior);
	    //printf("ldaerror=%.2f\n",error);
	  } else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("serd error is %.3f\n", error);

	//uncomment this section you want to compare with serd.R
	//PrintDouble( X, N, d, (char*)"terminal/SERD_tstX.txt");
	//PrintInt( y, N, (char*)"terminal/SERD_tsty.txt");
	//PrintInt( expValues, N , (char*)"terminal/SERD_expy.txt");
	//PrintDouble2( calValues, N , (char*)"terminal/SERD_calcValues.txt");
	//delete expValues;
	//delete calValues;

	delete temp_X;

	return error;
}

void dldaTrn(double** X, int* y, int N, int d, int* ind, double prior, DLDA* dlda)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d];
	double* mu_hat_1 = new double [d];
	double* temp = new double [d];
	double** cov_hat_0 = make_2D_matrix(d, d, dinit);
	double** cov_hat_1 = make_2D_matrix(d, d, dinit);
	double** pooled_cov = make_2D_matrix(d, d, dinit);
	double** inv_pooled_cov = make_2D_matrix(d, d, dinit);
	double** diag_pooled = make_2D_matrix(d,d, dinit);

	// calculate amount of class0 and class1
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	// create X0 and X1 and write values to them from X
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			for (j=0; j<d; j++)
				X_0[i_0][j] = X[i][ind[j]];
			i_0++;
		}
		else
		{
			for (j=0; j<d; j++)
				X_1[i_1][j] = X[i][ind[j]];
			i_1++;
		}
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector
	cov(X_0, N_0, d, cov_hat_0);
	cov(X_1, N_1, d, cov_hat_1);

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
	sum_two_matrix(cov_hat_0, cov_hat_1, d, d, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);
	diag_matrix(pooled_cov, d, diag_pooled);
	invert_matrix(diag_pooled, d, inv_pooled_cov);
	multiply_A_by_b (inv_pooled_cov, temp, d, d, (*dlda).a); // multiply inv{cov} by // hat{X1} - hat{X0}

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp); // hat{X1} + hat{X0}
	multiply_a_by_b ((*dlda).a, temp, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b

	if (prior==2.00){
	  (*dlda).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	} else {
	  (*dlda).b = -0.50*temp_b + log(prior/(1.00-prior));
	}

	//printf("\norig matrix is\n");
	//for ( i = 0; i < d; i++){
	//for( j = 0; j < d; j++ ){
	//printf("%.7f\t",pooled_cov[i][j]);
	//}
	//printf("\n");
	//}
	//	
	//printf("\ndiag matrix is\n");
	//for ( i = 0; i < d; i++){
	//for( j = 0; j < d; j++ ){
	//printf("%.7f\t",diag_pooled[i][j]);
	//}
	//printf("\n");
	//}
	//// uncomment this section if you want to write the classificator to .txt files
	//PrintDouble( X , N_0 + N_1 , d, (char*)"terminal/DLDA_X.txt");
	//PrintDouble( X_0 , N_0 , d, (char*)"terminal/DLDA_X_0.txt");
	//PrintDouble( X_1 , N_1 , d, (char*)"terminal/DLDA_X_1.txt");
	//PrintDouble2( mu_hat_0 , d, (char*)"terminal/DLDA_mu0.txt");
	//PrintDouble2( mu_hat_1 , d, (char*)"terminal/DLDA_mu1.txt");
	//PrintInt( y, N_0 + N_1 , (char*)"terminal/DLDA_y.txt");
	//PrintEDC( ((*dlda).a), (*dlda).b , d , (char*)"terminal/DLDA.txt");
	//
	//printf("before clearing the space\n\n");

	//printf("(*dlda).b is %.3f\n", (*dlda).b);
	//	
	//printf("dlda.a is\n");
	//for ( i = 0 ; i < d; i++){
	//printf("%.3f\t", (*dlda).a[i]);
	//}printf("\n");
	
	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete [] mu_hat_0;
	delete [] mu_hat_1;
	delete [] temp;
	delete_2D_matrix(d, d, cov_hat_0);
	delete_2D_matrix(d, d, cov_hat_1);
	delete_2D_matrix(d, d, pooled_cov);
	delete_2D_matrix(d, d, inv_pooled_cov);
	delete_2D_matrix(d, d, diag_pooled);

	return;
}

double dldaTst(double** X, int* y, int N, int d, int* ind, double prior, DLDA dlda)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;
	
	// uncomment next line if you want to compare with dlda.R
	//int* values = new int[N];

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, dlda.a, d, &temp);
		y_hat = ((temp+dlda.b) <= 0) ? 0 : 1; 

		// uncomment next line if you want to comapre with dlda.R
		//values[i] = y_hat;

		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}

	if (prior==2.00) {
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
		//printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("dlda error is %.3f\n", error);

	//uncomment this section you want to compare with dlda.R
	//PrintDouble( X, N, d, (char*)"terminal/DLDA_tstX.txt");
	//PrintInt( y, N, (char*)"terminal/DLDA_tsty.txt");
	//PrintInt( values, N , (char*)"terminal/DLDA_expy.txt");
	//printf("error is %f\n", error);
	//delete values;

	delete temp_X;

	return error;
}

void g13Trn(double** X, int* y, int N, int d, int* ind, double prior, G13* g13)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d];
	double* mu_hat_1 = new double [d];
	double* temp = new double [d];
	double* temp2 = new double [d];
	double** cov_hat_0 = make_2D_matrix(d, d, dinit);
	double** cov_hat_1 = make_2D_matrix(d, d, dinit);
	double** pooled_cov = make_2D_matrix(d, d, dinit);
	double** inv_pooled_cov = make_2D_matrix(d, d, dinit);
	double* Xtemp = new double [d];
	double ccc;

	// calculate amount of class0 and class1
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	// create X0 and X1 and write values to them from X
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			for (j=0; j<d; j++)
				X_0[i_0][j] = X[i][ind[j]];
			i_0++;
		}
		else
		{
			for (j=0; j<d; j++)
				X_1[i_1][j] = X[i][ind[j]];
			i_1++;
		}
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector
	cov(X_0, N_0, d, cov_hat_0);
	cov(X_1, N_1, d, cov_hat_1);

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
	sum_two_matrix(cov_hat_0, cov_hat_1, d, d, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);
	invert_matrix(pooled_cov, d, inv_pooled_cov);
	multiply_A_by_b (inv_pooled_cov, temp, d, d, temp2 ); // multiply inv{cov} by // hat{X1} - hat{X0}

	ccc = ((double) ((double)(N_0 + N_1 - 2 - d)/(double)(N_0 + N_1 - 2)));

	multiply_a_by_scalar( temp2, ccc, d, (*g13).a );

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp2); // hat{X1} + hat{X0}
	multiply_a_by_b ( temp2, (*g13).a, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b

	
	if (prior==2.00){
	  (*g13).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	} else {
	  (*g13).b = -0.50*temp_b + log(prior/(1.00-prior));
	}
	

	// uncomment this section if you want to write the classificator to .txt files
	//PrintDouble( X , N_0 + N_1 , d, (char*)"terminal/G13_X.txt");
	//PrintDouble( X_0 , N_0 , d, (char*)"terminal/G13_X_0.txt");
	//PrintDouble( X_1 , N_1 , d, (char*)"terminal/G13_X_1.txt");
	//PrintDouble2( mu_hat_0 , d, (char*)"terminal/G13_mu0.txt");
	//PrintDouble2( mu_hat_1 , d, (char*)"terminal/G13_mu1.txt");
	//PrintInt( y, N_0 + N_1 , (char*)"terminal/G13_y.txt");
	//PrintEDC( ((*g13).a), (*g13).b , d , (char*)"terminal/G13.txt");
	//	
	   
	//printf("g13.b is %.3f\n", (*g13).b);
	//	
	//for ( i = 0 ; i < d; i++){
	//printf("%.3f\t", g13->a[i]);
	//}printf("\n");
	
	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete [] mu_hat_0;
	delete [] mu_hat_1;
	delete [] temp;
	delete [] Xtemp;
	delete [] temp2;
	delete_2D_matrix(d, d, cov_hat_0);
	delete_2D_matrix(d, d, cov_hat_1);
	delete_2D_matrix(d, d, pooled_cov);
	delete_2D_matrix(d, d, inv_pooled_cov);

	//printf("after space clearing\n");
	//printf("g13.b is %f\n\n", (*g13).b);
	//
	//printf("g13.a is\n");
	//for ( i = 0 ; i<d; i++ ){
	//printf("%f\n", (*g13).a[i]);
	//}

	return;
}


double g13Tst(double** X, int* y, int N, int d, int* ind, double prior, G13 g13)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;

	// uncomment next line if you want to compare with dlda.R
	//int* values = new int[N];
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];


	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, g13.a, d, &temp);
		y_hat = ((temp+g13.b) <= 0) ? 0 : 1; 

		// uncomment next line if you want to comapre with dlda.R
		//values[i] = y_hat;

		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}
    
	if (prior==2.00){
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
	  //printf("ohoh prior = %.2f",prior);
	  //printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("g13 error is %.3f\n", error);

	//uncomment this section you want to compare with dlda.R
	//PrintDouble( X, N, d, (char*)"terminal/G13_tstX.txt");
	//PrintInt( y, N, (char*)"terminal/G13_tsty.txt");
	//PrintInt( values, N , (char*)"terminal/G13_expy.txt");
	//printf("error is %f\n", error);
	//delete values;

	delete temp_X;

	return error;
}


void edcTrn(double** X, int* y, int N, int d, int* ind, double prior, EDC* edc)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d];
	double* mu_hat_1 = new double [d];
	double* temp2 = new double [d];
	double temp_b = 0.00;

	double** identityMatrix;

	// calculate amount of class0 and class1
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	// create X0 and X1 and write values to them from X
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			for (j=0; j<d; j++)
				X_0[i_0][j] = X[i][ind[j]];
			i_0++;
		}
		else
		{
			for (j=0; j<d; j++)
				X_1[i_1][j] = X[i][ind[j]];
			i_1++;
		}
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector

	sum_two_vector(mu_hat_1, mu_hat_0, d, 1.00, -1.00, 1.00, (*edc).a); // hat{X1} - hat{X0}
	sum_two_matrix(&mu_hat_0, &mu_hat_1, 1, d, 1.00, 1.00, 1.00, &temp2); // hat{X1} + hat{X0}
	multiply_a_by_b (temp2, (*edc).a, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b

	if (prior==2.00){
	  (*edc).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	} else {
	  (*edc).b = -0.50*temp_b + log(prior/(1.00-prior));
	}

	//printf("edc.b is %.3f\n", (*edc).b);
	//
	//printf("edc.a is\n");
	//for( i = 0; i < d; i++ )
	//printf("%.3f\t", (*edc).a[i]);
	//printf("\n");
	

	//uncomment these lines if you want to compare it with code implemented in R
	//printf("before clearing space\n");
	//printf("(*edc).b is %f\n\n", (*edc).b);
	//for ( int k = 0; k < d ; k++ )
	//printf("%f\t%f\n", (*edc).a[k], temp[k]);
	//
	//PrintDouble( X , N_0 + N_1 , d, (char*)"terminal/X.txt");
	//PrintDouble( X_0 , N_0 , d, (char*)"terminal/X_0.txt");
	//PrintDouble( X_1 , N_1 , d, (char*)"terminal/X_1.txt");
	//PrintDouble2( mu_hat_0 , d, (char*)"terminal/mu0.txt");
	//PrintDouble2( mu_hat_1 , d, (char*)"terminal/mu1.txt");
	//PrintInt( y, N_0 + N_1 , (char*)"terminal/y.txt");
	//PrintEDC( ((*edc).a), (*edc).b , d , (char*)"terminal/EDC.txt");
	
	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete [] mu_hat_0;
	delete [] mu_hat_1;

	return;
}

double edcTst(double** X, int* y, int N, int d, int* ind, double prior, EDC edc)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;

	// uncomment next line if you want to compare with edc.R
	//int *values = new int [N];
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, edc.a, d, &temp);
		y_hat = ((temp+edc.b) <= 0) ? 0 : 1; 

		// uncomment next line if you want to compare with edc.R
		//values[i] = y_hat;

		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}

    
	if (prior==2.00) {
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
		//printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("edc error is %.3f\n", error);

	
	//// uncomment these lines if you want to compare with code implement in R
	//PrintDouble( X, N, d, (char*)"terminal/tstX.txt");
	//PrintInt( y, N, (char*)"terminal/tsty.txt");
	//PrintInt( values, N , (char*)"terminal/expy.txt");
	//printf("error is %f\n", error);
	//delete values;

	delete temp_X;


	return error;
}



void ldaTrn(double** X, int* y, int N, int d, int* ind, double prior, LDA* lda)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d];
	double* mu_hat_1 = new double [d];
	double* temp = new double [d];
	double** cov_hat_0 = make_2D_matrix(d, d, dinit);
	double** cov_hat_1 = make_2D_matrix(d, d, dinit);
	double** pooled_cov = make_2D_matrix(d, d, dinit);
	double** inv_pooled_cov = make_2D_matrix(d, d, dinit);

	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	{
	  if (y[i]==0)
	    {
	      for (j=0; j<d; j++)
		X_0[i_0][j] = X[i][ind[j]];
	      i_0++;
	    }
	  else
	    {
	      for (j=0; j<d; j++)
		X_1[i_1][j] = X[i][ind[j]];
	      i_1++;
	    }
	}

	mean(X_0, N_0, d, mu_hat_0); //mu_hat_0 is a d dimension vector
	mean(X_1, N_1, d, mu_hat_1); //mu_hat_1 is a d dimension vector
	cov(X_0, N_0, d, cov_hat_0);
	cov(X_1, N_1, d, cov_hat_1);

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, -1.00, 1.00, &temp); // hat{X1} - hat{X0}
	sum_two_matrix(cov_hat_0, cov_hat_1, d, d, (double) (N_0-1), (double) (N_1-1), (double) (N_0+N_1-2), pooled_cov);
	invert_matrix(pooled_cov, d, inv_pooled_cov);
	multiply_A_by_b (inv_pooled_cov, temp, d, d, (*lda).a); // multiply inv{cov} by // hat{X1} - hat{X0}
															//saves the results as "a" in lda
	//multiply_A_by_b (inv_pooled_cov, mu_hat_0, d, d, temp0);
	
	/*
	 RULES: 
	 1.Note that type of any variable in a function should match the defintion of that variable.
	 e.g.  mu_hat_0 in above is defined to be double* mu_hat_0 and in 
	 mean(X_0, N_0, d, mu_hat_0) is also double* m. 
	 2. Now whenver we want to use mu_hat_0 (which is a vector double*) in situation where the matrix 
	 should be used (e.g. sum_two_matrix(double**, double**)) then we need to use &mu_hat_0 so we
	 can increment on element of mu_hat_0.
	 3. An example of rule 2 is also temp variable, which in above is defined to be double* temp. 
	 Note that in definition of sum_two_matrix, at the last position we have to have double**, so we 
	 need to use &temp. But note that in multiply_A_by_b (inv_pooled_cov, temp...) we need a 
	 double* in the second position, we we safely use "temp" itself. 
	 4. Now imagine you want to multpily inv_pooled_cov to mu_hat_0. In multiply_A_by_b we have 
	 double** in first, and double* in the second position. So we use inv_pooled_cov, which is 
	 a matrix, and then mu_hat_0, which is a double*.
	 5. The last position in multiply_A_by_b is just a vector (double*). Also temp0 is defined to be a vector
	 so we safely use "temp0" in multiply_A_by_b (inv_pooled_cov, mu_hat_0, d, d, temp0);
	 6. Below you will see that in multiply_a_by_b we have &temp_b. This is exactly as rule 2.
	 We have defined temp_b as a scalar but we need a vector (double*) in the last position of 
	 multiply_a_by_b so we use &temp_b. But then below that we want to use ( -0.50*temp_b) as a double scalar
	 so we use it itself "temp_b".
	 7. Note that in multiply_A_by_b (inv_pooled_cov, temp, d, d, (*lda).a) we have (*lda).a because we
	 have defined double* a in LDA structure so a matches the last element of multiply_A_by_b, which is 
	 double*.
	 8. In (*lda).b = -0.50*temp_b, (*lda).b is a double and on the right side we have double too.
	 */

	sum_two_matrix(&mu_hat_1, &mu_hat_0, 1, d, 1.00, 1.00, 1.00, &temp); // hat{X1} + hat{X0}
	  
	multiply_a_by_b ((*lda).a, temp, d, &temp_b); //saves (hat{X1} + hat{X0})inv{cov}(hat{X1} - hat{X0}) as temp_b
	if (prior==2.00){
	  (*lda).b = -0.50*temp_b + log((double)N_1/(double)N_0);
	} else {
	  //printf("We are here 2");
	  (*lda).b = -0.50*temp_b + log(prior/(1.00-prior));
	}

	//printf("lda.b is %f\n", (*lda).b);
	//
	//printf("lda.a is\n");
	//
	//for( i = 0 ; i < d; i++ )
	//printf("%f\t", lda->a[i]);
	//
	//printf("\n");
	

	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete [] mu_hat_0;
	delete [] mu_hat_1;
	delete [] temp;
	delete_2D_matrix(d, d, cov_hat_0);
	delete_2D_matrix(d, d, cov_hat_1);
	delete_2D_matrix(d, d, pooled_cov);
	delete_2D_matrix(d, d, inv_pooled_cov);

	return;
}

double ldaTst(double** X, int* y, int N, int d, int* ind, double prior, LDA lda)
{
	int i, j;
	double temp = 0.00;
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	temp_X = new double [d];

	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];

		multiply_a_by_b (temp_X, lda.a, d, &temp);
		y_hat = ((temp+lda.b) <= 0) ? 0 : 1; 
		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}
    
	if (prior==2.00)
	{
	  error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
		//printf("ldaerror=%.2f\n",error);
	} else {
	  error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}

	//printf("lda error is %.3f\n", error);

	delete temp_X;

	return error;
}

/////////////////////////////////////////////////////////////////////////////////
//QDA codes

void qdaTrn(double** X, int* y, int N, int d, int* ind, double prior, QDA* qda)
{
	int i, j;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double dinit = 0.00;
	double temp_b = 0.00;
	double temp_c = 0.00;
	double det0=0;
	double det1=0;
	double** X_0;
	double** X_1;
	double* mu_hat_0 = new double [d];
	double* mu_hat_1 = new double [d];
	double* temp0 = new double [d];
	double* temp1 = new double [d];
	double** cov_hat_0 = make_2D_matrix(d, d, dinit);
	double** cov_hat_1 = make_2D_matrix(d, d, dinit);
	double** inv_cov_0 = make_2D_matrix(d, d, dinit);
	double** inv_cov_1 = make_2D_matrix(d, d, dinit);
	
	
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);
	
	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			for (j=0; j<d; j++)
				X_0[i_0][j] = X[i][ind[j]];
			i_0++;
		}
		else
		{
			for (j=0; j<d; j++)
				X_1[i_1][j] = X[i][ind[j]];
			i_1++;
		}
	}
	
	mean(X_0, N_0, d, mu_hat_0);
	mean(X_1, N_1, d, mu_hat_1);
	cov(X_0, N_0, d, cov_hat_0);
	cov(X_1, N_1, d, cov_hat_1);
	
	invert_matrix(cov_hat_0, d, inv_cov_0);
	invert_matrix(cov_hat_1, d, inv_cov_1);
	multiply_A_by_b(inv_cov_0, mu_hat_0, d, d, temp0);
	multiply_A_by_b(inv_cov_1, mu_hat_1, d, d, temp1);
	sum_two_matrix(inv_cov_0, inv_cov_1, d, d, 0.5, -0.5, 1.00, (*qda).a);
	sum_two_vector(temp1, temp0, d, 1.00, -1.00, 1.00, (*qda).b);
	multiply_a_by_b (mu_hat_0, temp0, d, &temp_b);
	multiply_a_by_b (mu_hat_1, temp1, d, &temp_c);
	det0=determinant(cov_hat_0,d);
	det1=determinant(cov_hat_1,d);
	
	if (prior==2.00){
		//printf("We are here 1");
		(*qda).c = 0.50*temp_b - 0.50*temp_c - 0.5*log(det1/det0) + log((double)N_1/(double)N_0);
	} else {
		//printf("We are here 2");
		(*qda).c = 0.50*temp_b - 0.50*temp_c - 0.5*log(det1/det0) + log(prior/(1.00-prior));
	}
	
	
	
	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);
	delete mu_hat_0;
	delete mu_hat_1;
	delete temp0;
	delete temp1;
	delete_2D_matrix(d, d, cov_hat_0);
	delete_2D_matrix(d, d, cov_hat_1);
	delete_2D_matrix(d, d, inv_cov_0);
	delete_2D_matrix(d, d, inv_cov_1);
	
	return;
}


double qdaTst(double** X, int* y, int N, int d, int* ind, double prior, QDA qda)
{
	int i, j;
	double temp = 0.00;
	double temp1 = 0.00;
	double* temp_v = new double [d];
	double* temp_X;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double error = 0.00;
	int y_hat;
	int N_0, N_1;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	temp_X = new double [d];
	
	error = 0.00;
	for (i=0; i<N; i++)
	{
		for (j=0; j<d; j++)
			temp_X[j] = X[i][ind[j]];
		
		multiply_A_by_b (qda.a, temp_X, d, d, temp_v); // A^T %*% X
		multiply_a_by_b (temp_v, temp_X, d, &temp); // X %*% A^T %*% X
		multiply_a_by_b (temp_X, qda.b, d, &temp1);
		y_hat = ((temp+temp1+qda.c) <= 0) ? 0 : 1; 
		if (y_hat != y[i] && y[i]==0) {
			temp_error0++;
		}
		
		if (y_hat != y[i] && y[i]==1) {
			temp_error1++;
		}
		
	}
    
	if (prior==2.00)
	{
		error = (1.00/(double)N)*temp_error0+(1.00/(double)N)*temp_error1;
		//printf("ohoh prior = %.2f",prior);
	} else {
		error = (1.00-prior)*temp_error0/N_0+prior*temp_error1/N_1;
	}
	
	
	delete temp_X;
	
	return error;
}


void divide_a_to_k_groups ( double* x, int d, int k, int m, double* x_temp){

  // this function will split the matrix X or vector x into k groups
  // with m elements, and write them into x_temp, required for Serd

  if ( d < m*k ){
    printf("for Serdobolski classifier k*m exceeds feature size :(, please, repick\n");
    printf("d is %d\tm is %d\tk is %d\tm*k is %d\n", d, m, k, m*k);
  }

  if ( d / m != k ){
    printf("for Serdobolski classifier k*m is not the same as d, the rest is cutoff\n");
  }

  for ( int kk = 0; kk < k ; kk++ ){
    x_temp[kk] = 0;  // nullify the value just in case

    for ( int mm = kk*m ; mm < (kk+1)*m ; mm++ ) {
      x_temp[kk] += x[ mm ];
    }

    x_temp[kk] = (double) (x_temp[kk]/m);
  } 

  return;
  
}

void divide_A_to_k_groups ( double** X, int N, int d, int k, int m, double** x_temp){

  // this function will split the matrix X or vector x into k groups
  // with m elements, and write them into x_temp, required for Serd

  if ( d < m*k ){
    printf("for Serdobolski classifier k*m exceeds feature size :(, please, repick\n");
    printf("d is %d\tm is %d\tk is %d\tm*k is %d\n", d, m, k, m*k);
  }

  if ( d / m != k ){
    printf("for Serdobolski classifier k*m is not the same as d, the rest is cutoff\n");
  }

  for ( int i = 0 ; i < N; i++ ){

    for ( int kk = 0; kk < k ; kk++ ){
      x_temp[i][kk] = 0;  // nullify the value just in case

    for ( int mm = kk*m ; mm < (kk+1)*m ; mm++ ) {
      x_temp[i][kk] += X[i][ mm ];
    }

    x_temp[i][kk] = (double) (x_temp[i][kk]/((double)(m)));
    } 

  }

  return;
  
}


////////////////////////////////////////////////////////////////////////////////
//SVM code

//SVM training
//   svm_tr(X1, X2, n1, n2, p, ind, kernel_type, subdata, subcl)
//     double **X1,**X2 -> data, training, X1[i][j] means sample i feature j of class 0, X2 class 1
//     int n1,n2 -> number of observations in X1 and X2
// 	   int *ind -> index, the feature set to be used
//     int p -> feature set size
//	   svm_node **subdata: training data saved in libsvm structure, will be used in model construction/support vector buiding
//	   svm_problem *subcl: training data saved in libsvm structure, will be used in model construction/support vector buiding
//	   svm_model *svm_tr: trained svm classifier


struct svm_model *svmTrn(double** X, int* y, int N, int d, int* ind, int kernel_type, svm_node** subdata, svm_problem* subcl)
{
	int i,j;
	double dinit = 0.00;
	int N_0 = 0;
	int N_1 = 0;
	int i_0 = 0;
	int i_1 = 0;
	double** X_0;
	double** X_1;

	svm_node *subdata_tr;
	struct svm_parameter param;
	struct svm_model *model;

	model = NULL;

	///////////////////////////////////////////////
	// default values for SVM
	param.svm_type = C_SVC;
	param.degree = 3;
	param.gamma = 0;	// 1/k
	param.coef0 = 1;
	param.nu = 0.5;
	param.cache_size = 40;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	//set kernel type
	param.kernel_type = kernel_type;
	if (param.kernel_type == POLY)
	{
		param.degree = 6;
	}

	param.gamma = 1.0/d;

	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	X_0 = make_2D_matrix(N_0, d, dinit);
	X_1 = make_2D_matrix(N_1, d, dinit);

	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			for (j=0; j<d; j++)
				X_0[i_0][j] = X[i][ind[j]];
			i_0++;
		}
		else
		{
			for (j=0; j<d; j++)
				X_1[i_1][j] = X[i][ind[j]];
			i_1++;
		}
	}

	//set data
	if((*subdata=Malloc(svm_node,(N_0+N_1)*(d+1)))==NULL)
		return(NULL);
	subdata_tr = *subdata;

	subcl->l = N_0+N_1;
	//allocate memory for data vector
	subcl->x = Malloc(svm_node*, subcl->l);
	subcl->y = Malloc(double, subcl->l);

	// set the pointer
	for (i=0; i<N_0+N_1; i++)
		subcl->x[i] = &subdata_tr[i*(d+1)];


	//set training data
	for(i=0; i<N_0; i++)
	{
		for(j=0; j<d; j++)
		{
			subdata_tr[i*(d+1)+j].value = X_0[i][j];
			subdata_tr[i*(d+1)+j].index = j+1;
		}
		subdata_tr[i*(d+1)+d].index = -1;
		subcl->y[i] = 0;
	}
	for(i=0; i<N_1; i++)
	{
		for(j=0; j<d; j++)
		{
			subdata_tr[(i+N_0)*(d+1)+j].value = X_1[i][j];
			subdata_tr[(i+N_0)*(d+1)+j].index = j+1;
		}
		subdata_tr[(i+N_0)*(d+1)+d].index = -1;
		subcl->y[i+N_0] = 1;
	}

	model = svm_train(subcl, &param);

	delete_2D_matrix(N_0, d, X_0);
	delete_2D_matrix(N_1, d, X_1);

	return model;
}


//SVM testing
//   svm_tst(x, y, n, p, ind, model)
//     double **x, -> testing data
//	   int *y, -> label, testing
//     int n -> number of observations in x
// 	   int *ind -> index, the feature set to be used
//     int p -> feature set size
//	   svm_model *model -> SVM model

double svmTst(double** X, int* y, int N, int d, int* ind, double prior, svm_model* model)
{
	int i, j;
	double er0 = 0;
	double er1 = 0;
	double error = 0.00;
	int N_0, N_1;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	svm_node *subdata_ts;
	svm_problem subcl_ts;

	//set data
	if((subdata_ts = Malloc(svm_node,(d+1)))==NULL)
		return(1);

	subcl_ts.l = 1;
	//allocate memory for data vector
	subcl_ts.x = Malloc(svm_node*, subcl_ts.l);
	subcl_ts.y = Malloc(double, subcl_ts.l);

	// set the pointer
	subcl_ts.x[0] = &subdata_ts[0];


	//set test data
	for(i=0; i<N; i++)
	{
		for(j=0; j<d; j++)
		{
			subdata_ts[j].value = X[i][ind[j]];
			subdata_ts[j].index = j+1;
		}
		subdata_ts[d].index = -1;
		subcl_ts.y[0] = y[i];
        if (y[i]==0) {
		if ((int) svm_predict(model,subcl_ts.x[0]) != (int)subcl_ts.y[0])
			er0++;
		}
		if (y[i]==1) {
		if ((int) svm_predict(model,subcl_ts.x[0]) != (int)subcl_ts.y[0])
			er1++;
		}
	}
	
	
	if (prior==2.00)
	{
		error = (1.00/(double)N)*er0+(1.00/(double)N)*er1;
		//printf("ohoh prior = %.2f",prior);
	} else {
		error = (1.00-prior)*er0/N_0+prior*er1/N_1;
	}



	//free memory
	free(subcl_ts.x);
	free(subcl_ts.y);
	free(subdata_ts);

	return error;
}

//SVM destroy
//SVM finishing subroutine, release memory
// lda_destroy(model, subdata, subcl)
//	   svm_node **subdata: training data saved in libsvm structure, used in model construction/support vector buiding
//	   svm_problem *subcl: training data saved in libsvm structure, used in model construction/support vector buiding
//	   svm_model *svm_tr: trained svm classifier

void svmDestroy(svm_model* model, svm_node* subdata, svm_problem* subcl)
{
	svm_destroy_model(model);
	//free memory
	free(subcl->x);
	free(subcl->y);
	free(subdata);

}


typedef float Qfloat;
typedef signed char schar;
#ifndef min
template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
#endif
#ifndef max
template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#endif
template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
template <class S, class T> inline void clone(T*& dst, S* src, int n)
{
	dst = new T[n];
	memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
#define INF HUGE_VAL
#ifndef Malloc
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#endif

//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Cache
{
public:
	Cache(int l,int size);
	~Cache();

	// request data [0,len)
	// return some position p where [p,len) need to be filled
	// (p >= len if nothing needs to be filled)
	int get_data(const int index, Qfloat **data, int len);
	void swap_index(int i, int j);	// future_option
private:
	int l;
	int size;
	struct head_t
	{
		head_t *prev, *next;	// a cicular list
		Qfloat *data;
		int len;		// data[0,len) is cached in this entry
	};

	head_t* head;
	head_t lru_head;
	void lru_delete(head_t *h);
	void lru_insert(head_t *h);
};

Cache::Cache(int l_,int size_):l(l_),size(size_)
{
	head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
	size /= sizeof(Qfloat);
	size -= l * sizeof(head_t) / sizeof(Qfloat);
	lru_head.next = lru_head.prev = &lru_head;
}

Cache::~Cache()
{
	for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
		free(h->data);
	free(head);
}

void Cache::lru_delete(head_t *h)
{
	// delete from current location
	h->prev->next = h->next;
	h->next->prev = h->prev;
}

void Cache::lru_insert(head_t *h)
{
	// insert to last position
	h->next = &lru_head;
	h->prev = lru_head.prev;
	h->prev->next = h;
	h->next->prev = h;
}

int Cache::get_data(const int index, Qfloat **data, int len)
{
	head_t *h = &head[index];
	if(h->len) lru_delete(h);
	int more = len - h->len;

	if(more > 0)
	{
		// free old space
		while(size < more)
		{
			head_t *old = lru_head.next;
			lru_delete(old);
			free(old->data);
			size += old->len;
			old->data = 0;
			old->len = 0;
		}

		// allocate new space
		h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
		size -= more;
		swap(h->len,len);
	}

	lru_insert(h);
	*data = h->data;
	return len;
}

void Cache::swap_index(int i, int j)
{
	if(i==j) return;

	if(head[i].len) lru_delete(&head[i]);
	if(head[j].len) lru_delete(&head[j]);
	swap(head[i].data,head[j].data);
	swap(head[i].len,head[j].len);
	if(head[i].len) lru_insert(&head[i]);
	if(head[j].len) lru_insert(&head[j]);

	if(i>j) swap(i,j);
	for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
	{
		if(h->len > i)
		{
			if(h->len > j)
				swap(h->data[i],h->data[j]);
			else
			{
				// give up
				lru_delete(h);
				free(h->data);
				size += h->len;
				h->data = 0;
				h->len = 0;
			}
		}
	}
}

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class Kernel {
public:
	Kernel(int l, svm_node * const * x, const svm_parameter& param);
	virtual ~Kernel();

	static double k_function(const svm_node *x, const svm_node *y,
			const svm_parameter& param);
	virtual Qfloat *get_Q(int column, int len) const = 0;
	virtual void swap_index(int i, int j) const	// no so const...
	{
		swap(x[i],x[j]);
		if(x_square) swap(x_square[i],x_square[j]);
	}
protected:

	double (Kernel::*kernel_function)(int i, int j) const;

private:
	const svm_node **x;
	double *x_square;

	// svm_parameter
	const int kernel_type;
	const double degree;
	const double gamma;
	const double coef0;

	static double dot(const svm_node *px, const svm_node *py);
	double kernel_linear(int i, int j) const
	{
		return dot(x[i],x[j]);
	}
	double kernel_poly(int i, int j) const
	{
		return pow(gamma*dot(x[i],x[j])+coef0,degree);
	}
	double kernel_rbf(int i, int j) const
	{
		return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
	}
	double kernel_sigmoid(int i, int j) const
	{
		return tanh(gamma*dot(x[i],x[j])+coef0);
	}
};

Kernel::Kernel(int l, svm_node * const * x_, const svm_parameter& param)
:kernel_type(param.kernel_type), degree(param.degree),
 gamma(param.gamma), coef0(param.coef0)
{
	switch(kernel_type)
	{
	case LINEAR:
		kernel_function = &Kernel::kernel_linear;
		break;
	case POLY:
		kernel_function = &Kernel::kernel_poly;
		break;
	case RBF:
		kernel_function = &Kernel::kernel_rbf;
		break;
	case SIGMOID:
		kernel_function = &Kernel::kernel_sigmoid;
		break;
	}

	clone(x,x_,l);

	if(kernel_type == RBF)
	{
		x_square = new double[l];
		for(int i=0;i<l;i++)
			x_square[i] = dot(x[i],x[i]);
	}
	else
		x_square = 0;
}

Kernel::~Kernel()
{
	delete[] x;
	delete[] x_square;
}

double Kernel::dot(const svm_node *px, const svm_node *py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}
	}
	return sum;
}

double Kernel::k_function(const svm_node *x, const svm_node *y,
		const svm_parameter& param)
{
	switch(param.kernel_type)
	{
	case LINEAR:
		return dot(x,y);
	case POLY:
		return pow(param.gamma*dot(x,y)+param.coef0,param.degree);
	case RBF:
	{
		double sum = 0;
		while(x->index != -1 && y->index !=-1)
		{
			if(x->index == y->index)
			{
				double d = x->value - y->value;
				sum += d*d;
				++x;
				++y;
			}
			else
			{
				if(x->index > y->index)
				{
					sum += y->value * y->value;
					++y;
				}
				else
				{
					sum += x->value * x->value;
					++x;
				}
			}
		}

		while(x->index != -1)
		{
			sum += x->value * x->value;
			++x;
		}

		while(y->index != -1)
		{
			sum += y->value * y->value;
			++y;
		}

		return exp(-param.gamma*sum);
	}
	case SIGMOID:
		return tanh(param.gamma*dot(x,y)+param.coef0);
	default:
		return 0;	/* Unreachable */
	}
}

// Generalized SMO+SVMlight algorithm
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + b^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, b, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping criterion
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
	Solver() {};
	virtual ~Solver() {};

	struct SolutionInfo {
		double obj;
		double rho;
		double upper_bound_p;
		double upper_bound_n;
		double r;	// for Solver_NU
	};

	void Solve(int l, const Kernel& Q, const double *b_, const schar *y_,
			double *alpha_, double Cp, double Cn, double eps,
			SolutionInfo* si, int shrinking);
protected:
	int active_size;
	schar *y;
	double *G;		// gradient of objective function
	enum { LOWER_BOUND, UPPER_BOUND, FREE };
	char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
	double *alpha;
	const Kernel *Q;
	double eps;
	double Cp,Cn;
	double *b;
	int *active_set;
	double *G_bar;		// gradient, if we treat free variables as 0
	int l;
	bool unshrinked;	// XXX

	double get_C(int i)
	{
		return (y[i] > 0)? Cp : Cn;
	}
	void update_alpha_status(int i)
	{
		if(alpha[i] >= get_C(i))
			alpha_status[i] = UPPER_BOUND;
		else if(alpha[i] <= 0)
			alpha_status[i] = LOWER_BOUND;
		else alpha_status[i] = FREE;
	}
	bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
	bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
	bool is_free(int i) { return alpha_status[i] == FREE; }
	void swap_index(int i, int j);
	void reconstruct_gradient();
	virtual int select_working_set(int &i, int &j);
	virtual double calculate_rho();
	virtual void do_shrinking();
};

void Solver::swap_index(int i, int j)
{
	Q->swap_index(i,j);
	swap(y[i],y[j]);
	swap(G[i],G[j]);
	swap(alpha_status[i],alpha_status[j]);
	swap(alpha[i],alpha[j]);
	swap(b[i],b[j]);
	swap(active_set[i],active_set[j]);
	swap(G_bar[i],G_bar[j]);
}

void Solver::reconstruct_gradient()
{
	// reconstruct inactive elements of G from G_bar and free variables

	if(active_size == l) return;

	int i;
	for(i=active_size;i<l;i++)
		G[i] = G_bar[i] + b[i];

	for(i=0;i<active_size;i++)
		if(is_free(i))
		{
			const Qfloat *Q_i = Q->get_Q(i,l);
			double alpha_i = alpha[i];
			for(int j=active_size;j<l;j++)
				G[j] += alpha_i * Q_i[j];
		}
}

void Solver::Solve(int l, const Kernel& Q, const double *b_, const schar *y_,
		double *alpha_, double Cp, double Cn, double eps,
		SolutionInfo* si, int shrinking)
{
	this->l = l;
	this->Q = &Q;
	clone(b, b_,l);
	clone(y, y_,l);
	clone(alpha,alpha_,l);
	this->Cp = Cp;
	this->Cn = Cn;
	this->eps = eps;
	unshrinked = false;

	// initialize alpha_status
	{
		alpha_status = new char[l];
		for(int i=0;i<l;i++)
			update_alpha_status(i);
	}

	// initialize active set (for shrinking)
	{
		active_set = new int[l];
		for(int i=0;i<l;i++)
			active_set[i] = i;
		active_size = l;
	}

	// initialize gradient
	{
		G = new double[l];
		G_bar = new double[l];
		int i;
		for(i=0;i<l;i++)
		{
			G[i] = b[i];
			G_bar[i] = 0;
		}
		for(i=0;i<l;i++)
			if(!is_lower_bound(i))
			{
				Qfloat *Q_i = Q.get_Q(i,l);
				double alpha_i = alpha[i];
				int j;
				for(j=0;j<l;j++)
					G[j] += alpha_i*Q_i[j];
				if(is_upper_bound(i))
					for(j=0;j<l;j++)
						G_bar[j] += get_C(i) * Q_i[j];
			}
	}

	// optimization step

	int iter = 0;
	int counter = min(l,1000)+1;

	while(1)
	{
		// show progress and do shrinking

		if(--counter == 0)
		{
			counter = min(l,1000);
			if(shrinking) do_shrinking();
			//info(".");
			info_flush();
		}

		int i,j;
		if(select_working_set(i,j)!=0)
		{
			// reconstruct the whole gradient
			reconstruct_gradient();
			// reset active set size and check
			active_size = l;
			//info("*");
			info_flush();
			if(select_working_set(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}

		++iter;

		// update alpha[i] and alpha[j], handle bounds carefully

		const Qfloat *Q_i = Q.get_Q(i,active_size);
		const Qfloat *Q_j = Q.get_Q(j,active_size);

		double C_i = get_C(i);
		double C_j = get_C(j);

		double old_alpha_i = alpha[i];
		double old_alpha_j = alpha[j];

		if(y[i]!=y[j])
		{
			double delta = (-G[i]-G[j])/max(Q_i[i]+Q_j[j]+2*Q_i[j],(Qfloat)0);
			double diff = alpha[i] - alpha[j];
			alpha[i] += delta;
			alpha[j] += delta;

			if(diff > 0)
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = diff;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = -diff;
				}
			}
			if(diff > C_i - C_j)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = C_i - diff;
				}
			}
			else
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = C_j + diff;
				}
			}
		}
		else
		{
			double delta = (G[i]-G[j])/max(Q_i[i]+Q_j[j]-2*Q_i[j],(Qfloat)0);
			double sum = alpha[i] + alpha[j];
			alpha[i] -= delta;
			alpha[j] += delta;
			if(sum > C_i)
			{
				if(alpha[i] > C_i)
				{
					alpha[i] = C_i;
					alpha[j] = sum - C_i;
				}
			}
			else
			{
				if(alpha[j] < 0)
				{
					alpha[j] = 0;
					alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(alpha[j] > C_j)
				{
					alpha[j] = C_j;
					alpha[i] = sum - C_j;
				}
			}
			else
			{
				if(alpha[i] < 0)
				{
					alpha[i] = 0;
					alpha[j] = sum;
				}
			}
		}

		// update G

		double delta_alpha_i = alpha[i] - old_alpha_i;
		double delta_alpha_j = alpha[j] - old_alpha_j;

		for(int k=0;k<active_size;k++)
		{
			G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
		}

		// update alpha_status and G_bar

		{
			bool ui = is_upper_bound(i);
			bool uj = is_upper_bound(j);
			update_alpha_status(i);
			update_alpha_status(j);
			int k;
			if(ui != is_upper_bound(i))
			{
				Q_i = Q.get_Q(i,l);
				if(ui)
					for(k=0;k<l;k++)
						G_bar[k] -= C_i * Q_i[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_i * Q_i[k];
			}

			if(uj != is_upper_bound(j))
			{
				Q_j = Q.get_Q(j,l);
				if(uj)
					for(k=0;k<l;k++)
						G_bar[k] -= C_j * Q_j[k];
				else
					for(k=0;k<l;k++)
						G_bar[k] += C_j * Q_j[k];
			}
		}
	}

	// calculate rho

	si->rho = calculate_rho();

	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;i++)
			v += alpha[i] * (G[i] + b[i]);

		si->obj = v/2;
	}

	// put back the solution
	{
		for(int i=0;i<l;i++)
			alpha_[active_set[i]] = alpha[i];
	}

	// juggle everything back
	/*{
		for(int i=0;i<l;i++)
			while(active_set[i] != i)
				swap_index(i,active_set[i]);
				// or Q.swap_index(i,active_set[i]);
	}*/

	si->upper_bound_p = Cp;
	si->upper_bound_n = Cn;

	//	info("\noptimization finished, #iter = %d\n",iter);

	delete[] b;
	delete[] y;
	delete[] alpha;
	delete[] alpha_status;
	delete[] active_set;
	delete[] G;
	delete[] G_bar;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;		// max { -grad(f)_i * d | y_i*d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;		// max { -grad(f)_i * d | y_i*d = -1 }
	int Gmax2_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y = +1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y = -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax2)
				{
					Gmax2 = -G[i];
					Gmax2_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax1)
				{
					Gmax1 = G[i];
					Gmax1_idx = i;
				}
			}
		}
	}

	if(Gmax1+Gmax2 < eps)
		return 1;

	out_i = Gmax1_idx;
	out_j = Gmax2_idx;
	return 0;
}

void Solver::do_shrinking()
{
	int i,j,k;
	if(select_working_set(i,j)!=0) return;
	double Gm1 = -y[j]*G[j];
	double Gm2 = y[i]*G[i];

	// shrink

	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm1) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || -(Gm1 + Gm2) > eps*10) return;

	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm2) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm1) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver::calculate_rho()
{
	double r;
	int nr_free = 0;
	double ub = INF, lb = -INF, sum_free = 0;
	for(int i=0;i<active_size;i++)
	{
		double yG = y[i]*G[i];

		if(is_lower_bound(i))
		{
			if(y[i] > 0)
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else if(is_upper_bound(i))
		{
			if(y[i] < 0)
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else
		{
			++nr_free;
			sum_free += yG;
		}
	}

	if(nr_free>0)
		r = sum_free/nr_free;
	else
		r = (ub+lb)/2;

	return r;
}

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU : public Solver
{
public:
	Solver_NU() {}
	void Solve(int l, const Kernel& Q, const double *b, const schar *y,
			double *alpha, double Cp, double Cn, double eps,
			SolutionInfo* si, int shrinking)
	{
		this->si = si;
		Solver::Solve(l,Q,b,y,alpha,Cp,Cn,eps,si,shrinking);
	}
private:
	SolutionInfo *si;
	int select_working_set(int &i, int &j);
	double calculate_rho();
	void do_shrinking();
};

int Solver_NU::select_working_set(int &out_i, int &out_j)
{
	// return i,j which maximize -grad(f)^T d , under constraint
	// if alpha_i == C, d != +1
	// if alpha_i == 0, d != -1

	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	int Gmax2_idx = -1;

	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	int Gmax3_idx = -1;

	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }
	int Gmax4_idx = -1;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)	// y == +1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax1)
				{
					Gmax1 = -G[i];
					Gmax1_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax2)
				{
					Gmax2 = G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y == -1
		{
			if(!is_upper_bound(i))	// d = +1
			{
				if(-G[i] > Gmax3)
				{
					Gmax3 = -G[i];
					Gmax3_idx = i;
				}
			}
			if(!is_lower_bound(i))	// d = -1
			{
				if(G[i] > Gmax4)
				{
					Gmax4 = G[i];
					Gmax4_idx = i;
				}
			}
		}
	}

	if(max(Gmax1+Gmax2,Gmax3+Gmax4) < eps)
		return 1;

	if(Gmax1+Gmax2 > Gmax3+Gmax4)
	{
		out_i = Gmax1_idx;
		out_j = Gmax2_idx;
	}
	else
	{
		out_i = Gmax3_idx;
		out_j = Gmax4_idx;
	}
	return 0;
}

void Solver_NU::do_shrinking()
{
	double Gmax1 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = +1 }
	double Gmax2 = -INF;	// max { -grad(f)_i * d | y_i = +1, d = -1 }
	double Gmax3 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = +1 }
	double Gmax4 = -INF;	// max { -grad(f)_i * d | y_i = -1, d = -1 }

	int k;
	for(k=0;k<active_size;k++)
	{
		if(!is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] > Gmax1) Gmax1 = -G[k];
			}
			else	if(-G[k] > Gmax3) Gmax3 = -G[k];
		}
		if(!is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] > Gmax2) Gmax2 = G[k];
			}
			else	if(G[k] > Gmax4) Gmax4 = G[k];
		}
	}

	double Gm1 = -Gmax2;
	double Gm2 = -Gmax1;
	double Gm3 = -Gmax4;
	double Gm4 = -Gmax3;

	for(k=0;k<active_size;k++)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] >= Gm1) continue;
			}
			else	if(-G[k] >= Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] >= Gm2) continue;
			}
			else	if(G[k] >= Gm4) continue;
		}
		else continue;

		--active_size;
		swap_index(k,active_size);
		--k;	// look at the newcomer
	}

	// unshrink, check all variables again before final iterations

	if(unshrinked || max(-(Gm1+Gm2),-(Gm3+Gm4)) > eps*10) return;

	unshrinked = true;
	reconstruct_gradient();

	for(k=l-1;k>=active_size;k--)
	{
		if(is_lower_bound(k))
		{
			if(y[k]==+1)
			{
				if(-G[k] < Gm1) continue;
			}
			else	if(-G[k] < Gm3) continue;
		}
		else if(is_upper_bound(k))
		{
			if(y[k]==+1)
			{
				if(G[k] < Gm2) continue;
			}
			else	if(G[k] < Gm4) continue;
		}
		else continue;

		swap_index(k,active_size);
		active_size++;
		++k;	// look at the newcomer
	}
}

double Solver_NU::calculate_rho()
{
	int nr_free1 = 0,nr_free2 = 0;
	double ub1 = INF, ub2 = INF;
	double lb1 = -INF, lb2 = -INF;
	double sum_free1 = 0, sum_free2 = 0;

	for(int i=0;i<active_size;i++)
	{
		if(y[i]==+1)
		{
			if(is_lower_bound(i))
				ub1 = min(ub1,G[i]);
			else if(is_upper_bound(i))
				lb1 = max(lb1,G[i]);
			else
			{
				++nr_free1;
				sum_free1 += G[i];
			}
		}
		else
		{
			if(is_lower_bound(i))
				ub2 = min(ub2,G[i]);
			else if(is_upper_bound(i))
				lb2 = max(lb2,G[i]);
			else
			{
				++nr_free2;
				sum_free2 += G[i];
			}
		}
	}

	double r1,r2;
	if(nr_free1 > 0)
		r1 = sum_free1/nr_free1;
	else
		r1 = (ub1+lb1)/2;

	if(nr_free2 > 0)
		r2 = sum_free2/nr_free2;
	else
		r2 = (ub2+lb2)/2;

	si->r = (r1+r2)/2;
	return (r1-r2)/2;
}

//
// Q matrices for various formulations
//
class SVC_Q: public Kernel
{
public:
	SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
	:Kernel(prob.l, prob.x, param)
	{
		clone(y,y_,prob.l);
		cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	}

	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(int j=start;j<len;j++)
				data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
		}
		return data;
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		Kernel::swap_index(i,j);
		swap(y[i],y[j]);
	}

	~SVC_Q()
	{
		delete[] y;
		delete cache;
	}
private:
	schar *y;
	Cache *cache;
};

class ONE_CLASS_Q: public Kernel
{
public:
	ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param)
	{
		cache = new Cache(prob.l,(int)(param.cache_size*(1<<20)));
	}

	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int start;
		if((start = cache->get_data(i,&data,len)) < len)
		{
			for(int j=start;j<len;j++)
				data[j] = (Qfloat)(this->*kernel_function)(i,j);
		}
		return data;
	}

	void swap_index(int i, int j) const
	{
		cache->swap_index(i,j);
		Kernel::swap_index(i,j);
	}

	~ONE_CLASS_Q()
	{
		delete cache;
	}
private:
	Cache *cache;
};

class SVR_Q: public Kernel
{
public:
	SVR_Q(const svm_problem& prob, const svm_parameter& param)
	:Kernel(prob.l, prob.x, param)
	{
		l = prob.l;
		cache = new Cache(l,(int)(param.cache_size*(1<<20)));
		sign = new schar[2*l];
		index = new int[2*l];
		for(int k=0;k<l;k++)
		{
			sign[k] = 1;
			sign[k+l] = -1;
			index[k] = k;
			index[k+l] = k;
		}
		buffer[0] = new Qfloat[2*l];
		buffer[1] = new Qfloat[2*l];
		next_buffer = 0;
	}

	void swap_index(int i, int j) const
	{
		swap(sign[i],sign[j]);
		swap(index[i],index[j]);
	}

	Qfloat *get_Q(int i, int len) const
	{
		Qfloat *data;
		int real_i = index[i];
		if(cache->get_data(real_i,&data,l) < l)
		{
			for(int j=0;j<l;j++)
				data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
		}

		// reorder and copy
		Qfloat *buf = buffer[next_buffer];
		next_buffer = 1 - next_buffer;
		schar si = sign[i];
		for(int j=0;j<len;j++)
			buf[j] = si * sign[j] * data[index[j]];
		return buf;
	}

	~SVR_Q()
	{
		delete cache;
		delete[] sign;
		delete[] index;
		delete[] buffer[0];
		delete[] buffer[1];
	}
private:
	int l;
	Cache *cache;
	schar *sign;
	int *index;
	mutable int next_buffer;
	Qfloat* buffer[2];
};

//
// construct and solve various formulations
//
static void solve_c_svc(
		const svm_problem *prob, const svm_parameter* param,
		double *alpha, Solver::SolutionInfo* si, double Cp, double Cn)
{
	int l = prob->l;
	double *minus_ones = new double[l];
	schar *y = new schar[l];

	int i;

	for(i=0;i<l;i++)
	{
		alpha[i] = 0;
		minus_ones[i] = -1;
		if(prob->y[i] > 0) y[i] = +1; else y[i]=-1;
	}

	Solver s;
	s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y,
			alpha, Cp, Cn, param->eps, si, param->shrinking);

	double sum_alpha=0;
	for(i=0;i<l;i++)
		sum_alpha += alpha[i];

	//info("nu = %f\n", sum_alpha/(param->C*prob->l));

	for(i=0;i<l;i++)
		alpha[i] *= y[i];

	delete[] minus_ones;
	delete[] y;
}

static void solve_nu_svc(
		const svm_problem *prob, const svm_parameter *param,
		double *alpha, Solver::SolutionInfo* si)
{
	int i;
	int l = prob->l;
	double nu = param->nu;

	schar *y = new schar[l];

	for(i=0;i<l;i++)
		if(prob->y[i]>0)
			y[i] = +1;
		else
			y[i] = -1;

	double sum_pos = nu*l/2;
	double sum_neg = nu*l/2;

	for(i=0;i<l;i++)
		if(y[i] == +1)
		{
			alpha[i] = min(1.0,sum_pos);
			sum_pos -= alpha[i];
		}
		else
		{
			alpha[i] = min(1.0,sum_neg);
			sum_neg -= alpha[i];
		}

	double *zeros = new double[l];

	for(i=0;i<l;i++)
		zeros[i] = 0;

	Solver_NU s;
	s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,
			alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
	double r = si->r;

	//info("C = %f\n",1/r);

	for(i=0;i<l;i++)
		alpha[i] *= y[i]/r;

	si->rho /= r;
	si->obj /= (r*r);
	si->upper_bound_p = 1/r;
	si->upper_bound_n = 1/r;

	delete[] y;
	delete[] zeros;
}

static void solve_one_class(
		const svm_problem *prob, const svm_parameter *param,
		double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double *zeros = new double[l];
	schar *ones = new schar[l];
	int i;

	int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound

	for(i=0;i<n;i++)
		alpha[i] = 1;
	alpha[n] = param->nu * prob->l - n;
	for(i=n+1;i<l;i++)
		alpha[i] = 0;

	for(i=0;i<l;i++)
	{
		zeros[i] = 0;
		ones[i] = 1;
	}

	Solver s;
	s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
			alpha, 1.0, 1.0, param->eps, si, param->shrinking);

	delete[] zeros;
	delete[] ones;
}

static void solve_epsilon_svr(
		const svm_problem *prob, const svm_parameter *param,
		double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
	int i;

	for(i=0;i<l;i++)
	{
		alpha2[i] = 0;
		linear_term[i] = param->p - prob->y[i];
		y[i] = 1;

		alpha2[i+l] = 0;
		linear_term[i+l] = param->p + prob->y[i];
		y[i+l] = -1;
	}

	Solver s;
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
			alpha2, param->C, param->C, param->eps, si, param->shrinking);

	double sum_alpha = 0;
	for(i=0;i<l;i++)
	{
		alpha[i] = alpha2[i] - alpha2[i+l];
		sum_alpha += fabs(alpha[i]);
	}
	//info("nu = %f\n",sum_alpha/(param->C*l));

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}

static void solve_nu_svr(
		const svm_problem *prob, const svm_parameter *param,
		double *alpha, Solver::SolutionInfo* si)
{
	int l = prob->l;
	double C = param->C;
	double *alpha2 = new double[2*l];
	double *linear_term = new double[2*l];
	schar *y = new schar[2*l];
	int i;

	double sum = C * param->nu * l / 2;
	for(i=0;i<l;i++)
	{
		alpha2[i] = alpha2[i+l] = min(sum,C);
		sum -= alpha2[i];

		linear_term[i] = - prob->y[i];
		y[i] = 1;

		linear_term[i+l] = prob->y[i];
		y[i+l] = -1;
	}

	Solver_NU s;
	s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
			alpha2, C, C, param->eps, si, param->shrinking);

	//info("epsilon = %f\n",-si->r);

	for(i=0;i<l;i++)
		alpha[i] = alpha2[i] - alpha2[i+l];

	delete[] alpha2;
	delete[] linear_term;
	delete[] y;
}

//
// decision_function
//
struct decision_function
{
	double *alpha;
	double rho;
};

decision_function svm_train_one(
		const svm_problem *prob, const svm_parameter *param,
		double Cp, double Cn)
{
	double *alpha = Malloc(double,prob->l);
	Solver::SolutionInfo si;
	switch(param->svm_type)
	{
	case C_SVC:
		solve_c_svc(prob,param,alpha,&si,Cp,Cn);
		break;
	case NU_SVC:
		solve_nu_svc(prob,param,alpha,&si);
		break;
	case ONE_CLASS:
		solve_one_class(prob,param,alpha,&si);
		break;
	case EPSILON_SVR:
		solve_epsilon_svr(prob,param,alpha,&si);
		break;
	case NU_SVR:
		solve_nu_svr(prob,param,alpha,&si);
		break;
	}

	//	info("obj = %f, rho = %f\n",si.obj,si.rho);

	// output SVs

	int nSV = 0;
	int nBSV = 0;
	for(int i=0;i<prob->l;i++)
	{
		if(fabs(alpha[i]) > 0)
		{
			++nSV;
			if(prob->y[i] > 0)
			{
				if(fabs(alpha[i]) >= si.upper_bound_p)
					++nBSV;
			}
			else
			{
				if(fabs(alpha[i]) >= si.upper_bound_n)
					++nBSV;
			}
		}
	}

	//	info("nSV = %d, nBSV = %d\n",nSV,nBSV);

	decision_function f;
	f.alpha = alpha;
	f.rho = si.rho;
	return f;
}


//
// Interface functions
//
svm_model *svm_train(const svm_problem *prob, const svm_parameter *param)
{
	svm_model *model = Malloc(svm_model,1);
	model->param = *param;
	model->free_sv = 0;	// XXX

	if(param->svm_type == ONE_CLASS ||
			param->svm_type == EPSILON_SVR ||
			param->svm_type == NU_SVR)
	{
		// regression or one-class-svm
		model->nr_class = 2;
		model->label = NULL;
		model->nSV = NULL;
		model->sv_coef = Malloc(double *,1);
		decision_function f = svm_train_one(prob,param,0,0);
		model->rho = Malloc(double,1);
		model->rho[0] = f.rho;

		int nSV = 0;
		int i;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0) ++nSV;
		model->l = nSV;
		model->SV = Malloc(svm_node *,nSV);
		model->sv_coef[0] = Malloc(double,nSV);
		int j = 0;
		for(i=0;i<prob->l;i++)
			if(fabs(f.alpha[i]) > 0)
			{
				model->SV[j] = prob->x[i];
				model->sv_coef[0][j] = f.alpha[i];
				++j;
			}

		free(f.alpha);
	}
	else
	{
		// classification
		// find out the number of classes
		int l = prob->l;
		int max_nr_class = 16;
		int nr_class = 0;
		int *label = Malloc(int,max_nr_class);
		int *count = Malloc(int,max_nr_class);
		int *index = Malloc(int,l);

		int i;
		for(i=0;i<l;i++)
		{
			int this_label = (int)prob->y[i];
			int j;
			for(j=0;j<nr_class;j++)
				if(this_label == label[j])
				{
					++count[j];
					break;
				}
			index[i] = j;
			if(j == nr_class)
			{
				if(nr_class == max_nr_class)
				{
					max_nr_class *= 2;
					label = (int *)realloc(label,max_nr_class*sizeof(int));
					count = (int *)realloc(count,max_nr_class*sizeof(int));
				}
				label[nr_class] = this_label;
				count[nr_class] = 1;
				++nr_class;
			}
		}

		// group training data of the same class

		int *start = Malloc(int,nr_class);
		start[0] = 0;
		for(i=1;i<nr_class;i++)
			start[i] = start[i-1]+count[i-1];

		svm_node **x = Malloc(svm_node *,l);

		for(i=0;i<l;i++)
		{
			x[start[index[i]]] = prob->x[i];
			++start[index[i]];
		}

		start[0] = 0;
		for(i=1;i<nr_class;i++)
			start[i] = start[i-1]+count[i-1];

		// calculate weighted C

		double *weighted_C = Malloc(double, nr_class);
		for(i=0;i<nr_class;i++)
			weighted_C[i] = param->C;
		for(i=0;i<param->nr_weight;i++)
		{
			int j;
			for(j=0;j<nr_class;j++)
				if(param->weight_label[i] == label[j])
					break;
			if(j == nr_class)
				fprintf(stderr,"warning: class label %d specified in weight is not found\n", param->weight_label[i]);
			else
				weighted_C[j] *= param->weight[i];
		}

		// train n*(n-1)/2 models

		bool *nonzero = Malloc(bool,l);
		for(i=0;i<l;i++)
			nonzero[i] = false;
		decision_function *f = Malloc(decision_function,nr_class*(nr_class-1)/2);

		int p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				svm_problem sub_prob;
				int si = start[i], sj = start[j];
				int ci = count[i], cj = count[j];
				sub_prob.l = ci+cj;
				sub_prob.x = Malloc(svm_node *,sub_prob.l);
				sub_prob.y = Malloc(double,sub_prob.l);
				int k;
				for(k=0;k<ci;k++)
				{
					sub_prob.x[k] = x[si+k];
					sub_prob.y[k] = +1;
				}
				for(k=0;k<cj;k++)
				{
					sub_prob.x[ci+k] = x[sj+k];
					sub_prob.y[ci+k] = -1;
				}

				f[p] = svm_train_one(&sub_prob,param,weighted_C[i],weighted_C[j]);
				for(k=0;k<ci;k++)
					if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
						nonzero[si+k] = true;
				for(k=0;k<cj;k++)
					if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
						nonzero[sj+k] = true;
				free(sub_prob.x);
				free(sub_prob.y);
				++p;
			}

		// build output

		model->nr_class = nr_class;

		model->label = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
			model->label[i] = label[i];

		model->rho = Malloc(double,nr_class*(nr_class-1)/2);
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			model->rho[i] = f[i].rho;

		int total_sv = 0;
		int *nz_count = Malloc(int,nr_class);
		model->nSV = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
		{
			int nSV = 0;
			for(int j=0;j<count[i];j++)
				if(nonzero[start[i]+j])
				{
					++nSV;
					++total_sv;
				}
			model->nSV[i] = nSV;
			nz_count[i] = nSV;
		}

		//		info("Total nSV = %d\n",total_sv);

		model->l = total_sv;
		model->SV = Malloc(svm_node *,total_sv);
		p = 0;
		for(i=0;i<l;i++)
			if(nonzero[i]) model->SV[p++] = x[i];

		int *nz_start = Malloc(int,nr_class);
		nz_start[0] = 0;
		for(i=1;i<nr_class;i++)
			nz_start[i] = nz_start[i-1]+nz_count[i-1];

		model->sv_coef = Malloc(double *,nr_class-1);
		for(i=0;i<nr_class-1;i++)
			model->sv_coef[i] = Malloc(double,total_sv);

		p = 0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				// classifier (i,j): coefficients with
				// i are in sv_coef[j-1][nz_start[i]...],
				// j are in sv_coef[i][nz_start[j]...]

				int si = start[i];
				int sj = start[j];
				int ci = count[i];
				int cj = count[j];

				int q = nz_start[i];
				int k;
				for(k=0;k<ci;k++)
					if(nonzero[si+k])
						model->sv_coef[j-1][q++] = f[p].alpha[k];
				q = nz_start[j];
				for(k=0;k<cj;k++)
					if(nonzero[sj+k])
						model->sv_coef[i][q++] = f[p].alpha[ci+k];
				++p;
			}

		free(label);
		free(count);
		free(index);
		free(start);
		free(x);
		free(weighted_C);
		free(nonzero);
		for(i=0;i<nr_class*(nr_class-1)/2;i++)
			free(f[i].alpha);
		free(f);
		free(nz_count);
		free(nz_start);
	}
	return model;
}

double svm_predict(const svm_model *model, const svm_node *x)
{
	if(model->param.svm_type == ONE_CLASS ||
			model->param.svm_type == EPSILON_SVR ||
			model->param.svm_type == NU_SVR)
	{
		double *sv_coef = model->sv_coef[0];
		double sum = 0;
		for(int i=0;i<model->l;i++)
			sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
		sum -= model->rho[0];
		if(model->param.svm_type == ONE_CLASS)
			return (sum>0)?1:-1;
		else
			return sum;
	}
	else
	{
		int i;
		int nr_class = model->nr_class;
		int l = model->l;

		double *kvalue = Malloc(double,l);
		for(i=0;i<l;i++)
			kvalue[i] = Kernel::k_function(x,model->SV[i],model->param);

		int *start = Malloc(int,nr_class);
		start[0] = 0;
		for(i=1;i<nr_class;i++)
			start[i] = start[i-1]+model->nSV[i-1];

		int *vote = Malloc(int,nr_class);
		for(i=0;i<nr_class;i++)
			vote[i] = 0;
		int p=0;
		for(i=0;i<nr_class;i++)
			for(int j=i+1;j<nr_class;j++)
			{
				double sum = 0;
				int si = start[i];
				int sj = start[j];
				int ci = model->nSV[i];
				int cj = model->nSV[j];

				int k;
				double *coef1 = model->sv_coef[j-1];
				double *coef2 = model->sv_coef[i];
				for(k=0;k<ci;k++)
					sum += coef1[si+k] * kvalue[si+k];
				for(k=0;k<cj;k++)
					sum += coef2[sj+k] * kvalue[sj+k];
				sum -= model->rho[p++];
				if(sum > 0)
					++vote[i];
				else
					++vote[j];
			}

		int vote_max_idx = 0;
		for(i=1;i<nr_class;i++)
			if(vote[i] > vote[vote_max_idx])
				vote_max_idx = i;
		free(kvalue);
		free(start);
		free(vote);
		return model->label[vote_max_idx];
	}
}

const char *svm_type_table[] =
{
		"c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL
};

const char *kernel_type_table[]=
{
		"linear","polynomial","rbf","sigmoid",NULL
};

int svm_save_model(const char *model_file_name, const svm_model *model)
{
	FILE *fp = fopen(model_file_name,"w");
	if(fp==NULL) return -1;

	const svm_parameter& param = model->param;

	fprintf(fp,"svm_type %s\n", svm_type_table[param.svm_type]);
	fprintf(fp,"kernel_type %s\n", kernel_type_table[param.kernel_type]);

	if(param.kernel_type == POLY)
		fprintf(fp,"degree %g\n", param.degree);

	if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
		fprintf(fp,"gamma %g\n", param.gamma);

	if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
		fprintf(fp,"coef0 %g\n", param.coef0);

	int nr_class = model->nr_class;
	int l = model->l;
	fprintf(fp, "nr_class %d\n", nr_class);
	fprintf(fp, "total_sv %d\n",l);

	{
		fprintf(fp, "rho");
		for(int i=0;i<nr_class*(nr_class-1)/2;i++)
			fprintf(fp," %g",model->rho[i]);
		fprintf(fp, "\n");
	}

	if(model->label)
	{
		fprintf(fp, "label");
		for(int i=0;i<nr_class;i++)
			fprintf(fp," %d",model->label[i]);
		fprintf(fp, "\n");
	}

	if(model->nSV)
	{
		fprintf(fp, "nr_sv");
		for(int i=0;i<nr_class;i++)
			fprintf(fp," %d",model->nSV[i]);
		fprintf(fp, "\n");
	}

	fprintf(fp, "SV\n");
	const double * const *sv_coef = model->sv_coef;
	const svm_node * const *SV = model->SV;

	for(int i=0;i<l;i++)
	{
		for(int j=0;j<nr_class-1;j++)
			fprintf(fp, "%.16g ",sv_coef[j][i]);

		const svm_node *p = SV[i];
		while(p->index != -1)
		{
			fprintf(fp,"%d:%.8g ",p->index,p->value);
			p++;
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return 0;
}

int get_sv(const svm_model *model)
{
	return model->l;
}

svm_model *svm_load_model(const char *model_file_name)
{
	FILE *fp = fopen(model_file_name,"rb");
	if(fp==NULL) return NULL;

	// read parameters

	svm_model *model = Malloc(svm_model,1);
	svm_parameter& param = model->param;
	model->rho = NULL;
	model->label = NULL;
	model->nSV = NULL;

	char cmd[81];
	while(1)
	{
		fscanf(fp,"%80s",cmd);

		if(strcmp(cmd,"svm_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;svm_type_table[i];i++)
			{
				if(strcmp(svm_type_table[i],cmd)==0)
				{
					param.svm_type=i;
					break;
				}
			}
			if(svm_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown svm type.\n");
				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"kernel_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;kernel_type_table[i];i++)
			{
				if(strcmp(kernel_type_table[i],cmd)==0)
				{
					param.kernel_type=i;
					break;
				}
			}
			if(kernel_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown kernel function.\n");
				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"degree")==0)
			fscanf(fp,"%lf",&param.degree);
		else if(strcmp(cmd,"gamma")==0)
			fscanf(fp,"%lf",&param.gamma);
		else if(strcmp(cmd,"coef0")==0)
			fscanf(fp,"%lf",&param.coef0);
		else if(strcmp(cmd,"nr_class")==0)
			fscanf(fp,"%d",&model->nr_class);
		else if(strcmp(cmd,"total_sv")==0)
			fscanf(fp,"%d",&model->l);
		else if(strcmp(cmd,"rho")==0)
		{
			int n = model->nr_class * (model->nr_class-1)/2;
			model->rho = Malloc(double,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%lf",&model->rho[i]);
		}
		else if(strcmp(cmd,"label")==0)
		{
			int n = model->nr_class;
			model->label = Malloc(int,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%d",&model->label[i]);
		}
		else if(strcmp(cmd,"nr_sv")==0)
		{
			int n = model->nr_class;
			model->nSV = Malloc(int,n);
			for(int i=0;i<n;i++)
				fscanf(fp,"%d",&model->nSV[i]);
		}
		else if(strcmp(cmd,"SV")==0)
		{
			while(1)
			{
				int c = getc(fp);
				if(c==EOF || c=='\n') break;
			}
			break;
		}
		else
		{
			fprintf(stderr,"unknown text in model file\n");
			free(model->rho);
			free(model->label);
			free(model->nSV);
			free(model);
			return NULL;
		}
	}

	// read sv_coef and SV

	int elements = 0;
	long pos = ftell(fp);

	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
		case '\n':
			// count the '-1' element
		case ':':
			++elements;
			break;
		case EOF:
			goto out;
		default:
			;
		}
	}
	out:
	fseek(fp,pos,SEEK_SET);

	int m = model->nr_class - 1;
	int l = model->l;
	model->sv_coef = Malloc(double *,m);
	int i;
	for(i=0;i<m;i++)
		model->sv_coef[i] = Malloc(double,l);
	model->SV = Malloc(svm_node*,l);
	svm_node *x_space = Malloc(svm_node,elements);

	int j=0;
	for(i=0;i<l;i++)
	{
		model->SV[i] = &x_space[j];
		for(int k=0;k<m;k++)
			fscanf(fp,"%lf",&model->sv_coef[k][i]);
		while(1)
		{
			int c;
			do {
				c = getc(fp);
				if(c=='\n') goto out2;
			} while(isspace(c));
			ungetc(c,fp);
			fscanf(fp,"%d:%lf",&(x_space[j].index),&(x_space[j].value));
			++j;
		}
		out2:
		x_space[j++].index = -1;
	}

	fclose(fp);

	model->free_sv = 1;	// XXX
	return model;
}

void svm_destroy_model(svm_model* model)
{
	if(model->free_sv)
		free((void *)(model->SV[0]));
	for(int i=0;i<model->nr_class-1;i++)
		free(model->sv_coef[i]);
	free(model->SV);
	free(model->sv_coef);
	free(model->rho);
	free(model->label);
	free(model->nSV);
	free(model);
}

const char *svm_check_parameter(const svm_problem *prob, const svm_parameter *param)
{
	// svm_type

	int svm_type = param->svm_type;
	if(svm_type != C_SVC &&
			svm_type != NU_SVC &&
			svm_type != ONE_CLASS &&
			svm_type != EPSILON_SVR &&
			svm_type != NU_SVR)
		return "unknown svm type";

	// kernel_type

	int kernel_type = param->kernel_type;
	if(kernel_type != LINEAR &&
			kernel_type != POLY &&
			kernel_type != RBF &&
			kernel_type != SIGMOID)
		return "unknown kernel type";

	// cache_size,eps,C,nu,p,shrinking

	if(param->cache_size <= 0)
		return "cache_size <= 0";

	if(param->eps <= 0)
		return "eps <= 0";

	if(svm_type == C_SVC ||
			svm_type == EPSILON_SVR ||
			svm_type == NU_SVR)
		if(param->C <= 0)
			return "C <= 0";

	if(svm_type == NU_SVC ||
			svm_type == ONE_CLASS ||
			svm_type == NU_SVR)
		if(param->nu < 0 || param->nu > 1)
			return "nu < 0 or nu > 1";

	if(svm_type == EPSILON_SVR)
		if(param->p < 0)
			return "p < 0";

	if(param->shrinking != 0 &&
			param->shrinking != 1)
		return "shrinking != 0 and shrinking != 1";


	// check whether nu-svc is feasible

	if(svm_type == NU_SVC)
	{
		int l = prob->l;
		int max_nr_class = 16;
		int nr_class = 0;
		int *label = Malloc(int,max_nr_class);
		int *count = Malloc(int,max_nr_class);

		int i;
		for(i=0;i<l;i++)
		{
			int this_label = (int)prob->y[i];
			int j;
			for(j=0;j<nr_class;j++)
				if(this_label == label[j])
				{
					++count[j];
					break;
				}
			if(j == nr_class)
			{
				if(nr_class == max_nr_class)
				{
					max_nr_class *= 2;
					label = (int *)realloc(label,max_nr_class*sizeof(int));
					count = (int *)realloc(count,max_nr_class*sizeof(int));
				}
				label[nr_class] = this_label;
				count[nr_class] = 1;
				++nr_class;
			}
		}

		for(i=0;i<nr_class;i++)
		{
			int n1 = count[i];
			for(int j=i+1;j<nr_class;j++)
			{
				int n2 = count[j];
				if(param->nu*(n1+n2)/2 > min(n1,n2))
				{
					free(label);
					free(count);
					return "specified nu is infeasible";
				}
			}
		}
	}

	return NULL;
}

void info(char *fmt,...)
{
	va_list ap;
	va_start(ap,fmt);
	vprintf(fmt,ap);
	va_end(ap);
}
void info_flush()
{
	fflush(stdout);
}

int myComp(const void* a, const void* b) {
  Edge* a1 = (Edge*)a;
  Edge* b1 = (Edge*)b;
  return (*a1).weight < (*b1).weight;
}

// Creates a graph with V vertices and E edges
Graph* createGraph(int V, int E)
{
  //Graph* graph = (Graph*) malloc( sizeof(struct Graph) );
  Graph* graph = new Graph;
  (*graph).V = V;
  (*graph).E = E;
  
  //graph->edge = (struct Edge*) malloc( graph->E * sizeof( struct Edge ) );
  (*graph).edge = new Edge [(*graph).E];
 
  return graph;
}

int find(subset subsets[], int i)
{
  //printf("subset[i].parent is %d\ti is %d\n", subsets[i].parent, i);
    // find root and make root as parent of i (path compression)
    if (subsets[i].parent != i)
      subsets[i].parent = find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

void Union(subset subsets[], int x, int y)
{
  //printf("\nfind for x union\n");
  int xroot = find(subsets, x);

  //printf("\nfind for y union\n");
  int yroot = find(subsets, y);
  
  // Attach smaller rank tree under root of high rank tree
  // (Union by Rank)
  if (subsets[xroot].rank < subsets[yroot].rank)
    subsets[xroot].parent = yroot;
  else if (subsets[xroot].rank > subsets[yroot].rank)
    subsets[yroot].parent = xroot;
  
  // If ranks are same, then make one as root and increment
  // its rank by one
  else
    {
      subsets[yroot].parent = xroot;
      subsets[xroot].rank++;
    }
}

// The main function to construct MST using Kruskal's algorithm
//void KruskalMST(Graph* graph, double** totalCovMatrix, double** finalCovMatrix, int size)
void KruskalMST(Graph* graph, double** totalCovMatrix, Edge *result, int size)
{
  int V = graph->V;
  //Edge result[V];  // This will store the resultant MST
  int e = 0;  // An index variable, used for result[]
  int i = 0;  // An index variable, used for sorted edges

  bool isPresent = false;
  //printf("size is %d\n\n", size);

  int *availableRows = new int[size-1]();
  // start from row 1 till size
  for ( int i = 0 ; i < size-1 ; i++ ){
    availableRows[i] = i+1;
    //printf("i is %d\tavailableRows[i] is %d\n", i, availableRows[i]);
  }
  
  // Step 1:  Sort all the edges in non-decreasing order of their weight
  // If we are not allowed to change the given graph, we can create a copy of
  // array of edges
  qsort((*graph).edge, (*graph).E, sizeof((*graph).edge[0]), myComp);
  
  // Allocate memory for creating V ssubsets
  subset *subsets = new subset [V];
  
  // Create V subsets with single elements
  for (int v = 0; v < V; ++v)
    {
      subsets[v].parent = v;
      subsets[v].rank = 0;
    }
  
  // Number of edges to be taken is equal to V-1
  while (e < size - 1)
    {
      // Step 2: Pick the smallest edge. And increment the index
      // for next iteration
      Edge next_edge = graph->edge[i++];
      
      //printf("\nfind for x\n");
      int x = find(subsets, next_edge.src);

      //printf("\nfind for y\n");
      int y = find(subsets, next_edge.dest);
      
      // If including this edge does't cause cycle and the row was not used before, include it
      // in result and increment the index of result for next edge
      //printf("availableRows is\n");
      //for ( int i = 0 ; i < size-1 ; i++ ){
      //printf("%d  ",availableRows[i]);
      //}
      //printf("\nnext_edge.row is %d\nboolean value is %d\n\n",next_edge.row, isPresent);
      if (x != y && availableRows[next_edge.row-1] != 0)
        {
	  result[e++] = next_edge;
	  Union(subsets, x, y);
	  availableRows[next_edge.row-1] = 0;
	  //printf("edge.src is %d\tesde.dest s %d\nedge.row is %d\tedge.col is %d\nedgeweight is %.3f\n\n", next_edge.src, next_edge.dest, next_edge.src, next_edge.col, next_edge.weight);
	}
      // Else discard the next_edge
    }
  
  //for ( i = 0 ; i < e ; i++ ){
  //finalCovMatrix[result[i].row][result[i].col] =
  //totalCovMatrix[result[i].row][result[i].col];
  //}
  
  return;
}

//void maxSpanningTreeMatrix ( int size , double** totalCovMatrix , double** finalCovMatrix ) {
void maxSpanningTreeMatrix ( int size , double** totalCovMatrix , Edge *result ) {
  
  int i, j;   // simple counters
  
  int V = 2*size - 2;        // verticies
  int E = size*(size-1)/2;   // edges
  
  Graph* graph = createGraph(V,E);

  //*result = new Edge[V];
  
  int counter = 0;
  
  // initialize all edges in a graph
  for ( i = 1 ; i < size ; i++ ){
    for ( j = 0; j < i ; j++ ) {
      (*graph).edge[counter].src = i;
      (*graph).edge[counter].dest = size - 2 + j;

      (*graph).edge[counter].row = i;
      (*graph).edge[counter].col = j;
      
      (*graph).edge[counter].weight = fabs(totalCovMatrix[i][j]);
      counter++;
    }
  }
  
  //KruskalMST(graph, totalCovMatrix, finalCovMatrix, size);
  KruskalMST(graph, totalCovMatrix, result, size);
  
}

void ZarFOTT( double** corMat , int size , double** zarCovMat ){
  // take the maximum of each row and write to tempDouble
  // zarCovMat = tempDouble * tempDouble
  int temp_ind;
  double temp_value;
  double** tempDouble = make_2D_matrix(size, size, 0);
  for ( int i = 1; i < size ; i++ ){
    temp_value = 0;
    for ( int j = 0; j < i; j++ ){
      if ( corMat[i][j] > temp_value ){
	temp_value = corMat[i][j];
	temp_ind = j;
      }
    }
    tempDouble[i][temp_ind] = corMat[i][temp_ind];
  }

  //printf("tempDouble is\n");
  //for ( int i = 0 ; i < size ; i++ ){
  //for ( int j = 0; j < size ; j++ ){
  //printf("%.3f\t", tempDouble[i][j]);
  //} printf("\n");
  //}

  multiply_A_by_B( tempDouble , tempDouble , size, size, size, zarCovMat );
  

}


void svm_trueFalse_predict(const svm_model *model, const svm_node *x, double *est_value )
{
  if(model->param.svm_type == ONE_CLASS ||
     model->param.svm_type == EPSILON_SVR ||
     model->param.svm_type == NU_SVR)
    {
      printf("cannot calculate svm\n\n");
    }
  else
    {
      int i;
      int nr_class = model->nr_class;
      int l = model->l;
      
      double *kvalue = Malloc(double,l);
      for(i=0;i<l;i++)
	kvalue[i] = Kernel::k_function(x,model->SV[i],model->param);
      
      int *start = Malloc(int,nr_class);
      start[0] = 0;
      for(i=1;i<nr_class;i++)
	start[i] = start[i-1]+model->nSV[i-1];
      
      int *vote = Malloc(int,nr_class);
      for(i=0;i<nr_class;i++)
	vote[i] = 0;
      int p=0;
      for(i=0;i<nr_class;i++)
	for(int j=i+1;j<nr_class;j++)
	  {
	    double sum = 0;
	    int si = start[i];
	    int sj = start[j];
	    int ci = model->nSV[i];
	    int cj = model->nSV[j];
	    
	    int k;
	    double *coef1 = model->sv_coef[j-1];
	    double *coef2 = model->sv_coef[i];
	    for(k=0;k<ci;k++)
	      sum += coef1[si+k] * kvalue[si+k];
	    for(k=0;k<cj;k++)
	      sum += coef2[sj+k] * kvalue[sj+k];
	    sum -= model->rho[p++];
	    *est_value = sum;

	    if(sum > 0)
	      ++vote[i];
	    else
	      ++vote[j];
	  }


      //int vote_max_idx = 0;
      //for(i=1;i<nr_class;i++)
      //if(vote[i] > vote[vote_max_idx])
      //vote_max_idx = i;
      //
      //printf("printing the votes");
      //for ( int ii = 0 ; ii < nr_class ; ii++ ){
      //printf("%d\t", vote[ii]);
      //}
      //printf("\n\n");
      
      free(kvalue);
      free(start);
      free(vote);
    }
}
