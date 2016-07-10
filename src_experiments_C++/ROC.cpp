/*
 * classifiers.cpp
 *
 *  Created on: March 16, 2017
 *      Author: Daniyar Bakir
 */

#include "standardHeaders.h"
//#include "classifiers.cpp"
#include "matrixOperations.h"
#include "classifiers.h"
#include "print.h"
#include "random.h"
#include "ROC.h"

int true_positive ( int* mtrue, int* mest, int n, int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated
  // n     - number of samples in mtrue and mest
  // pos   - int that is considered as positive

  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){
    if ( mtrue[i] == pos && mest[i] == pos )
      a++;
    
  }

  return a;

}


int true_negative ( int* mtrue, int* mest, int n, int neg ) {

  // mtrue - vector of true
  // mest  - vector if estimated
  // n     - number of samples in mtrue and mest

  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){

    if ( mtrue[i] == neg && mest[i] == neg )
      a++;

  }

  return a;

}


int false_positive ( int* mtrue, int* mest, int n, int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated
  // n     - number of samples in mtrue and mest

  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){

    if ( mtrue[i] == (pos+1)%2 && mest[i] == pos )
      a++;

  }

  return a;

}


int false_negative ( int* mtrue, int* mest, int n, int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated
  // n     - number of samples in mtrue and mest

  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){

    if ( mtrue[i] == pos && mest[i] == (pos+1)%2 )
      a++;

  }

  return a;

}


int true_positive ( int mtrue, int mest, int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated

  if ( mtrue == pos && mest == pos )
    return 1;
  else
    return 0;

}

int true_negative ( int mtrue, int mest, int neg ) {

  // mtrue - vector of true
  // mest  - vector if estimated

  if ( mtrue == neg && mest == neg )
    return 1;
  else
    return 0;

}

int false_positive ( int mtrue, int mest, int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated

  if ( mtrue == (pos+1)%2 && mest == pos )
    return 1;
  else
    return 0;

}

int false_negative ( int mtrue, int mest , int pos ) {

  // mtrue - vector of true
  // mest  - vector if estimated

  if ( mtrue == pos && mest == (pos+1)%2 )
    return 1;
  else
    return 0;

}

trueFalse allResults ( int* mtrue, int* mest, int n ){

  // mtrue - vector of true
  // mest  - vector if estimated
  // n     - number of samples in mtrue and mest
  // res   - array of tp, tn, fp, fn

  trueFalse res;

  res.tp = 0;
  res.tn = 0;
  res.fp = 0;
  res.fn = 0;

  for ( int i = 0 ; i < n ; i++ ){


    if      ( mtrue[i] == 1 && mest[i] == 1 )
      res.tp++;
    else if ( mtrue[i] == 0 && mest[i] == 0 )
      res.tn++;
    else if ( mtrue[i] == 0 && mest[i] == 1 )
      res.fp++;
    else if ( mtrue[i] == 1 && mest[i] == 0 )
      res.fn++;

  }

  return res;

}


void svmTrueFalse (double** X, int* y, int N, int d, int* ind, svm_model* model, double* est_value, int* est_label ){
//void svmTrueFalse (double** X, int* y, int N, int d, int* ind, svm_model* model, int* est_label ){
  
  // this function will output the vector estimated values at a given threshold, but firstly I should find it somewhere in the code
  // X    -- input data
  // y    -- input label
  // N    -- amount of samples
  // d    -- data dimension
  // ind  -- best features
  // prior -- ???
  // model -- svm model
  
  int i, j;
  
  double er0 = 0;
  double er1 = 0;
  double error = 0.00;
  int N_0, N_1;

  //double *est_value = new double[N]();
  
  N_0 = 0;
  N_1 = 0;
  
  for (i=0; i<N; i++)
    if (y[i]==0)
      N_0++;
    else
      N_1++;
  
  svm_node *subdata_ts;
  svm_problem subcl_ts;

  //printf("svm initialization has been finished\n");
  
  //set data
  // be careful with this function since I changed it to the void instead of double
  if ((subdata_ts = Malloc(svm_node,(d+1)))==NULL)
    printf("could not create subdata_ts\n");
    //return 1;
  
  subcl_ts.l = 1;
  //allocate memory for data vector
  subcl_ts.x = Malloc(svm_node*, subcl_ts.l);
  subcl_ts.y = Malloc(double, subcl_ts.l);
  
  // set the pointer
  subcl_ts.x[0] = &subdata_ts[0];
  
  
  //set test data
  for (i=0; i<N; i++) {
    for(j=0; j<d; j++) {
      
      subdata_ts[j].value = X[i][ind[j]];
      subdata_ts[j].index = j+1;
    }
    subdata_ts[d].index = -1;
    subcl_ts.y[0] = y[i];
    
    //printf("svm predict is %d\n", (int) svm_predict(model, subcl_ts.x[0]));
    //printf("svm predict is %.3f\n", (double) svm_predict(model, subcl_ts.x[0]));
    
    est_label[i] = (int) svm_predict(model, subcl_ts.x[0]);
    //est_label[i] = (est_label[i]+1)%2;
    svm_trueFalse_predict(model, subcl_ts.x[0], &est_value[i]);

  }
  
  //free memory
  free(subcl_ts.x);
  free(subcl_ts.y);
  free(subdata_ts);

  //delete[] est_value;
  
}



/*
  void svm_trueFalse_predict(const svm_model *model, const svm_node *x, double* est_value )
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
  printf("sum is %.3f\tnr_class is %d\n", sum, model->nr_class);
  *est_value = sum;
  }
  
  free(kvalue);
  free(start);
  free(vote);
      }
  }
*/


void linClassTrueFalse (double** X, int* y, int N, int d, int* ind, linModel model, double* est_value, int* est_label) {
  
  // X         -- data
  // y         -- true class
  // N         -- number of samples
  // d         -- number of features
  // ind       -- best features
  // linModel  -- linear model of the classifier, e.g. MODEL.a and MODEL.b 
  // est_value -- estimated values
  // est_label -- estimated lables
  
  int i, j;
  double temp = 0.00;
  double* temp_X;
  double temp_error0 = 0.00;
  double temp_error1 = 0.00;
  double error = 0.00;
  int y_hat;
  int N_0, N_1;
  
  // uncomment next line if you want to compare with model.R
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
  for (i=0; i<N; i++) {
    for (j=0; j<d; j++) {
      temp_X[j] = X[i][ind[j]];
      //printf("%.3f\t", temp_X[j]);
    }
    //printf("\n\n");
    
    multiply_a_by_b (temp_X, model.a, d, &temp);

    est_value[i] = temp+model.b;
    //printf("est_value[%d] is %.3f\n", i, est_value[i]);
    //printf("est_value[%d] is %.3f\ttemp %.3f\n", i, est_value[i], temp);
    y_hat = ( est_value[i] <= 0) ? 0 : 1; 

    est_label[i] = y_hat;
    
    // uncomment next line if you want to comapre with model.R
    //values[i] = y_hat;
    
  }
  
  //printf("model error is %.3f\n", error);
  
  //uncomment this section you want to compare with model.R
  //PrintDouble( X, N, d, (char*)"terminal/MODEL_tstX.txt");
  //PrintInt( y, N, (char*)"terminal/MODEL_tsty.txt");
  //PrintInt( values, N , (char*)"terminal/MODEL_expy.txt");
  //printf("error is %f\n", error);
  //delete values;
  
  delete temp_X;

}


void serdTrueFalse ( double** X, int* y, int N, int d, int* ind, linModel serd, int k, int m, double *est_value, int *est_label )
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

		est_value[i] = temp+serd.b;
		
		est_label[i] = ( est_value[i] >= 0) ? 0 : 1; 

		//printf("est_value[%d] is %.3f\test_label[%d] is %d\n", i, est_value[i], i, est_label[i]);

		// uncomment next lines if you want to comapre with serd.R
		//expValues[i] = y_hat;
		//calValues[i] = temp+serd.b;

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

}

double calc_auc (double *f, int *mtrue, int *mest, int N, int N0, int N1 , int v0, int v1 ){

  // f      -- probability, they will be sordet in ascending order
  // mtrue  -- true label
  // mest   -- estimated label
  // N      -- number of samples
  // N0, N1 -- samples for label 0 and 1
  // v0, v1 -- label value for label 0 and label 1

  double *f_sorted = new double[N]();
  int *temp_ind = new int[N]();

  int *mtrue_sorted = new int[N]();
  int *mest_sorted  = new int[N]();

  int neg, pos; // label for positive and negative
  int tp_max, fp_max; // samples for positive and negative

  for ( int ii = 0 ; ii < N ; ii++ ){
    f_sorted[ii] = f[ii];
    temp_ind[ii] = ii;
  }

  bubble_sort( f_sorted, temp_ind, N, N );

  for ( int ii = 0 ; ii < N ; ii++ ){
    mtrue_sorted[ii] = mtrue[temp_ind[N-1-ii]];
    mest_sorted [ii] = mest [temp_ind[N-1-ii]];
    f_sorted    [ii] = f    [temp_ind[N-1-ii]];
  }

  if ( mest_sorted[0] == v0 ) { // if the negative eval[0] < 0 && elab[0] == 0
    neg = 0;
    pos = 1;
    fp_max = N0;
    tp_max = N1;
  } else {
    neg = 1;
    pos = 0;
    fp_max = N1;
    tp_max = N0;
  }
    
  /*

  printf("temp_ind is\n");
  for ( int ii = 0 ; ii < N ; ii++ ){
    printf("%d\t", temp_ind[ii]);
  }

  printf("\n\n");

  printf("fsorted is\n");
  for ( int ii = 0 ; ii < N ; ii++ ){
    printf("%.3f\t", f_sorted[ii]);
  }

  printf("\n\n");

  printf("mtrue_sorted is\n");
  for ( int ii = 0 ; ii < N ; ii++ ){
    printf("%d\t", mtrue_sorted[ii]);
  }

  printf("\n\n");

  printf("mest_sorted is\n");
  for ( int ii = 0 ; ii < N ; ii++ ){
    printf("%d\t", mest_sorted[ii]);
  }

  printf("\n\n");

  */

  

  int fp = 0, tp = 0;

  int fp0 = 0, tp0 = 0;


  double f0 = -DBL_MAX; 

  int i = 0;

  double A = 0.00;

  //printf("tp_max is %d\tfp_max is %d\n", tp_max, fp_max);

  while ( i < N ){

    if ( f_sorted[i] != f0 ){

      //printf("fp is %d\tfp0 is %d\ttp is %d\ttp0 is %d\ti is %d\t", fp, fp0, tp, tp0, i);
      //printf("mtrue is %d\tmest is %d\n", mtrue_sorted[i], mest_sorted[i]);

      //printf("fp/fpmax is %.3f\ntp");
      //A += trapezoid_area( (double) (fp/fp_max), (double) (fp0/fp_max) , (double) (tp/fp_max), (double) (tp0/tp_max) );
      //A += trapezoid_area( fp/fp_max, fp0/fp_max, tp/tp_max, tp0/tp_max);
      A += trapezoid_area( (double) fp, (double) fp0, (double) tp, (double) tp0);
      f0 = f_sorted[i];
      fp0 = fp;
      tp0 = tp;
    }

    //printf("true is %d\test is %d\ttp is %d\tfp is %d\n",mtrue_sorted[i]==pos,mest_sorted[i]==pos,tp,fp);

    if ( true_positive ( mtrue_sorted[i], mest_sorted[i], pos ) )
      tp++;
    else if ( false_positive ( mtrue_sorted[i], mest_sorted[i], pos ) )
      fp++;

    i++;

  }

  //A += trapezoid_area( (double) (fp_max/fp_max) , (double) (fp0/fp_max) , (double) (tp_max/tp_max) , (double) (tp0/tp_max) );
  //A += trapezoid_area( (fp_max/fp_max) , (fp0/fp_max) , (tp_max/tp_max) , (tp0/tp_max) );
  A += trapezoid_area( (double) fp_max , (double) fp0 , (double) tp_max , (double) tp0 );
  A /= (fp_max*tp_max);

  //PrintDouble_a( f_sorted, N, (char*)"eval_lsvm.txt");
  //PrintInt_a( mtrue_sorted, N, (char*)"true_labels.txt");
  //PrintInt_a( mest_sorted, N, (char*)"est_labels.txt");

  //printf("v.txt is\n");
  //for ( int ii = 0 ; ii < N ; ii++ )
  //printf( "ii is %d\tvalue is %.3f\test is %d\ttrue is %d\n", ii, f_sorted[ii], mest_sorted[ii], mtrue_sorted[ii] );

  //PrintDouble_a( f_sorted, N, (char*)"v.txt");
  //PrintInt_a( mtrue_sorted, N, (char*)"t.txt");
  //PrintInt_a( mest_sorted, N, (char*)"e.txt");

  delete[] f_sorted;
  delete[] temp_ind;
  delete[] mtrue_sorted;
  delete[] mest_sorted;

  return A;

}

double trapezoid_area ( double x0, double x1, double y0, double y1 ){

  return fabs(x1-x0) * (y0+y1)/2;

}

double aucroc2 ( double* ev, int* tl, int N, int n0, int n1 ){

  double *ev_sorted = new double[N]();
  int *temp_ind = new int[N]();

  int *tl_sorted = new int[N]();

  int* ind_pos;
  int* ind_neg;

  int neg, pos; // label for positive and negative
  int nn = 0, np = 0;

  int ipos = 0, ineg = 0;

  int zeros_in_neg = 0, zeros_in_pos = 0; // amount of samples of class 0 in eval < 0 and eval > 0
  int ones_in_neg = 0, ones_in_pos = 0;   // amount of samples of class 1 in eval < 0 and eval > 0

  for ( int ii = 0 ; ii < N ; ii++ ){
    ev_sorted[ii] = ev[ii];
    temp_ind[ii] = ii;
  }

  bubble_sort( ev_sorted, temp_ind, N, N );

  for ( int ii = 0 ; ii < N ; ii++ ){
    tl_sorted[ii] = tl[temp_ind[N-1-ii]];
    ev_sorted[ii] = ev[temp_ind[N-1-ii]];
  }

  zeros_in_neg = countLabels2 ( ev_sorted, tl_sorted, N, 0, -1 );
  zeros_in_pos = countLabels2 ( ev_sorted, tl_sorted, N, 0,  1 );

  ones_in_neg = countLabels2 ( ev_sorted, tl_sorted, N, 1, -1 );
  ones_in_pos = countLabels2 ( ev_sorted, tl_sorted, N, 1,  1 );

  neg = 0;
  pos = 1;
  nn = n0;
  np = n1;

  ind_pos = new int[np]();
  ind_neg = new int[nn]();

  ipos = 0;
  ineg = 0;

  for ( int i = 0 ; i < N ; i++ ){

    if ( tl_sorted[i] == neg ){
      ind_neg[ineg] = i;
      ineg++;
    } else if ( tl_sorted[i] == pos ){
      ind_pos[ipos] = i;
      ipos++;
    }
    
  }

  int isum = 0;
    
  for ( ipos = 0; ipos < np ; ipos++ ){
    for ( ineg = 0; ineg < nn; ineg++ ){

      //printf("ip=%d in=%d ep=%.3f en=%.3f\n",ipos,ineg,ev_sorted[ind_pos[ipos]],ev_sorted[ind_neg[ineg]]);
      if ( ev_sorted[ind_pos[ipos]] > ev_sorted[ind_neg[ineg]] )
	isum++;

    }
  }

  /*
  printf("neg is %d\tpos is %d\n",neg,pos);
  printf("ev_sorted is\n");
  for ( int i = 0 ; i < N ; i++ ) {
    printf("%.3f\t%d\n", ev_sorted[i], tl_sorted[i]);
  }
  */

  double auc = ( (double) isum ) / ( (double) (np*nn) );
  if ( auc < 0.5 )
    auc = 1-auc;

  
  delete[] ind_neg;
  delete[] ind_pos;
  delete[] tl_sorted;
  delete[] ev_sorted;
  delete[] temp_ind;

  return auc;

}


int countLabels ( int* el, int label, int n ){

  // count labels of "label" in el
  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){

    if ( el[i] == label )
      a++;
  }

  return a;

}




int countValues ( double* el, int value, int n){


  // count amount of positive or negative values
  // value = 1 for positive
  // value = -1 for negative
  int a = 0;

  for ( int i = 0 ; i < n ; i++ ){

    if ( value == 1 ){
      if ( el[i] >= 0 ) {
	a++;
      }
    } else if (value == -1){
      if ( el[i] <= 0 ) {
	a++;
      }
    } else {
      printf("only 1 and -1 can be used in countValues(...)\n");
    }

  }

  return a;

}

int countLabels2( double* eval, int* tlab, int n, int label, int w ){
  
  // this function counts how many values of "label" are in tlab when eval is pos
  // w is 1 or -1
  
  int a = 0;

  for ( int i = 0 ; i < n ; i++ ) {

    if ( w*eval[i] > 0 && tlab[i] == label )
      a++;
    
  }
  
  return a;
  
}

