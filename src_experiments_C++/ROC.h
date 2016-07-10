/*
 * ROC.h
 *
 *  Created on: March 17, 2017
 *      Author: Daniyar Bakir
 */

#include "standardHeaders.h"

struct linModel {
  double *a;
  double b;
};

struct trueFalse {
  int tp, tn, fp, fn;
};


int true_positive  (int* mtrue, int* mest, int n, int pos);
int true_negative  (int* mtrue, int* mest, int n, int neg);
int false_positive (int* mtrue, int* mest, int n, int pos);
int false_negative (int* mtrue, int* mest, int n, int pos);

int true_positive  (int mtrue, int mest, int pos);
int true_negative  (int mtrue, int mest, int neg);
int false_positive (int mtrue, int mest, int pos);
int false_negative (int mtrue, int mest, int pos);

trueFalse allResults( int* mtrue, int* mest, int N);

void svmTrueFalse (double** X, int* y, int N, int d, int* ind, svm_model* model, double* est_value, int* est_label );
//void svmTrueFalse (double** X, int* y, int N, int d, int* ind, svm_model* model, int* est_label );
void linClassTrueFalse (double** X, int* y, int N, int d, int* ind, linModel model, double* est_value, int* est_label);
void serdTrueFalse ( double** X, int* y, int N, int d, int* ind, linModel serd, int k, int m, double *est_value, int *est_label );

double trapezoid_area ( double x0, double x1, double y0, double y1 );
double calc_auc (double *f, int *mtrue, int *mest, int N, int tp_max, int fp_max, int neg, int pos );

double aucroc2 ( double* ev, int* tl, int N, int n0, int n1 );

int countLabels ( int* el, int label, int n );
int countValues ( double* el, int value, int n);

int countLabels2( double* eval, int* tlab, int n, int label, int w );

template <class T>
linModel convertLinModel (T orig_model, int d){

  // this function converts linear model such as DLDA, LDA, ZAR and others into a generic and ultimate linModel struct, so that a single linClassTst or something like it can be used easily.
  // orig_model   -- original model
  // d            -- number of features

  linModel new_model;

  new_model.a = new double[d]();

  for ( int i = 0 ; i < d ; i++ ){
    new_model.a[i] = orig_model.a[i];
  }

  new_model.b = orig_model.b;

  return new_model;

}

