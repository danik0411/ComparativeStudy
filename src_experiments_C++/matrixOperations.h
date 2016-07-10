/*
 * matrixOperations.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "standardHeaders.h"
#include "utilities.h"

#ifndef MATRIXOPERATIONS_H_
#define MATRIXOPERATIONS_H_

void diag_matrix(double** A, int n, double** diagA);

double determinant(double** A, int n);
void invert_matrix(double** A, int n, double** invA);
void chol_dc(double** A, int N, double** L);
void min_abs_a_b(double* a, double b, int n, double* min_value, int* min_index);
void bubble_sort(double* in, int* index, int length, int K);

void mean(double** X, int N, int d, double* m);
void mean(double* X, int N, int d, double* m);
void cov(double **X, int N, int d, double** S);

void sum_two_vector(double* a, double* b, int n, double factor_a, double factor_b, double factor_c, double* c);
void sum_two_matrix(double** A, double** B, int n1, int n2, double factor_A, double factor_B, double factor_C, double** C);

void multiply_a_by_scalar (double* a, double b, int n, double* c);

void multiply_a_by_b (double* a, double* b, int n, double* c);

void multiply_A_by_scalar (double** A, double b, int n1, int n2, double** B);
void multiply_A_by_b (double** A, double* b, int n1, int n2, double* c);
void multiply_A_by_B(double** A, double** B, int n1, int n2, int n3, double** C);

void transpose_2D_matrix( double** A, int d, double** B);

void copyA2B(double** A, int n, double** B);

void mean_min_dist(double** X, int* y, int N, int d, double* sig, int* ind);

double** make_2D_identity_matrix (int n1, int n2, double init0, double init1);

double** make_2D_matrix(int n1, int n2, double init);
void delete_2D_matrix(int n1, int n2, double** A);
double*** make_3D_dmatrix(int n1, int n2, int n3, double init);
void delete_3D_matrix(int n1, int n2, int n3, double*** A);

void read_2D_matrix_from_file(char* str, int n1, int n2, float** TT, FILE* ftp);

#endif /* MATRIXOPERATIONS_H_ */
