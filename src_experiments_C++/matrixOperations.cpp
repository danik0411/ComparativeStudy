/*
 * matrixOperations.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "matrixOperations.h"


double determinant(double** A, int n)
{
	int i, j, j1, j2;
	double det = 0;
	double** m = NULL;
	double dinit = 0.00;

	if (n < 1)
	{ /* Error */

	}
	else if (n == 1)
	{ /* Shouldn't get used */
		det = A[0][0];
	}
	else if (n == 2)
	{
		det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
	}
	else
	{
		det = 0;
		for (j1=0; j1<n; j1++)
		{
			m = make_2D_matrix(n-1, n-1, dinit);

			for (i=1; i<n; i++)
			{
				j2 = 0;
				for (j=0; j<n; j++)
				{
					if (j == j1)
						continue;
					m[i-1][j2] = A[i][j];
					j2++;
				}
			}

			det += pow(-1.0,1.0+j1+1.0)*A[0][j1] * determinant(m, n-1);

			delete_2D_matrix(n-1, n-1, m);
		}
	}

	return det;
}

//=========================================//
//         Diagonal elements of            //
//	     a square matrix               //
//=========================================//
void diag_matrix(double** A, int n, double** diagA)
{

  int i, j;

  for ( i = 0 ; i < n ; i++ ){
    for ( j = 0; j < n ; j++ ){
      if ( i == j )
	diagA[i][j] = A[i][j];
      else
	diagA[i][j] = 0.00;
    }
  }

  return;
}


//=========================================//
/*         Invert a square matrix          */
//=========================================//
void invert_matrix(double** A, int n, double** invA)
{

	int* indxc = new int [n];
	int* indxr = new int [n];
	int* ipiv = new int [n];
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, temp;
	int ii, jj;
	double dinit = 0.00;

	double** a;
	a = make_2D_matrix (n, n, dinit);
	copyA2B(A, n, a);

	// actual Gauss-jordan algorithm
	for (j=0; j<n; j++)
		ipiv[j] = 0;
	for (i=0; i<n; i++)
	{
		big = 0.0;
		for (j=0; j<n; j++)
			if (ipiv[j] != 1)
				for (k=0; k<n; k++)
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big)
						{
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol)
		{
			for (l=0; l<n; l++)
				SWAP(a[irow][l], a[icol][l])
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)
			printf("\ngaussj: Singular Matrix\n");

		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0; l<n; l++)
			a[icol][l] *= pivinv;
		for (ll=0; ll<n; ll++)
			if (ll != icol)
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0; l<n; l++)
					a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=n-1; l>=0; l--)
	{
		if (indxr[l] != indxc[l])
			for (k=0; k<n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}

	copyA2B(a, n, invA);


	delete [] ipiv;
	delete [] indxr;
	delete [] indxc;
	delete_2D_matrix (n, n, a);

	return;
}

//=================================//
/*     Cholesky Decomposition    */
//=================================//
void chol_dc(double **A, int N, double **L)
/* Given N-by-N positive-definite symmetric matrix A, find its Cholesky decomposition A=L*L'
 * such that L is lower triangular matrix
 *
 *  Formula:
 *     L_ii = sqrt(A_ii - Sum_{k=1}^{i-1} L_ik^2)
 *     L_ji = (A_ij - Sum_{k=1}^{i-1} L_ik * L_jk)/L_ii, where j>i
 *
 * Reference: Numerical Recipe in C
 *     www.ma.utexas.edu/documentation/nr/bookcpdf.html
 * or  www.nr.com
 */
{
	int i, j, k;

	for (i=0; i<N; i++)
	{
		/* // above diagon are all zeros */
		for (j=0; j<i; j++)
			L[j][i] = 0.0;

		/* 		// compute lower diagonal elements */
		for (j=i; j<N; j++)
		{
			for (L[j][i]=A[i][j], k=i-1; k>=0; k--)
				L[j][i] -= L[i][k]*L[j][k];

			if (i==j)
			{ /* // diagonal elements */
				if (L[i][i] == 0)
					printf("Cholesky Decomposition failed. A must be positive-definite and symmetric.\n");
				L[i][i] = sqrt(L[i][i]);
			}
			else /* // non-diagonal elements */
				L[j][i] = L[j][i]/L[i][i];
		}
	}

	return;
}

void min_abs_a_b(double* a, double b, int n, double* min_value, int* min_index)
{
	/*
Purpose:
finds the minimum of abs(a(i)-b) value in the array data.
Inputs:
data: array of elements
n: size of the array
Output:
	 */
	double min;
	double threshold = 0.00;
	int index = 0;
	min = fabs(a[0]-b);
	for (int i=1; i<n; i++)
	{
		if (fabs(a[i]-b)<(min-threshold))
		{
			min = fabs(a[i]-b);
			index = i;
		}
	}
	*min_value =  min;
	*min_index = index;

	return;
}

//==============================================//
/*     Bubble-sort: Find K largest numbers    */
//==============================================//
void bubble_sort(double* in, int* index, int length, int K)
{
	int i,j;
	double temp1;
	int temp2;

	for (j=0; j<K; j++)
	{
		for (i=length-1; i>j; i--)
		{
			if (in[i] > in[i-1])
			{
				temp1 = in[i-1];
				in[i-1] = in[i];
				in[i] = temp1;
				temp2 = index[i-1];
				index[i-1] = index[i];
				index[i] = temp2;
			}
		}
	}

	return;
}
//row wise mean
void mean(double **X,
		int N,
		int d,
		double* m) 
{
	int i, j;

	for (i=0; i<d; i++)
	{
		m[i] = 0.00;
		for (j=0; j<N; j++)
			m[i] += (double) X[j][i];
		m[i] /= N;
	}

	return;
}

void mean(double *X,
		int N,
		int d,
		double* m)
{
	int i;
	*m = 0.00;
	for (i=0; i<d; i++)
		*m += (double) X[i];
	*m /= N;

	return;
}


void cov(double **X,
		int N,
		int d,
		double** S)
{
	int i, j, k;
	double* m = new double [d];
	double temp = 0.00;

	mean(X, N, d, m);

	for (j=0; j<d; j++)
		for (k=0; k<d; k++)
		{
			temp= 0.0;
			for (i=0; i<N; i++)
				temp += (X[i][j]-m[j])*(X[i][k]-m[k]);

			S[j][k] = temp/(N-1);
		}

	delete [] m;

	return;
}


void sum_two_vector(double* a, double* b, int n, double factor_a, double factor_b, double factor_c, double* c)
{
	for (int i=0;i<n;i++)
		c[i] = (factor_a*a[i] + factor_b*b[i])/factor_c;

	return;
}

void sum_two_matrix(double** A, double** B, int n1, int n2, double factor_A, double factor_B, double factor_C, double** C)
{

	for (int i=0;i<n1;i++)
		for (int j=0;j<n2;j++)
			C[i][j] = (factor_A*A[i][j] + factor_B*B[i][j])/factor_C;

	return;
}

void multiply_a_by_scalar (double* a, double b, int n, double* c)
{
	/*
Purpose:
Multiply a[n] by scalar and store it at c[n]
	 */
	int i;

	for (i=0;i<n;i++)
		c[i] = a[i]*b;

	return;
}

void multiply_a_by_b (double* a, double* b, int n, double* c)
{
	/*
Purpose:
Multiply a[n] by b[n] and store it at c
	 */
	int i;
	double sum = 0.00;

	for (i=0;i<n;i++)
		sum += a[i]*b[i];
	*c = sum;

	return;
}

void multiply_A_by_scalar (double** A, double b, int n1, int n2, double** B)
{
	/*
Purpose:
Multiply A[n1xn2] by b[1xn2] and store it at c[1xn1]
	 */
	int i, j;

	for (i=0; i<n1; i++)
	  for (j=0; j<n2; j++)
	    B[i][j] = A[i][j]*b;

	return;
}

void multiply_A_by_b (double** A, double* b, int n1, int n2, double* c)
{
	/*
Purpose:
Multiply A[n1xn2] by b[1xn2] and store it at c[1xn1]
	 */
	int i, j;
	double sum = 0.00;

	for (i=0; i<n1; i++)
	{
		sum = 0.00;
		for (j=0; j<n2; j++)
			sum += A[i][j]*b[j];
		c[i] = sum;
	}

	return;
}

void multiply_A_by_B (double** A, double** B, int n1, int n2, int n3, double** C)
{
	/*
Purpose:
Multiply A[n1xn2] by B[n2xn3] and store it at C[n1xn3]
	 */
	int i, j, k;
	double sum = 0.00;

	for (i=0; i<n1; i++)
	{
		for (j=0; j<n3; j++)
		{
			sum = 0;
			for (k=0; k<n2; k++)
				sum = A[i][k]*B[k][j] + sum;

			C [i][j] = sum;
		}
	}

	return;
}

void transpose_2D_matrix( double** A, int d, double** B){
  int i,j;

  for ( i = 0; i < d; i++ ){
    for (j = 0; j < d ; j++ ){
      B[j][i] = A[i][j];
    }
  }

  return;

}

void copyA2B(double** A, int n, double** B)
{
	int i, j;

	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			B[i][j] = A[i][j];
		}
	}

	return;
}


void mean_min_dist(double** X, int* y, int N, int d, double* sig, int* ind)
{
	/*
	 * function: mean_min_dist
	 *
	 * summary : computes the mean of
	 *   the minimum distances from each point
	 *   to the others
	 */

	int i, j, k, count;
	double* min = new double [N];
	double dis, dm[2], v;

	/* compute mimimum distances */
	for (i=0; i<N; i++)
	{
		dm[0] = HUGE_VAL;
		dm[1] = HUGE_VAL;

		if (y[i]==0)
		{
			for (j=0; j<i; j++)
				if (y[i]==y[j])
				{
					dis = 0;
					for (k=0; k<d; k++)
					{
						v = X[i][ind[k]]-X[j][ind[k]];
						dis += v*v;
					}
					if (dis<dm[0])
						dm[0] = dis;
				}

			for (j=i+1; j<N; j++)
				if (y[i]==y[j])
				{
					dis = 0;
					for (k=0; k<d; k++)
					{
						v = X[i][ind[k]]-X[j][ind[k]];
						dis += v*v;
					}
					if (dis<dm[0])
						dm[0] = dis;
				}

			min[i] = sqrt(dm[0]);
		}
		else
		{
			for (j=0; j<i; j++)
				if (y[i]==y[j])
				{
					dis = 0;
					for (k=0; k<d; k++)
					{
						v = X[i][ind[k]]-X[j][ind[k]];
						dis += v*v;
					}
					if (dis<dm[1])
						dm[1] = dis;
				}

			for (j=i+1; j<N; j++)
				if (y[i]==y[j])
				{
					dis = 0;
					for (k=0; k<d; k++)
					{
						v = X[i][ind[k]]-X[j][ind[k]];
						dis += v*v;
					}
					if (dis<dm[1])
						dm[1] = dis;
				}

			min[i] = sqrt(dm[1]);
		}
	}

	/* find mean */
	dm[0] = 0;
	dm[1] = 0;
	count = 0;
	for (i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			dm[0] += min[i];
			count++;
		}
		else
		{
			dm[1] += min[i];
		}
	}
	sig[0] = dm[0]/count;
	sig[1] = dm[1]/(N-count);

	delete min;

	return;
}

double** make_2D_matrix (int n1, int n2, double init)
										{
	int i;

	double** A = new double* [n1];

	for(i=0; i<n1; i++)
		A[i] = new double[n2];

	for (i=0;i<n1;i++)
	{
		for (int j=0;j<n2;j++)
			A[i][j] = init;
	}

	return A;
										}

double** make_2D_identity_matrix (int n1, int n2, double init0, double init1)
{
	int i;

	double** A = new double* [n1];   // array of pointers

	for(i=0; i<n1; i++)
	  A[i] = new double[n2];         // array of values given as addresses

	for (i=0;i<n1;i++)
	{
	  for (int j=0;j<n2;j++){
	      A[i][j] = init0;
	  }
	}

	for (i=0;i<n1;i++)
	  A[i][i] = init1;

	return A;
}

void delete_2D_matrix (int n1, int n2, double** A)
{
	int i;

	for(i=0; i<n1; i++)
	{
		delete [] A[i];//delete the rows
	}

	delete [] A;//delete the column
}


double*** make_3D_matrix (int n1, int n2, int n3, double init)
										{
	int i;

	double*** A = new double** [n1];

	for(i=0; i<n1; i++)
		A[i] = make_2D_matrix (n2, n3, init);

	return A;
										}

void delete_3D_matrix (int n1, int n2, int n3, double*** A)
{
	int i;

	for(i=0; i<n1; i++)
		delete_2D_matrix (n2, n3, A[i]);

	delete [] A;
}
