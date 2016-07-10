/*
 * print.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "print.h"

void PrintDouble_A(double** data, int rows, int cols, char* name)
{
	int i, j;
	FILE *fp;

	fp = fopen(name, "w");
	if(fp==NULL)
	{
		printf("cannot open output file\n");
		return;
	}

	for (i=0; i<rows; i++)
	{
		for (j=0; j<cols; j++)
		{
		  //fprintf(fp, "%.3f & ", data[i][j]);
			fprintf(fp, "%.3f\t", data[i][j]);
		}
		//fprintf(fp, "\\\\\n");
		fprintf(fp, "\n");
	}

	fclose(fp);

	return;
}

void PrintDouble_a(double* data, int cols, char* name)
{
	int j;
	FILE *fp;

	fp = fopen(name, "w");
	if(fp==NULL)
	{
		printf("cannot open output file\n");
		return;
	}
	for (j=0; j<cols; j++)
	{
	  //printf("%.6f\n", data[j]);
	  //fprintf(fp, "%.6f\n", data[j]);
	  fprintf(fp, "%.6f\n", data[j]);
	}

	fclose(fp);

}

void PrintInt_a(int* data, int cols, char* name)
{
	int j;
	FILE *fp;

	fp = fopen(name, "w");
	if(fp==NULL)
	{
		printf("cannot open output file\n");
		return;
	}
	for (j=0; j<cols; j++)
	{
	  fprintf(fp, "%d\n", data[j]);
	}
	//fprintf(fp, "\n");
	fclose(fp);

}

void PrintDiscriminant( double *a, double b, int d, char* name){
	int j;
	FILE *fp;

	fp = fopen(name, "w");
	if(fp==NULL)
	{
		printf("cannot open output file\n");
		return;
	}
	fprintf(fp, "%.6f\n\n", b);
	for (j=0; j<d; j++)
	{
	  fprintf(fp, "%.6f", a[j]);
	  if ( j < d-1 )
	    fprintf(fp,"\t");
	}
	fprintf(fp, "\n");

	fclose(fp);

	return;

}

void a_double_toScreen(double* a, int n){


  for ( int i = 0 ; i < n ; i++ )
    printf("%.3f\n", a[i]);

  return;

}

void a_int_toScreen(int* a, int n){


  for ( int i = 0 ; i < n ; i++ )
    printf("%d\n", a[i]);

  return;

}
