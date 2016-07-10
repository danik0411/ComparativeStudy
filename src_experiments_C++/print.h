/*
 * print.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef PRINT_H_
#define PRINT_H_

#include "standardHeaders.h"


void PrintDouble_A(double** data, int rows, int cols, char* name);
void PrintDouble_a(double* data, int cols, char* name);
void PrintInt_a(int* data, int cols, char* name);
void PrintDiscriminant( double *a, double b, int d, char* name);

void a_double_toScreen(double* a, int n);
void a_int_toScreen(int* a, int n);


#endif /* PRINT_H_ */
