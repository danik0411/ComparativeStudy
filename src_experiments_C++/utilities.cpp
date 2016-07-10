/*
 * utilities.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "utilities.h"

/**  Error handler  **************************************************/

void erhand(char *err_msg)
/* Error handler */
{
	fprintf(stderr,"Run-time error:\n");
	fprintf(stderr,"%s\n", err_msg);
	fprintf(stderr,"Exiting to system.\n");
	exit(1);
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
