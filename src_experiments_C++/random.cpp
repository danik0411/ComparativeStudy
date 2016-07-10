/*
 * random.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "random.h"

void UniformArray(double* a, int n, long* seed)
{
	int i;

	/* Returns array of doubles taken from uniform distribution on [0,1] */
	assert(n > 0);
	assert(a != 0); /* Calling program must allocate memory */

	for(i=0; i<n; i++)
		a[i] = (double)ran2(seed);

	return;
}

void GaussianArray(double* a, int n, double mu, double sigma, long* seed)
{
	/*
     Return array of doubles taken from N~(mu, sigma^2).
     Calls the Gaussian generator based on numerical recipes
	 */
	int i;

	assert(n > 0);
	assert(sigma >= 0); /* non-negative variance */
	assert(a != 0); /* Calling program must allocate memory */

	for(i=0; i<n; i++)
		a[i] = mu + sigma*(double)randn(seed);

	return;
}

float ran1(long* idum)
/*
	     Returns a uniform random deviate between 0.0 and 1.0
	     (exclusive of the endpoint values). Call with idum a negative
	     integer to initialize; thereafter, do not alter idum between
	     successive deviates in a sequence. RNMX should approximate
	     the largest floating value that is less than 1.
 */
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy)
	{ /* Initialize. */
		if (-(*idum) < 1)
			*idum = 1; /* Be sure to prevent idum = 0. */
		else
			*idum = -(*idum);

		for (j=NTAB+7; j>=0; j--)
		{ /* Load the shuffle table (after 8 warm-ups). */
			k = (*idum)/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;

			if (*idum < 0)
				*idum += IM;

			if (j < NTAB)
				iv[j] = *idum;
		}
		iy=iv[0];
	}

	k = (*idum)/IQ; /* Start here when not initializing. */
	*idum = IA*(*idum-k*IQ)-IR*k; /* Compute idum=(IA*idum) % IM without over */

	if (*idum < 0)
		*idum += IM; /* flows by Schrage’s method. */

	j = iy/NDIV; /* Will be in the range 0..NTAB-1. */
	iy = iv[j]; /* Output previously stored value and refill the */
	iv[j] = *idum; /* shuffle table. */

	if ((temp=AM*iy) > RNMX)
		return RNMX; /* Because users don’t expect endpoint values. */
	else
		return temp;
}

float ran2(long* idum)
/*
       Long period (> 2 × 10^18) random number generator of L’Ecuyer
       with Bays-Durham shuffle and added safeguards. Returns a
       uniform random deviate between 0.0 and 1.0 (exclusive of
       the endpoint values). Call with idum a negative integer to
       initialize; thereafter, do not alter idum between successive
       deviates in a sequence. RNMX should approximate the largest floating
       value that is less than 1.
 */

{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0)
	{ /* Initialize. */
		if (-(*idum) < 1)
			*idum=1; /* Be sure to prevent idum = 0. */
		else
			*idum = -(*idum);

		idum2=(*idum);
		for (j=NTAB+7; j>=0; j--)
		{ /* Load the shuffle table (after 8 warm-ups). */
			k = (*idum)/IQ1;
			*idum = IA1*(*idum-k*IQ1)-k*IR1;

			if (*idum < 0)
				*idum += IM1;

			if (j < NTAB)
				iv[j] = *idum;
		}
		iy=iv[0];
	}

	k = (*idum)/IQ1; /* Start here when not initializing. */
	*idum = IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */

	if (*idum < 0)
		*idum += IM1; /* overflows by Schrage’s method. */

	k = idum2/IQ2;
	idum2 = IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise. */

	if (idum2 < 0)
		idum2 += IM2;

	j = iy/NDIV1; /* Will be in the range 0..NTAB-1. */
	iy = iv[j]-idum2; /* Here idum is shuffled, idum and idum2 are */
	iv[j] = *idum; /* combined to generate output. */


	if (iy < 1)
		iy += IMM1;

	if ((temp=AM*iy) > RNMX)
		return (float)RNMX; /* Because users don’t expect endpoint values. */
	else
		return temp;
}

////////////////////////////////////////////////////////////
//random point generation, from NR

/*    double ran2(long *idum) - generates random numbers uniformily         */
/*                              distributed over (0,1).  Random sequence is */
/*                              initialized by calling ran2 with *idum a    */
/*                              negative number.                            */

float randn(long* idum)
{
	static int iset = 0;
	static double gset;
	float fac, rsq, v1, v2;

	if (*idum < 0)
		iset = 0; /* Reinitialize. */

	if (iset == 0)
	{
		do
		{
			v1 = 2.0*ran2(idum)-1.0;
			v2 = 2.0*ran2(idum)-1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;

		return v2*fac;
	}
	else
	{
		iset = 0;

		return gset;
	}
}

void shuffleArray(int N, long* seed, int* array)
{
	/*
     Return a scrambled version of array a. Note that a[] is changed.
     Knuth's Shuffle is used here.
	 */

	int i;
	int SwapIndex;
	int temp;


	for(i=N-1; i>0; i--)
	{
		SwapIndex = (int)(ran2(seed)*(i+1)); /* i+1 so that includes boundaries */
		temp = array[i];
		array[i] = array[SwapIndex];
		array[SwapIndex] = temp;
	}

	return;
}

double ncdf(double x)
{
	const double a[] = {3.16112374387056560e00,
			1.13864154151050156e02,
			3.77485237685302021e02,
			3.20937758913846947e03,
			1.85777706184603153e-1};
	const double b[] = {2.36012909523441209e01,
			2.44024637934444173e02,
			1.28261652607737228e03,
			2.84423683343917062e03};
	const double c[] = {5.64188496988670089e-1,
			8.88314979438837594e00,
			6.61191906371416295e01,
			2.98635138197400131e02,
			8.81952221241769090e02,
			1.71204761263407058e03,
			2.05107837782607147e03,
			1.23033935479799725e03,
			2.15311535474403846e-8};
	const double d[] = {1.57449261107098347e01,
			1.17693950891312499e02,
			5.37181101862009858e02,
			1.62138957456669019e03,
			3.29079923573345963e03,
			4.36261909014324716e03,
			3.43936767414372164e03,
			1.23033935480374942e03};
	const double p[] = {3.05326634961232344e-1,
			3.60344899949804439e-1,
			1.25781726111229246e-1,
			1.60837851487422766e-2,
			6.58749161529837803e-4,
			1.63153871373020978e-2};
	const double q[] = {2.56852019228982242e00,
			1.87295284992346047e00,
			5.27905102951428412e-1,
			6.05183413124413191e-2,
			2.33520497626869185e-3};
	const double xbreak = 0.46875;
	const double pi = acos(-1);
	double y,z, xnum, xden, r, del, e;
	int i;

	x = x/sqrt(2);
	y = fabs(x);

	/* evaluate erf for |x|<=0.46875 */

	if (fabs(x) <= xbreak)
	{
		z = y*y;
		xnum = a[4]*z;
		xden = z;

		for (i=0;i<3;i++)
		{
			xnum = (xnum+a[i])*z;
			xden = (xden+b[i])*z;
		}

		r = x*(xnum+a[3])/(xden+b[3]);
	}

	/* evaluate erf for 0.46875<=|x|<=4.0 */
	else if ((fabs(x) > xbreak) && (fabs(x) <= 4.0))
	{
		xnum = c[8]*y;
		xden = y;

		for (i=0; i<7; i++)
		{
			xnum = (xnum+c[i])*y;
			xden = (xden+d[i])*y;
		}

		r = (xnum+c[7])/(xden+d[7]);
		z = floor(y*16.0)/16.0;
		del = (y-z)*(y+z);
		r = exp(-z*z)*exp(-del)*r;

		if (x>0)
			r = 1-r;
		else
			r = r-1;
	}

	/* evaluate erf for |x|>4.0 */
	else
	{
		z = 1.0/(y*y);
		xnum = p[5]*z;
		xden = z;

		for (i=0; i<4; i++)
		{
			xnum = (xnum+p[i])*z;
			xden = (xden+q[i])*z;
		}

		r = z*(xnum+p[4])/(xden+q[4]);
		r = ((1/sqrt(pi))-r)/y;
		z = floor(y*16.0)/16.0;
		del = (y-z)*(y+z);
		r = exp(-z*z)*exp(-del)*r;

		if (x>0)
			r = 1-r;
		else
			r = r-1;
	}

	e = 0.5*(1+r);
	if (e>1)
		e=1;

	return e;
}
