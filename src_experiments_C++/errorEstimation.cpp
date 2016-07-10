/*
 * errorEstimation.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#include "errorEstimation.h"

double ldaBolster(double** X, int* y, int N, int d, int* ind, LDA lda)
{
	int i, j, k;
	double temp_error = 0.00;
	double error = 0.00;
	double x, norm = 0.00, v;
	double* sig = new double [2];
	double cp[] = {1,0.67448425,1.17741394,1.53816223,1.83213806,2.08601379, \
			2.31260045,2.51908941,2.710005,2.88839621,3.05643871, \
			3.21574102,3.36753874,3.51279834,3.65229692,3.78666860, \
			3.91643959,4.04205182,4.16388078,4.28224868,4.39743449, \
			4.50968163,4.61920387,4.72619084,4.83081033,4.93321291,\
			5.03356934,5.13183594,5.22827148,5.32318115,5.41625977,\
			5.50781250,5.59783936,5.68634033,5.77362061,5.85968018,\
			5.94421387,6.02783203,6.11022949,6.19140625,6.27166748,\
			6.35101318,6.42913818,6.50665283,6.58294678,6.65832520,\
			6.73309326,6.80694580,6.87988281,6.95220947,7.02392578,\
			7.09472656,7.16491699,7.23419189,7.30316162,7.37121582,\
			7.43865967,7.50549316,7.57202148,7.63763428,7.70294189,\
			7.76763916,7.83172607,7.89520264,7.95837402,8.02093506,\
			8.08288574,8.14453125,8.20587158,8.26660156,8.32672119};


	mean_min_dist(X, y, N, d, sig, ind);

	for (i=0; i<d; i++)
		norm += lda.a[i]*lda.a[i];
	norm = sqrt(norm);

	for(i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			v = lda.b;

			for (j=0; j<d; j++)
				v += X[i][ind[j]]*lda.a[j];

			x = fabs(v)/(norm*sig[0]/cp[d]);

			if (v>0)
				temp_error = temp_error + ncdf(x); /* wrong side */
			else
				temp_error = temp_error + 1 - ncdf(x); /* right side */
		}
		else
		{
			v = lda.b;

			for (j=0; j<d; j++)
				v += X[i][ind[j]]*lda.a[j];

			x = fabs(v)/(norm*sig[1]/cp[d]);

			if (v<0)
				temp_error = temp_error + ncdf(x); /* wrong side */
			else
				temp_error = temp_error + 1 - ncdf(x); /* right side */
		}
	}

	error = temp_error/N;

	delete sig;

	return error;
}


//this function is not yet working
double qdaBolster(double** X, int* y, int N, int d, int* ind, LDA lda)
{
	int i, j, k;
	double temp_error = 0.00;
	double error = 0.00;
	double x, norm = 0.00, v;
	double* sig = new double [2];
	double cp[] = {1,0.67448425,1.17741394,1.53816223,1.83213806,2.08601379, \
		2.31260045,2.51908941,2.710005,2.88839621,3.05643871, \
		3.21574102,3.36753874,3.51279834,3.65229692,3.78666860, \
		3.91643959,4.04205182,4.16388078,4.28224868,4.39743449, \
		4.50968163,4.61920387,4.72619084,4.83081033,4.93321291,\
		5.03356934,5.13183594,5.22827148,5.32318115,5.41625977,\
		5.50781250,5.59783936,5.68634033,5.77362061,5.85968018,\
		5.94421387,6.02783203,6.11022949,6.19140625,6.27166748,\
		6.35101318,6.42913818,6.50665283,6.58294678,6.65832520,\
		6.73309326,6.80694580,6.87988281,6.95220947,7.02392578,\
		7.09472656,7.16491699,7.23419189,7.30316162,7.37121582,\
		7.43865967,7.50549316,7.57202148,7.63763428,7.70294189,\
		7.76763916,7.83172607,7.89520264,7.95837402,8.02093506,\
	8.08288574,8.14453125,8.20587158,8.26660156,8.32672119};
	
	
	mean_min_dist(X, y, N, d, sig, ind);
	
	for (i=0; i<d; i++)
		norm += lda.a[i]*lda.a[i];
	norm = sqrt(norm);
	
	for(i=0; i<N; i++)
	{
		if (y[i]==0)
		{
			v = lda.b;
			
			for (j=0; j<d; j++)
				v += X[i][ind[j]]*lda.a[j];
			
			x = fabs(v)/(norm*sig[0]/cp[d]);
			
			if (v>0)
				temp_error = temp_error + ncdf(x); /* wrong side */
			else
				temp_error = temp_error + 1 - ncdf(x); /* right side */
		}
		else
		{
			v = lda.b;
			
			for (j=0; j<d; j++)
				v += X[i][ind[j]]*lda.a[j];
			
			x = fabs(v)/(norm*sig[1]/cp[d]);
			
			if (v<0)
				temp_error = temp_error + ncdf(x); /* wrong side */
			else
				temp_error = temp_error + 1 - ncdf(x); /* right side */
		}
	}
	
	error = temp_error/N;
	
	delete sig;
	
	return error;
}


double lsvmBolster(double** X, int *y, int N, int d, int* ind, int N_bolster, long* seed)
{
	int i, j, k;
	double temp_error = 0.00;
	double error = 0.00;
	double dinit = 0.00;
	double** X_bolster = make_2D_matrix(1, d, dinit);
	double* sig = new double [2];
	int* ind_temp;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	double cp[] = {1,0.67448425,1.17741394,1.53816223,1.83213806,2.08601379, \
			2.31260045,2.51908941,2.710005,2.88839621,3.05643871, \
			3.21574102,3.36753874,3.51279834,3.65229692,3.78666860, \
			3.91643959,4.04205182,4.16388078,4.28224868,4.39743449, \
			4.50968163,4.61920387,4.72619084,4.83081033,4.93321291,\
			5.03356934,5.13183594,5.22827148,5.32318115,5.41625977,\
			5.50781250,5.59783936,5.68634033,5.77362061,5.85968018,\
			5.94421387,6.02783203,6.11022949,6.19140625,6.27166748,\
			6.35101318,6.42913818,6.50665283,6.58294678,6.65832520,\
			6.73309326,6.80694580,6.87988281,6.95220947,7.02392578,\
			7.09472656,7.16491699,7.23419189,7.30316162,7.37121582,\
			7.43865967,7.50549316,7.57202148,7.63763428,7.70294189,\
			7.76763916,7.83172607,7.89520264,7.95837402,8.02093506,\
			8.08288574,8.14453125,8.20587158,8.26660156,8.32672119};

	mean_min_dist(X, y, N, d, sig, ind);

	ind_temp = new int [d];
	for (i=0; i<d; i++)
		ind_temp[i] = i;

	svm = svmTrn(X, y, N, d, ind, 0, &subdata_svm, &subcl_svm);

	//generate points
	for (i=0; i<N; i++)
		for (j=0; j<N_bolster; j++)
		{
			for (k=0; k<d; k++)
				if (y[i]==0)
					X_bolster[0][k] = X[i][ind[k]] + randn(seed)*sig[0]/cp[d];
				else
					X_bolster[0][k] = X[i][ind[k]] + randn(seed)*sig[1]/cp[d];

			temp_error += svmTst(X_bolster, y+i, 1, d, ind_temp,  2.00, svm);
		}

	error = temp_error/N_bolster/N;

	delete ind_temp;
	delete sig;
	delete_2D_matrix(1, d, X_bolster);
	svmDestroy(svm, subdata_svm, &subcl_svm);

	return error;
}

double ksvmBolster(double** X, int *y, int N, int d, int* ind, int N_bolster, long* seed)
{
	int i, j, k;
	double temp_error = 0.00;
	double error = 0.00;
	double dinit = 0.00;
	double** X_bolster = make_2D_matrix(1, d, dinit);
	double* sig = new double [2];
	int* ind_temp;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	double cp[] = {1,0.67448425,1.17741394,1.53816223,1.83213806,2.08601379, \
			2.31260045,2.51908941,2.710005,2.88839621,3.05643871, \
			3.21574102,3.36753874,3.51279834,3.65229692,3.78666860, \
			3.91643959,4.04205182,4.16388078,4.28224868,4.39743449, \
			4.50968163,4.61920387,4.72619084,4.83081033,4.93321291,\
			5.03356934,5.13183594,5.22827148,5.32318115,5.41625977,\
			5.50781250,5.59783936,5.68634033,5.77362061,5.85968018,\
			5.94421387,6.02783203,6.11022949,6.19140625,6.27166748,\
			6.35101318,6.42913818,6.50665283,6.58294678,6.65832520,\
			6.73309326,6.80694580,6.87988281,6.95220947,7.02392578,\
			7.09472656,7.16491699,7.23419189,7.30316162,7.37121582,\
			7.43865967,7.50549316,7.57202148,7.63763428,7.70294189,\
			7.76763916,7.83172607,7.89520264,7.95837402,8.02093506,\
			8.08288574,8.14453125,8.20587158,8.26660156,8.32672119};

	mean_min_dist(X, y, N, d, sig, ind);

	ind_temp = new int [d];
	for (i=0; i<d; i++)
		ind_temp[i] = i;

	svm = svmTrn(X, y, N, d, ind, 2, &subdata_svm, &subcl_svm);

	//generate points
	for (i=0; i<N; i++)
		for (j=0; j<N_bolster; j++)
		{
			for (k=0; k<d; k++)
				if (y[i]==0)
					X_bolster[0][k] = X[i][ind[k]] + randn(seed)*sig[0]/cp[d];
				else
					X_bolster[0][k] = X[i][ind[k]] + randn(seed)*sig[1]/cp[d];

			temp_error += svmTst(X_bolster, y+i, 1, d, ind_temp, 2.00, svm);
		}

	error = temp_error/N_bolster/N;

	delete ind_temp;
	delete sig;
	delete_2D_matrix(1, d, X_bolster);
	svmDestroy(svm, subdata_svm, &subcl_svm);

	return error;
}

double ldaLOO(double** X, int* y, int N, int D, int d, double prior)
{
	int i, j, k;
	double dinit = 0.00;
	double** X_trn = make_2D_matrix(N-1, D, dinit);
	int* y_trn = new int [N-1];
	double** X_tst = make_2D_matrix(1, D, dinit);
	double temp_error = 0.00;
	double error = 0.00;
	int* best_features;

	LDA lda;

	for (i=0; i<N; i++)
	{
		for (k=0; k<D; k++)
			X_tst[0][k] = X[i][k];

		for (j=0; j<i; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j][k] = X[j][k];
			y_trn[j] = y[j];
		}
		for (j=i+1; j<N; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j-1][k] = X[j][k];
			y_trn[j-1] = y[j];
		}

		best_features = new int [d];
		featureSelection(X_trn, y_trn, N-1, D, d, best_features);

		lda.a = new double [d];
		ldaTrn(X_trn, y_trn, N-1, d, best_features, prior, &lda);
		temp_error += ldaTst(X_tst, y+i, 1, d, best_features,  prior, lda);

		delete lda.a;
		delete best_features;
	}

	error = temp_error/N;

	delete_2D_matrix(N-1, D, X_trn);
	delete y_trn;
	delete_2D_matrix(1, D, X_tst);

	return error;
}

double qdaLOO(double** X, int* y, int N, int D, int d, double prior)
{
	int i, j, k;
	double dinit = 0.00;
	double** X_trn = make_2D_matrix(N-1, D, dinit);
	int* y_trn = new int [N-1];
	double** X_tst = make_2D_matrix(1, D, dinit);
	double temp_error = 0.00;
	double error = 0.00;
	int* best_features;
	
	QDA qda;
	
	for (i=0; i<N; i++)
	{
		for (k=0; k<D; k++)
			X_tst[0][k] = X[i][k];
		
		for (j=0; j<i; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j][k] = X[j][k];
			y_trn[j] = y[j];
		}
		for (j=i+1; j<N; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j-1][k] = X[j][k];
			y_trn[j-1] = y[j];
		}
		
		best_features = new int [d];
		featureSelection(X_trn, y_trn, N-1, D, d, best_features);
		
		qda.a = make_2D_matrix(d, d, dinit);
		qda.b = new double [d];
		qdaTrn(X_trn, y_trn, N-1, d, best_features, prior, &qda);
		temp_error += qdaTst(X_tst, y+i, 1, d, best_features,  prior, qda);
		
		delete qda.b;
		delete_2D_matrix(d, d, qda.a);
		delete best_features;
	}
	
	error = temp_error/N;
	
	delete_2D_matrix(N-1, D, X_trn);
	delete y_trn;
	delete_2D_matrix(1, D, X_tst);
	
	return error;
}


double lsvmLOO(double** X, int* y, int N, int D, int d)
{
	int i, j, k;
	double dinit = 0.00;
	double** X_trn = make_2D_matrix(N-1, D, dinit);
	int* y_trn = new int [N-1];
	double** X_tst = make_2D_matrix(1, D, dinit);
	double temp_error = 0.00;
	double error = 0.00;
	int* best_features;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	for (i=0; i<N; i++)
	{
		for (k=0; k<D; k++)
			X_tst[0][k] = X[i][k];

		for (j=0; j<i; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j][k] = X[j][k];
			y_trn[j] = y[j];
		}
		for (j=i+1; j<N; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j-1][k] = X[j][k];
			y_trn[j-1] = y[j];
		}

		best_features = new int [d];
		featureSelection(X_trn, y_trn, N-1, D, d, best_features);

		svm = svmTrn(X_trn, y_trn, N-1, d, best_features, 0, &subdata_svm, &subcl_svm);
		temp_error += svmTst(X_tst, y+i, 1, d, best_features,  2.00, svm);

		svmDestroy(svm, subdata_svm, &subcl_svm);
		delete best_features;
	}

	error = temp_error/N;

	delete_2D_matrix(N-1, D, X_trn);
	delete y_trn;
	delete_2D_matrix(1, D, X_tst);

	return error;
}

double ksvmLOO(double** X, int* y, int N, int D, int d)
{
	int i, j, k;
	double dinit = 0.00;
	double** X_trn = make_2D_matrix(N-1, D, dinit);
	int* y_trn = new int [N-1];
	double** X_tst = make_2D_matrix(1, D, dinit);
	double temp_error = 0.00;
	double error = 0.00;
	int* best_features;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	for (i=0; i<N; i++)
	{
		for (k=0; k<D; k++)
			X_tst[0][k] = X[i][k];

		for (j=0; j<i; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j][k] = X[j][k];
			y_trn[j] = y[j];
		}
		for (j=i+1; j<N; j++)
		{
			for (k=0; k<D; k++)
				X_trn[j-1][k] = X[j][k];
			y_trn[j-1] = y[j];
		}

		best_features = new int [d];
		featureSelection(X_trn, y_trn, N-1, D, d, best_features);

		svm = svmTrn(X_trn, y_trn, N-1, d, best_features, 2, &subdata_svm, &subcl_svm);
		temp_error += svmTst(X_tst, y+i, 1, d, best_features,  2.00, svm);

		svmDestroy(svm, subdata_svm, &subcl_svm);
		delete best_features;
	}

	error = temp_error/N;

	delete_2D_matrix(N-1, D, X_trn);
	delete y_trn;
	delete_2D_matrix(1, D, X_tst);

	return error;
}

double ldaCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r;
	int l_0, l_1;
	int N_0, N_1;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double temp_error00 = 0.00;
	double temp_error11 = 0.00;
	double error = 0.00;

	double** X_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_trn;
	int* y_tst0;
	int* y_tst1;
	
	int l00 = 0;
	int l11 = 0;
	double err0 = 0.00;
	double err1 = 0.00;

	int* shuffled_0;
	int* shuffled_1;

	int* best_features;

	LDA lda;

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;


	/* size of sample to leave out */
	l_0 = int(floor(N_0/K));
	l_1 = int(floor(N_1/K));

	X_trn = make_2D_matrix(N-l_0-l_1, D, dinit);
	y_trn = new int [N-l_0-l_1];
	X_tst0 = make_2D_matrix(l_0, D, dinit);
	X_tst1 = make_2D_matrix(l_1, D, dinit);	
	y_tst0 = new int [l_0];
	y_tst1 = new int [l_1];

	shuffled_0 = new int [N_0];
	shuffled_1 = new int [N_1];

	for (i=0; i<N_0; i++)
		shuffled_0[i] = i;

	for (i=0; i<N_1; i++)
		shuffled_1[i] = i;

	for (r=0; r<R; r++)
	{
		shuffleArray(N_0, seed, shuffled_0);
		shuffleArray(N_1, seed, shuffled_1);

		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_0[i]][j];

				y_trn[i] = y[shuffled_0[i]];
			}

			for (i=k*l_0; i<(k+1)*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_tst0[i-k*l_0][j] = X[shuffled_0[i]][j];

				y_tst0[i-k*l_0] = y[shuffled_0[i]];
				l00++;
			}

			for (i=(k+1)*l_0; i<N_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_0][j] = X[shuffled_0[i]][j];

				y_trn[i-l_0] = y[shuffled_0[i]];
			}
			

			for (i=0; i<k*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i][j] = X[N_0+shuffled_1[i]][j];

				y_trn[N_0-l_0+i] = y[N_0+shuffled_1[i]];
			}

			for (i=k*l_1; i<(k+1)*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_tst1[i-k*l_1][j] = X[N_0+shuffled_1[i]][j];

				y_tst1[i-k*l_1] = y[N_0+shuffled_1[i]];
				l11++;
			}

			for (i=(k+1)*l_1; i<N_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i-l_1][j] = X[N_0+shuffled_1[i]][j];

				y_trn[N_0-l_0+i-l_1] = y[N_0+shuffled_1[i]];
			}

			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_0-l_1, D, d, best_features);

			lda.a = new double [d];
			ldaTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, prior, &lda);
			//printf("Here\n");
			//temp_error0 += ldaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, lda);//pior is set to 2.00 so it simply devides by total samples
			//temp_error1 += ldaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, lda);
			err0 = ldaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, lda);
			err0=(double)l_0*err0;
			temp_error00 += err0;
			err1 = ldaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, lda);
			err1=(double)l_1*err1;
			temp_error11 += err1;

			delete lda.a;
			delete best_features;
			//printf("N_0=%d\n",N_0);
			//printf("N_1=%d\n",N_1);
			//printf("N=%d\n",N);
		}
	}
	

	 
	if (prior==2.00)
	{
		
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error00/l00+((double)N_1/(double)N)*temp_error11/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error00/l00+(prior)*temp_error11/l11;
	}
	
	//printf("errorSep2=%.2f\n",error);

	delete_2D_matrix(N-l_0-l_1, D, X_trn);
	delete_2D_matrix(l_0, D, X_tst0);
	delete_2D_matrix(l_1, D, X_tst1);
	delete y_trn;
	delete y_tst0;
	delete y_tst1;
	delete shuffled_0;
	delete shuffled_1;

	return error;
}



//This is the Mixture CV
double ldaCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r, i0, i1;
	int l_01;
	int l_0, l_1;
	double temp_error0 = 0.00;
	double err0 = 0.00;
	double temp_error1 = 0.00;
	double err1 = 0.00;
	double error = 0.00;
	int l00 = 0;
	int l11 = 0;
	
	double** X_trn;
	int* y_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_tst0;
	int* y_tst1;
	
	int* shuffled_01;
	
	int* best_features;
	int N_0, N_1;
	
	LDA lda;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_01 = int(floor(N/K));
	//l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_01, D, dinit);
	y_trn = new int [N-l_01];
	
	shuffled_01 = new int [N];
	//shuffled_1 = new int [N_1];
	
	for (i=0; i<N; i++)
		shuffled_01[i] = i;
	
	//for (i=0; i<N_1; i++)
		//shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N, seed, shuffled_01);
		//shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_01; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_01[i]][j];
				
				y_trn[i] = y[shuffled_01[i]];
			}
			
			l_0 = 0;
			l_1 = 0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
				if (y[shuffled_01[i]]==0)
					l_0++;
				else
					l_1++;
			}
			
			if (l_0>0) {
			X_tst0 = make_2D_matrix(l_0, D, dinit);
			y_tst0 = new int [l_0];
			}
			if (l_1>0) {
			X_tst1 = make_2D_matrix(l_1, D, dinit);	
			y_tst1 = new int [l_1];
			}
			
			i0=0;
			i1=0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
	            
				if (y[shuffled_01[i]]==0) {
					for (j=0; j<D; j++)
						X_tst0[i0][j] = X[shuffled_01[i]][j];
				
					y_tst0[i0] = y[shuffled_01[i]];
					l00++;
					i0++;
				} else {
					for (j=0; j<D; j++)
						X_tst1[i1][j] = X[shuffled_01[i]][j];
					
					y_tst1[i1] = y[shuffled_01[i]];
					l11++;
					i1++;
				}
			}
			
			for (i=(k+1)*l_01; i<N; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_01][j] = X[shuffled_01[i]][j];
				
				y_trn[i-l_01] = y[shuffled_01[i]];
			}
			
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_01, D, d, best_features);
			
			lda.a = new double [d];
			ldaTrn(X_trn, y_trn, N-l_01, d, best_features, prior, &lda);
			if (l_0>0) {
			 err0 = ldaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, lda);
				err0=(double)l_0*err0;
				temp_error0 += err0;
			}
			if (l_1>0) {
			err1 = ldaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, lda);
				err1=(double)l_1*err1;
				temp_error1 += err1;
			}

			
			delete lda.a;
			delete best_features;
			if (l_0>0) {
			delete_2D_matrix(l_0, D, X_tst0);
			delete y_tst0;
			}
			if (l_1>0) {
			delete_2D_matrix(l_1, D, X_tst1);
			delete y_tst1;
			}
			/*
			printf("l_0=%d\n",l_0);
			printf("l_1=%d\n",l_1);
			printf("i0=%d\n",i0);
			printf("i1=%d\n",i1);
			printf("++++++\n",i1);
			*/
			
		}
	}
	
	if (prior==2.00)
	{
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error0/l00+((double)N_1/(double)N)*temp_error1/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error0/l00+(prior)*temp_error1/l11;
	}
	// printf("errorMix=%.2f\n",error);
	
	delete_2D_matrix(N-l_01, D, X_trn);
	delete y_trn;
	delete shuffled_01;
	//delete shuffled_1;
	
	return error;
}



/////////////////////
double qdaCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r;
	int l_0, l_1;
	int N_0, N_1;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double temp_error00 = 0.00;
	double temp_error11 = 0.00;
	double error = 0.00;
	
	int l00 = 0;
	int l11 = 0;
	double err0 = 0.00;
	double err1 = 0.00;
	
	double** X_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_trn;
	int* y_tst0;
	int* y_tst1;
	
	int* shuffled_0;
	int* shuffled_1;
	
	int* best_features;	
	
	QDA qda;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_0 = int(floor(N_0/K));
	l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_0-l_1, D, dinit);
	X_tst0 = make_2D_matrix(l_0, D, dinit);
	X_tst1 = make_2D_matrix(l_1, D, dinit);
	y_trn = new int [N-l_0-l_1];
	y_tst0 = new int [l_0];
	y_tst1 = new int [l_1];
	
	shuffled_0 = new int [N_0];
	shuffled_1 = new int [N_1];
	
	for (i=0; i<N_0; i++)
		shuffled_0[i] = i;
	
	for (i=0; i<N_1; i++)
		shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N_0, seed, shuffled_0);
		shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
				for (i=0; i<k*l_0; i++)
				{
					for (j=0; j<D; j++)
						X_trn[i][j] = X[shuffled_0[i]][j];
					
					y_trn[i] = y[shuffled_0[i]];
				}
				
				for (i=k*l_0; i<(k+1)*l_0; i++)
				{
					for (j=0; j<D; j++)
						X_tst0[i-k*l_0][j] = X[shuffled_0[i]][j];
					
					y_tst0[i-k*l_0] = y[shuffled_0[i]];
					l00++;
				}
				
				for (i=(k+1)*l_0; i<N_0; i++)
				{
					for (j=0; j<D; j++)
						X_trn[i-l_0][j] = X[shuffled_0[i]][j];
					
					y_trn[i-l_0] = y[shuffled_0[i]];
				}
				
				
				for (i=0; i<k*l_1; i++)
				{
					for (j=0; j<D; j++)
						X_trn[N_0-l_0+i][j] = X[N_0+shuffled_1[i]][j];
					
					y_trn[N_0-l_0+i] = y[N_0+shuffled_1[i]];
				}
				
				for (i=k*l_1; i<(k+1)*l_1; i++)
				{
					for (j=0; j<D; j++)
						X_tst1[i-k*l_1][j] = X[N_0+shuffled_1[i]][j];
					
					y_tst1[i-k*l_1] = y[N_0+shuffled_1[i]];
					l11++;
				}
				
				for (i=(k+1)*l_1; i<N_1; i++)
				{
					for (j=0; j<D; j++)
						X_trn[N_0-l_0+i-l_1][j] = X[N_0+shuffled_1[i]][j];
					
					y_trn[N_0-l_0+i-l_1] = y[N_0+shuffled_1[i]];
			}
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_0-l_1, D, d, best_features);
			
			qda.a = make_2D_matrix(d, d, dinit);
			qda.b = new double [d];
			qdaTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, prior, &qda);
			//temp_error0 += qdaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, qda);
			//temp_error1 += qdaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, qda);
			err0 = qdaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, qda);
			err0=(double)l_0*err0;
			temp_error00 += err0;
			err1 = qdaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, qda);
			err1=(double)l_1*err1;
			temp_error11 += err1;
			
				
			delete qda.b;
			delete_2D_matrix(d, d, qda.a);
			delete best_features;
		}
	}
	
	

	
	/* This computes the error similaly to what follows with temp_error00 and temp_error11. 
	if (prior==2.00)
	{
		//printf("We are here 1");
		error = (1.00-((double)N_1/(double)N))*temp_error0/K/R+((double)N_1/(double)N)*temp_error1/K/R;
	} else {
		error = (1.00-prior)*temp_error0/K/R+(prior)*temp_error1/K/R;
	}
	 
	
	 printf("errorSep1=%.2f\n",error);
	 */
	
	//The error is computed in this way to be consistent with Mixture sampling
	if (prior==2.00)
	{
		
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error00/l00+((double)N_1/(double)N)*temp_error11/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error00/l00+(prior)*temp_error11/l11;
	}
	//printf("errorSep2=%.2f\n",error);
	
		delete_2D_matrix(N-l_0-l_1, D, X_trn);
		delete_2D_matrix(l_0, D, X_tst0);
		delete_2D_matrix(l_1, D, X_tst1);
		delete y_trn;
		delete y_tst0;
		delete y_tst1;
		delete shuffled_0;
		delete shuffled_1;
	
	return error;
}



//This is the incorrect CV
double qdaCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r, i0, i1;
	int l_01;
	int l_0, l_1;
	int l00 = 0;
	int l11 = 0;
	double temp_error0 = 0.00;
	double err0 = 0.00;
	double temp_error1 = 0.00;
	double err1 = 0.00;
	double error = 0.00;
	
	double** X_trn;
	int* y_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_tst0;
	int* y_tst1;
	
	int* shuffled_01;
	
	int* best_features;
	int N_0, N_1;

	
	QDA qda;
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	//printf("N_0=%d\n",N_0);
	//printf("N_1=%d\n",N_1);
	
	
	/* size of sample to leave out */
	l_01 = int(floor(N/K));
	//l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_01, D, dinit);
	y_trn = new int [N-l_01];

	
	shuffled_01 = new int [N];
	//shuffled_1 = new int [N_1];
	
	for (i=0; i<N; i++)
		shuffled_01[i] = i;
	
	//for (i=0; i<N_1; i++)
	//shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N, seed, shuffled_01);
		//shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_01; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_01[i]][j];
				
				y_trn[i] = y[shuffled_01[i]];
			}
			
			l_0 = 0;
			l_1 = 0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
				if (y[shuffled_01[i]]==0)
					l_0++;
				else
					l_1++;
			}
			
			if (l_0>0) {
				X_tst0 = make_2D_matrix(l_0, D, dinit);
				y_tst0 = new int [l_0];
			}
			if (l_1>0) {
				X_tst1 = make_2D_matrix(l_1, D, dinit);	
				y_tst1 = new int [l_1];
			}
			
			i0=0;
			i1=0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
				
				if (y[shuffled_01[i]]==0) {
					for (j=0; j<D; j++)
						X_tst0[i0][j] = X[shuffled_01[i]][j];
					
					y_tst0[i0] = y[shuffled_01[i]];
					l00++;
					i0++;
				} else {
					for (j=0; j<D; j++)
						X_tst1[i1][j] = X[shuffled_01[i]][j];
					
					y_tst1[i1] = y[shuffled_01[i]];
					l11++;
					i1++;
				}
			}
			
			for (i=(k+1)*l_01; i<N; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_01][j] = X[shuffled_01[i]][j];
				
				y_trn[i-l_01] = y[shuffled_01[i]];
			}
			
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_01, D, d, best_features);
			
			qda.a = make_2D_matrix(d, d, dinit);
			qda.b = new double [d];
			qdaTrn(X_trn, y_trn, N-l_01, d, best_features, prior, &qda);
			if (l_0>0) {
				err0 = qdaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, qda);
				err0=(double)l_0*err0;
				temp_error0 += err0;
			}
			if (l_1>0) {
				err1 = qdaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, qda);
				//printf("err1=%.2f\n",err1);
				err1=(double)l_1*err1;
				//printf("err1.2=%.2f\n",err1);
				temp_error1 += err1;
		    }
			
			
			
			
			delete qda.b;
			delete_2D_matrix(d, d, qda.a);
			delete best_features;
			if (l_0>0) {
				delete_2D_matrix(l_0, D, X_tst0);
				delete y_tst0;
			}
			if (l_1>0) {
				delete_2D_matrix(l_1, D, X_tst1);
				delete y_tst1;
		    }
			
		}
	}
	
	//int KKK;
	//printf("l11=%d\n",l11);
	//KKK= N_1*R;
	//printf("N_1*R=%d\n",KKK);
	//double llll;
	//llll=temp_error0/l00;
	//printf("temp_error0/l00=%.2f\n",llll);
	
	if (prior==2.00)
	{
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error0/l00+((double)N_1/(double)N)*temp_error1/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error0/l00+(prior)*temp_error1/l11;
	}
	
	//printf("errorMix=%.2f\n",error);

	
	delete_2D_matrix(N-l_01, D, X_trn);
	delete y_trn;
	delete shuffled_01;
	//delete shuffled_1;
	
	return error;
}


/////////////////



double lsvmCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r;
	int l_0, l_1;
	int N_0, N_1;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double temp_error00 = 0.00;
	double temp_error11 = 0.00;
	double error = 0.00;
	
	double** X_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_trn;
	int* y_tst0;
	int* y_tst1;
	
	int l00 = 0;
	int l11 = 0;
	double err0 = 0.00;
	double err1 = 0.00;
	
	int* shuffled_0;
	int* shuffled_1;
	
	int* best_features;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_0 = int(floor(N_0/K));
	l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_0-l_1, D, dinit);
	y_trn = new int [N-l_0-l_1];
	X_tst0 = make_2D_matrix(l_0, D, dinit);
	X_tst1 = make_2D_matrix(l_1, D, dinit);	
	y_tst0 = new int [l_0];
	y_tst1 = new int [l_1];
	
	shuffled_0 = new int [N_0];
	shuffled_1 = new int [N_1];
	
	for (i=0; i<N_0; i++)
		shuffled_0[i] = i;
	
	for (i=0; i<N_1; i++)
		shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N_0, seed, shuffled_0);
		shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_0[i]][j];
				
				y_trn[i] = y[shuffled_0[i]];
			}
			
			for (i=k*l_0; i<(k+1)*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_tst0[i-k*l_0][j] = X[shuffled_0[i]][j];
				
				y_tst0[i-k*l_0] = y[shuffled_0[i]];
				l00++;
			}
			
			for (i=(k+1)*l_0; i<N_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_0][j] = X[shuffled_0[i]][j];
				
				y_trn[i-l_0] = y[shuffled_0[i]];
			}
			
			
			for (i=0; i<k*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i][j] = X[N_0+shuffled_1[i]][j];
				
				y_trn[N_0-l_0+i] = y[N_0+shuffled_1[i]];
			}
			
			for (i=k*l_1; i<(k+1)*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_tst1[i-k*l_1][j] = X[N_0+shuffled_1[i]][j];
				
				y_tst1[i-k*l_1] = y[N_0+shuffled_1[i]];
				l11++;
			}
			
			for (i=(k+1)*l_1; i<N_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i-l_1][j] = X[N_0+shuffled_1[i]][j];
				
				y_trn[N_0-l_0+i-l_1] = y[N_0+shuffled_1[i]];
			}
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_0-l_1, D, d, best_features);
			
			
			
			svm = svmTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, 0, &subdata_svm, &subcl_svm);
			//printf("Here\n");
			//temp_error0 += ldaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, lda);//pior is set to 2.00 so it simply devides by total samples
			//temp_error1 += ldaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, lda);
			err0 = svmTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, svm);
			err0=(double)l_0*err0;
			temp_error00 += err0;
			err1 = svmTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, svm);
			err1=(double)l_1*err1;
			temp_error11 += err1;

			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
		}
	}

	if (prior==2.00)
	{
		
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error00/l00+((double)N_1/(double)N)*temp_error11/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error00/l00+(prior)*temp_error11/l11;
	}
	


	delete_2D_matrix(N-l_0-l_1, D, X_trn);
	delete_2D_matrix(l_0, D, X_tst0);
	delete_2D_matrix(l_1, D, X_tst1);
	delete y_trn;
	delete y_tst0;
	delete y_tst1;
	delete shuffled_0;
	delete shuffled_1;

	return error;
}



double lsvmCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r, i0, i1;
	int l_01;
	int l_0, l_1;
	double temp_error0 = 0.00;
	double err0 = 0.00;
	double temp_error1 = 0.00;
	double err1 = 0.00;
	double error = 0.00;
	int l00 = 0;
	int l11 = 0;
	
	double** X_trn;
	int* y_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_tst0;
	int* y_tst1;
	
	int* shuffled_01;
	
	int* best_features;
	int N_0, N_1;
	
	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure
	
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_01 = int(floor(N/K));
	//l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_01, D, dinit);
	y_trn = new int [N-l_01];
	
	shuffled_01 = new int [N];
	//shuffled_1 = new int [N_1];
	
	for (i=0; i<N; i++)
		shuffled_01[i] = i;
	
	//for (i=0; i<N_1; i++)
	//shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N, seed, shuffled_01);
		//shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_01; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_01[i]][j];
				
				y_trn[i] = y[shuffled_01[i]];
			}
			
			l_0 = 0;
			l_1 = 0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
				if (y[shuffled_01[i]]==0)
					l_0++;
				else
					l_1++;
			}
			
			if (l_0>0) {
				X_tst0 = make_2D_matrix(l_0, D, dinit);
				y_tst0 = new int [l_0];
			}
			if (l_1>0) {
				X_tst1 = make_2D_matrix(l_1, D, dinit);	
				y_tst1 = new int [l_1];
			}
			
			i0=0;
			i1=0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
	            
				if (y[shuffled_01[i]]==0) {
					for (j=0; j<D; j++)
						X_tst0[i0][j] = X[shuffled_01[i]][j];
					
					y_tst0[i0] = y[shuffled_01[i]];
					l00++;
					i0++;
				} else {
					for (j=0; j<D; j++)
						X_tst1[i1][j] = X[shuffled_01[i]][j];
					
					y_tst1[i1] = y[shuffled_01[i]];
					l11++;
					i1++;
				}
			}
			
			for (i=(k+1)*l_01; i<N; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_01][j] = X[shuffled_01[i]][j];
				
				y_trn[i-l_01] = y[shuffled_01[i]];
			}
			
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_01, D, d, best_features);
			
			svm = svmTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, 0, &subdata_svm, &subcl_svm);

			if (l_0>0) {
				err0 = svmTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, svm);
				err0=(double)l_0*err0;
				temp_error0 += err0;
			}
			if (l_1>0) {
				err1 = svmTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, svm);
				err1=(double)l_1*err1;
				temp_error1 += err1;
			}
			
			
			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
			if (l_0>0) {
				delete_2D_matrix(l_0, D, X_tst0);
				delete y_tst0;
			}
			if (l_1>0) {
				delete_2D_matrix(l_1, D, X_tst1);
				delete y_tst1;
			}
			/*
			 printf("l_0=%d\n",l_0);
			 printf("l_1=%d\n",l_1);
			 printf("i0=%d\n",i0);
			 printf("i1=%d\n",i1);
			 printf("++++++\n",i1);
				 */
				
			}
		}

	
	if (prior==2.00)
	{
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error0/l00+((double)N_1/(double)N)*temp_error1/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error0/l00+(prior)*temp_error1/l11;
	}
	// printf("errorMix=%.2f\n",error);
	
	delete_2D_matrix(N-l_01, D, X_trn);
	delete y_trn;
	delete shuffled_01;
	//delete shuffled_1;
	
	return error;
}

////////////////////////////////////////////////////////////////////////////////////////////////////







double ksvmCVkFold(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r;
	int l_0, l_1;
	int N_0, N_1;
	double temp_error0 = 0.00;
	double temp_error1 = 0.00;
	double temp_error00 = 0.00;
	double temp_error11 = 0.00;
	double error = 0.00;
	
	double** X_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_trn;
	int* y_tst0;
	int* y_tst1;
	
	int l00 = 0;
	int l11 = 0;
	double err0 = 0.00;
	double err1 = 0.00;
	
	int* shuffled_0;
	int* shuffled_1;
	
	int* best_features;
    
	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure
    
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_0 = int(floor(N_0/K));
	l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_0-l_1, D, dinit);
	y_trn = new int [N-l_0-l_1];
	X_tst0 = make_2D_matrix(l_0, D, dinit);
	X_tst1 = make_2D_matrix(l_1, D, dinit);
	y_tst0 = new int [l_0];
	y_tst1 = new int [l_1];
	
	shuffled_0 = new int [N_0];
	shuffled_1 = new int [N_1];
	
	for (i=0; i<N_0; i++)
		shuffled_0[i] = i;
	
	for (i=0; i<N_1; i++)
		shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N_0, seed, shuffled_0);
		shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_0[i]][j];
				
				y_trn[i] = y[shuffled_0[i]];
			}
			
			for (i=k*l_0; i<(k+1)*l_0; i++)
			{
				for (j=0; j<D; j++)
					X_tst0[i-k*l_0][j] = X[shuffled_0[i]][j];
				
				y_tst0[i-k*l_0] = y[shuffled_0[i]];
				l00++;
			}
			
			for (i=(k+1)*l_0; i<N_0; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_0][j] = X[shuffled_0[i]][j];
				
				y_trn[i-l_0] = y[shuffled_0[i]];
			}
			
			
			for (i=0; i<k*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i][j] = X[N_0+shuffled_1[i]][j];
				
				y_trn[N_0-l_0+i] = y[N_0+shuffled_1[i]];
			}
			
			for (i=k*l_1; i<(k+1)*l_1; i++)
			{
				for (j=0; j<D; j++)
					X_tst1[i-k*l_1][j] = X[N_0+shuffled_1[i]][j];
				
				y_tst1[i-k*l_1] = y[N_0+shuffled_1[i]];
				l11++;
			}
			
			for (i=(k+1)*l_1; i<N_1; i++)
			{
				for (j=0; j<D; j++)
					X_trn[N_0-l_0+i-l_1][j] = X[N_0+shuffled_1[i]][j];
				
				y_trn[N_0-l_0+i-l_1] = y[N_0+shuffled_1[i]];
			}
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_0-l_1, D, d, best_features);
			
			
			
			svm = svmTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, 2, &subdata_svm, &subcl_svm);
			//printf("Here\n");
			//temp_error0 += ldaTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, lda);//pior is set to 2.00 so it simply devides by total samples
			//temp_error1 += ldaTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, lda);
			err0 = svmTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, svm);
			err0=(double)l_0*err0;
			temp_error00 += err0;
			err1 = svmTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, svm);
			err1=(double)l_1*err1;
			temp_error11 += err1;
            
			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
		}
	}
    
	if (prior==2.00)
	{
		
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error00/l00+((double)N_1/(double)N)*temp_error11/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error00/l00+(prior)*temp_error11/l11;
	}
	
    
    
	delete_2D_matrix(N-l_0-l_1, D, X_trn);
	delete_2D_matrix(l_0, D, X_tst0);
	delete_2D_matrix(l_1, D, X_tst1);
	delete y_trn;
	delete y_tst0;
	delete y_tst1;
	delete shuffled_0;
	delete shuffled_1;
    
	return error;
}



double ksvmCVkFold_Mixture(double** X, int* y, int N, int D, int d, int K, int R, double prior, long* seed)
{
	double dinit = 0.00;
	int i, j, k, r, i0, i1;
	int l_01;
	int l_0, l_1;
	double temp_error0 = 0.00;
	double err0 = 0.00;
	double temp_error1 = 0.00;
	double err1 = 0.00;
	double error = 0.00;
	int l00 = 0;
	int l11 = 0;
	
	double** X_trn;
	int* y_trn;
	double** X_tst0;
	double** X_tst1;
	int* y_tst0;
	int* y_tst1;
	
	int* shuffled_01;
	
	int* best_features;
	int N_0, N_1;
	
	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure
	
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	
	/* size of sample to leave out */
	l_01 = int(floor(N/K));
	//l_1 = int(floor(N_1/K));
	
	X_trn = make_2D_matrix(N-l_01, D, dinit);
	y_trn = new int [N-l_01];
	
	shuffled_01 = new int [N];
	//shuffled_1 = new int [N_1];
	
	for (i=0; i<N; i++)
		shuffled_01[i] = i;
	
	//for (i=0; i<N_1; i++)
	//shuffled_1[i] = i;
	
	for (r=0; r<R; r++)
	{
		shuffleArray(N, seed, shuffled_01);
		//shuffleArray(N_1, seed, shuffled_1);
		
		for (k=0; k<K; k++)
		{
			for (i=0; i<k*l_01; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i][j] = X[shuffled_01[i]][j];
				
				y_trn[i] = y[shuffled_01[i]];
			}
			
			l_0 = 0;
			l_1 = 0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
				if (y[shuffled_01[i]]==0)
					l_0++;
				else
					l_1++;
			}
			
			if (l_0>0) {
				X_tst0 = make_2D_matrix(l_0, D, dinit);
				y_tst0 = new int [l_0];
			}
			if (l_1>0) {
				X_tst1 = make_2D_matrix(l_1, D, dinit);
				y_tst1 = new int [l_1];
			}
			
			i0=0;
			i1=0;
			for (i=k*l_01; i<(k+1)*l_01; i++)
			{
	            
				if (y[shuffled_01[i]]==0) {
					for (j=0; j<D; j++)
						X_tst0[i0][j] = X[shuffled_01[i]][j];
					
					y_tst0[i0] = y[shuffled_01[i]];
					l00++;
					i0++;
				} else {
					for (j=0; j<D; j++)
						X_tst1[i1][j] = X[shuffled_01[i]][j];
					
					y_tst1[i1] = y[shuffled_01[i]];
					l11++;
					i1++;
				}
			}
			
			for (i=(k+1)*l_01; i<N; i++)
			{
				for (j=0; j<D; j++)
					X_trn[i-l_01][j] = X[shuffled_01[i]][j];
				
				y_trn[i-l_01] = y[shuffled_01[i]];
			}
			
			
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N-l_01, D, d, best_features);
			
			svm = svmTrn(X_trn, y_trn, N-l_0-l_1, d, best_features, 2, &subdata_svm, &subcl_svm);
            
			if (l_0>0) {
				err0 = svmTst(X_tst0, y_tst0, l_0, d, best_features,  2.00, svm);
				err0=(double)l_0*err0;
				temp_error0 += err0;
			}
			if (l_1>0) {
				err1 = svmTst(X_tst1, y_tst1, l_1, d, best_features,  2.00, svm);
				err1=(double)l_1*err1;
				temp_error1 += err1;
			}
			
			
			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
			if (l_0>0) {
				delete_2D_matrix(l_0, D, X_tst0);
				delete y_tst0;
			}
			if (l_1>0) {
				delete_2D_matrix(l_1, D, X_tst1);
				delete y_tst1;
			}
			/*
			 printf("l_0=%d\n",l_0);
			 printf("l_1=%d\n",l_1);
			 printf("i0=%d\n",i0);
			 printf("i1=%d\n",i1);
			 printf("++++++\n",i1);
             */
            
        }
    }
    
	
	if (prior==2.00)
	{
		//printf("We are here 1: Prir=%.2f",prior);
		error = ((double)N_0/(double)N)*temp_error0/l00+((double)N_1/(double)N)*temp_error1/l11;
	} else {
		//printf("We are here 2: Prir=%.2f",prior);
		error = (1.00-prior)*temp_error0/l00+(prior)*temp_error1/l11;
	}
	// printf("errorMix=%.2f\n",error);
	
	delete_2D_matrix(N-l_01, D, X_trn);
	delete y_trn;
	delete shuffled_01;
	//delete shuffled_1;
	
	return error;
}














double ldaBoot632(double** X, int* y, int N, int D, int d, int B, double resub, double prior, long* seed)
{
	int min_s = d-1; // min number of samples
	double dinit = 0.00;
	int b, i, j, k, N_tst, N_tst_error, ind, P_0, P_1, N_0, N_1;
	double temp_error = 0.00;
	double error = 0.00;

	int* P;
	double** X_trn;
	double** X_tst;
	int* y_trn;

	int* best_features;

	LDA lda;

	P = new int [N];
	X_trn = make_2D_matrix(N, D, dinit);
	y_trn = new int [N];
	X_tst = make_2D_matrix(1, D, dinit);

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	for (b=0; b<B; b++)
	{
		for (i=0; i<N; i++)
			P[i] = 0;

		for (i=0; i<N; i++)
		{
			ind = (int)(ran2(seed)*N); /* i+1 so that includes boundaries */
			for (j=0; j<D; j++)
				X_trn[i][j] = X[ind][j];
			y_trn[i] = y[ind];
			P[ind]=1;
		}

		P_0 = 0;
		P_1 = 0;
		for (i=0; i<N; i++)
		{
			if (y[i]==0 && P[i]!=0)
				P_0++;
			else if (y[i]==1 && P[i]!=0)
				P_1++;
		}

		if (P_0>min_s && P_1>min_s)
		{
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N, D, d, best_features);

			lda.a = new double [d];

			ldaTrn(X_trn, y_trn, N, d, best_features, prior, &lda);

			N_tst = 0;
			N_tst_error = 0;
			for (i=0; i<N; i++)
			{
				if (P[i]==0)
				{
					N_tst++;
					for (j=0; j<D; j++)
						X_tst[0][j] = X[i][j];

					N_tst_error += (int)ldaTst(X_tst, y+i, 1, d, best_features,  prior, lda);
				}
			}

			temp_error += (double) N_tst_error/N_tst;
			delete lda.a;
			delete best_features;
		}
		else
		{
			if (P_0<=min_s && P_1<=min_s)
				temp_error += P_1>P_0 ? (double)(N_0-P_0)/(N-P_0-P_1) : (double)(N_1-P_1)/(N-P_0-P_1);
			else if (P_0<=min_s)
				temp_error += (double) (N_0-P_0)/(N-P_0-P_1);
			else
				temp_error += (double) (N_1-P_1)/(N-P_0-P_1);
		}
	}

	error = 0.632*temp_error/B + 0.368*resub;

	delete_2D_matrix(N, D, X_trn);
	delete_2D_matrix(1, D, X_tst);
	delete y_trn;
	delete P;

	return error;
}

//this function is not yet working
double qdaBoot632(double** X, int* y, int N, int D, int d, int B, double resub, double prior, long* seed)
{
	int min_s = d-1; // min number of samples
	double dinit = 0.00;
	int b, i, j, k, N_tst, N_tst_error, ind, P_0, P_1, N_0, N_1;
	double temp_error = 0.00;
	double error = 0.00;
	
	int* P;
	double** X_trn;
	double** X_tst;
	int* y_trn;
	
	int* best_features;
	
	LDA lda;
	
	P = new int [N];
	X_trn = make_2D_matrix(N, D, dinit);
	y_trn = new int [N];
	X_tst = make_2D_matrix(1, D, dinit);
	
	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;
	
	for (b=0; b<B; b++)
	{
		for (i=0; i<N; i++)
			P[i] = 0;
		
		for (i=0; i<N; i++)
		{
			ind = (int)(ran2(seed)*N); /* i+1 so that includes boundaries */
			for (j=0; j<D; j++)
				X_trn[i][j] = X[ind][j];
			y_trn[i] = y[ind];
			P[ind]=1;
		}
		
		P_0 = 0;
		P_1 = 0;
		for (i=0; i<N; i++)
		{
			if (y[i]==0 && P[i]!=0)
				P_0++;
			else if (y[i]==1 && P[i]!=0)
				P_1++;
		}
		
		if (P_0>min_s && P_1>min_s)
		{
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N, D, d, best_features);
			
			lda.a = new double [d];
			
			ldaTrn(X_trn, y_trn, N, d, best_features, prior, &lda);
			
			N_tst = 0;
			N_tst_error = 0;
			for (i=0; i<N; i++)
			{
				if (P[i]==0)
				{
					N_tst++;
					for (j=0; j<D; j++)
						X_tst[0][j] = X[i][j];
					
					N_tst_error += (int)ldaTst(X_tst, y+i, 1, d, best_features,  prior, lda);
				}
			}
			
			temp_error += (double) N_tst_error/N_tst;
			delete lda.a;
			delete best_features;
		}
		else
		{
			if (P_0<=min_s && P_1<=min_s)
				temp_error += P_1>P_0 ? (double)(N_0-P_0)/(N-P_0-P_1) : (double)(N_1-P_1)/(N-P_0-P_1);
			else if (P_0<=min_s)
				temp_error += (double) (N_0-P_0)/(N-P_0-P_1);
			else
				temp_error += (double) (N_1-P_1)/(N-P_0-P_1);
		}
	}
	
	error = 0.632*temp_error/B + 0.368*resub;
	
	delete_2D_matrix(N, D, X_trn);
	delete_2D_matrix(1, D, X_tst);
	delete y_trn;
	delete P;
	
	return error;
}


double lsvmBoot632(double** X, int* y, int N, int D, int d, int B, double resub, long* seed)
{
	int min_s = 0; // min number of samples
	double dinit = 0.00;
	int b, i, j, k, N_tst, N_tst_error, ind, P_0, P_1, N_0, N_1;
	double temp_error = 0.00;
	double error = 0.00;

	int* P;
	double** X_trn;
	double** X_tst;
	int* y_trn;

	int* best_features;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	P = new int [N];
	X_trn = make_2D_matrix(N, D, dinit);
	y_trn = new int [N];
	X_tst = make_2D_matrix(1, D, dinit);

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	for (b=0; b<B; b++)
	{
		for (i=0; i<N; i++)
			P[i] = 0;

		for (i=0; i<N; i++)
		{
			ind = (int)(ran2(seed)*N); /* i+1 so that includes boundaries */
			for (j=0; j<D; j++)
				X_trn[i][j] = X[ind][j];
			y_trn[i] = y[ind];
			P[ind]=1;
		}

		P_0 = 0;
		P_1 = 0;
		for (i=0; i<N; i++)
		{
			if (y[i]==0 && P[i]!=0)
				P_0++;
			else if (y[i]==1 && P[i]!=0)
				P_1++;
		}

		if (P_0>min_s && P_1>min_s)
		{
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N, D, d, best_features);

			svm = svmTrn(X_trn, y_trn, N, d, best_features, 0, &subdata_svm, &subcl_svm);

			N_tst = 0;
			N_tst_error = 0;
			for (i=0; i<N; i++)
			{
				if (P[i]==0)
				{
					N_tst++;
					for (j=0; j<D; j++)
						X_tst[0][j] = X[i][j];

					N_tst_error += (int)svmTst(X_tst, y+i, 1, d, best_features,  2.00, svm);
				}
			}

			temp_error += (double) N_tst_error/N_tst;

			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
		}
		else
		{
			if (P_0<=min_s)
				temp_error += (double) N_0/(N-P_1);
			else
				temp_error += (double) N_1/(N-P_0);
		}
	}

	error = 0.632*temp_error/B + 0.368*resub;

	delete_2D_matrix(N, D, X_trn);
	delete_2D_matrix(1, D, X_tst);
	delete y_trn;
	delete P;

	return error;
}

double ksvmBoot632(double** X, int* y, int N, int D, int d, int B, double resub, long* seed)
{
	int min_s = 0; // min number of samples
	double dinit = 0.00;
	int b, i, j, k, N_tst, N_tst_error, ind, P_0, P_1, N_0, N_1;
	double temp_error = 0.00;
	double error = 0.00;

	int* P;
	double** X_trn;
	double** X_tst;
	int* y_trn;

	int* best_features;

	svm_model *svm;	//svm training model
	svm_node *subdata_svm;	//svm training data
	svm_problem subcl_svm;	//svm training data structure

	P = new int [N];
	X_trn = make_2D_matrix(N, D, dinit);
	y_trn = new int [N];
	X_tst = make_2D_matrix(1, D, dinit);

	N_0 = 0;
	N_1 = 0;
	for (i=0; i<N; i++)
		if (y[i]==0)
			N_0++;
		else
			N_1++;

	for (b=0; b<B; b++)
	{
		for (i=0; i<N; i++)
			P[i] = 0;

		for (i=0; i<N; i++)
		{
			ind = (int)(ran2(seed)*N); /* i+1 so that includes boundaries */
			for (j=0; j<D; j++)
				X_trn[i][j] = X[ind][j];
			y_trn[i] = y[ind];
			P[ind]=1;
		}

		P_0 = 0;
		P_1 = 0;
		for (i=0; i<N; i++)
		{
			if (y[i]==0 && P[i]!=0)
				P_0++;
			else if (y[i]==1 && P[i]!=0)
				P_1++;
		}

		if (P_0>min_s && P_1>min_s)
		{
			best_features = new int [d];
			featureSelection(X_trn, y_trn, N, D, d, best_features);

			svm = svmTrn(X_trn, y_trn, N, d, best_features, 2, &subdata_svm, &subcl_svm);

			N_tst = 0;
			N_tst_error = 0;
			for (i=0; i<N; i++)
			{
				if (P[i]==0)
				{
					N_tst++;
					for (j=0; j<D; j++)
						X_tst[0][j] = X[i][j];

					N_tst_error += (int)svmTst(X_tst, y+i, 1, d, best_features,  2.00, svm);
				}
			}

			temp_error += (double) N_tst_error/N_tst;

			svmDestroy(svm, subdata_svm, &subcl_svm);
			delete best_features;
		}
		else
		{
			if (P_0<=min_s)
				temp_error += (double) N_0/(N-P_1);
			else
				temp_error += (double) N_1/(N-P_0);
		}
	}

	error = 0.632*temp_error/B + 0.368*resub;

	delete_2D_matrix(N, D, X_trn);
	delete_2D_matrix(1, D, X_tst);
	delete y_trn;
	delete P;

	return error;
}
