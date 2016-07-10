#include "featureSelection.h"
#include "print.h"


int featureSelection(double** X, int* y, int N, int D, int d, int* best_features)
{

	int* feature_rank;
	double* feature_temp_class0;
	double* feature_temp_class1;
	double* feature_quality;

	int i, j, k, l;
	int N_class0=0, N_class1=0, i_class0=0, i_class1=0;

	for (i=0; i<N; i++)
	  if (y[i]==0) {
			N_class0++;
	  }
	  else if(y[i] == 1) {
			N_class1++;
	  }

	feature_rank = new int [D];
	feature_temp_class0 = new double [N_class0];
	feature_temp_class1 = new double [N_class1];
	feature_quality = new double [D];

	for (i=0; i<D; i++)
		feature_rank[i] = i;

	for (i=0; i<D; i++)
	{
		for (j=0; j<N; j++)
		{
			if (y[j]==0)
			{
				feature_temp_class0[i_class0] = X[j][i];
				i_class0++;
			}
			else
			{
				feature_temp_class1[i_class1] = X[j][i];
				i_class1++;
			}
		}
		i_class0 = 0;
		i_class1 = 0;

		//evaluate quality of feature
		feature_quality[i]=(double)fabs((double)t_test(feature_temp_class0, feature_temp_class1, N_class0, N_class1));
		
	}

	bubble_sort(feature_quality, feature_rank, D, d);

	for (i=0; i<d; i++)
		best_features[i] = feature_rank[i];


	delete [] feature_rank;
	delete [] feature_temp_class0;
	delete [] feature_temp_class1;
	delete [] feature_quality;

	return 0;
}

//return the t-test results
double t_test(double* feature_class0, double* feature_class1, int n_class0, int n_class1)
{
	int k;
	double mean_tr[2], var_tr[2];
	double value;

	//compute the mean and standard deviation of selected training data
	mean_tr[0]=0;
	var_tr[0]=0;
	mean_tr[1]=0;
	var_tr[1]=0;

	for (k=0; k<n_class0; k++)
		mean_tr[0]+=feature_class0[k];

	for (k=0; k<n_class1; k++)
		mean_tr[1]+=feature_class1[k];

	mean_tr[0]/=n_class0;
	mean_tr[1]/=n_class1;

	for (k=0; k<n_class0; k++)
		var_tr[0]+=(feature_class0[k]-mean_tr[0])*(feature_class0[k]-mean_tr[0]);

	for (k=0; k<n_class1; k++)
		var_tr[1]+=(feature_class1[k]-mean_tr[1])*(feature_class1[k]-mean_tr[1]);

	if (n_class0<2)
	{
		var_tr[1]/=(n_class1-1);
		var_tr[0]=var_tr[1];
	}
	else if (n_class1<2)
	{
		var_tr[0]/=(n_class0-1);
		var_tr[1]=var_tr[0];
	}
	else
	{
		var_tr[1]/=(n_class1-1);
		var_tr[0]/=(n_class0-1);
	}

	if ((var_tr[0]+var_tr[1])==0 && mean_tr[0]==mean_tr[1])
		value=0;
	else if ((var_tr[0]+var_tr[1])==0)
	{
		if (mean_tr[0]>mean_tr[1])
			value=65535;
		else
			value=-65535;
	}
	else
		value=(mean_tr[0]-mean_tr[1])/sqrt((var_tr[0]*(double)(n_class0-1)+var_tr[1]*(double)(n_class1-1))/(double)(n_class0+n_class1-2)*(1/(double)n_class0+1/(double)n_class1));


	//printf("n0 is %d\tn1 is %d\nmean0 is %.3f\tmean1 is %.3f\nvar0 is %.3f\tvar1 is %.3f\n", n_class0, n_class1, mean_tr[0], mean_tr[1], var_tr[0], var_tr[1]);

	return(value);
}
