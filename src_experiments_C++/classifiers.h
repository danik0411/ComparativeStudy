/*
 * classifiers.h
 *
 *  Created on: Feb 1, 2014
 *      Author: Yousef, Amin, Chao
 */

#ifndef CLASSIFIERS_H_
#define CLASSIFIERS_H_

#include "standardHeaders.h"
#include "matrixOperations.h"

struct regparam {
  double *gamma;
  double optGamma;
  int amount;
};

struct ZAR {
  double* a;
  double b;
};

struct SERD {
  double *a;
  double b;
};

struct DLDA {
  double *a;
  double b;
};

struct G13 {
  double *a;
  double b;
};

struct EDC {
  double *a;
  double b;
};

struct RLDA {
  double *a;
  double b;
};

struct LDA {
  double *a;
  double b;
};

struct QDA {
  double **a;
  double *b;
  double c;
};

int findIndOfMinElement( double* inputArray, int lastElement);
void estimateOptimumGamma(double* x0, double* x1, double** pooled_cov, regparam* gammas, int d, int n0, int n1 );

void rldaTrn(double** X, int* y, int N, int d, int* ind, double prior, regparam* gammas, RLDA* rlda);
double rldaTst(double** X, int* y, int N, int d, int* ind, double prior, RLDA rlda);

void zarTrn( double** X, int*, int, int, int*, double, ZAR* zar);
double zarTst(double** X, int* y, int N, int d, int* ind, double prior, ZAR lda);

void divide_a_to_k_groups ( double* x, int d, int k, int m, double* x_temp);
void divide_A_to_k_groups ( double** x, int N, int d, int k, int m, double** x_temp);

void serdTrn(double** X, int* y, int N, int d, int* ind, double prior, SERD* serd, int k, int m);
void serdCovTrn(double** X, int* y, int N, int d, int* ind, double prior, SERD* serd, int k, int m);
double serdTst(double** X, int* y, int N, int d, int* ind, double prior, SERD serd, int k, int m);

void dldaTrn(double** X, int* y, int N, int d, int* ind, double prior, DLDA* dlda);
double dldaTst(double** X, int* y, int N, int d, int* ind, double prior, DLDA dlda);

void g13Trn(double** X, int* y, int N, int d, int* ind, double prior, G13* g13);
double g13Tst(double** X, int* y, int N, int d, int* ind, double prior, G13 g13);

void edcTrn(double** X, int* y, int N, int d, int* ind, double prior, EDC* edc);
double edcTst(double** X, int* y, int N, int d, int* ind, double prior, EDC edc);

void ldaTrn(double** X, int* y, int N, int d, int* ind, double prior, LDA* lda);
double ldaTst(double** X, int* y, int N, int d, int* ind, double prior, LDA lda);

void qdaTrn(double** X, int* y, int N, int d, int* ind, double prior, QDA* qda);
double qdaTst(double** X, int* y, int N, int d, int* ind, double prior, QDA qda);
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

struct svm_model *svmTrn(double** X, int* y, int N, int d, int* ind, int kernel_type, struct svm_node** subdata, struct svm_problem* subcl);
double svmTst(double** X, int* y, int N, int d, int* ind, double prior, svm_model* model);
void svmDestroy(struct svm_model*, struct svm_node*, struct svm_problem*);

struct Edge {
  int src, dest;
  int row, col;
  double weight;
};

struct Graph {
  int V, E;
  Edge* edge;
};

struct subset {
  int parent;
  int rank;
};

int myComp( const void*, const void*);
Graph* createGraph(int V, int E);
int find(subset subsets[], int i);
void Union(subset subsets[], int x, int y);
//void KruskalMST(Graph* graph, double** totalCovMatrix, double** finalCovMatrix, int size);
//void maxSpanningTreeMatrix ( int size , double** totalCovMatrix , double** finalCovMatrix );
void KruskalMST(Graph* graph, double** totalCovMatrix, Edge *result, int size);
void maxSpanningTreeMatrix ( int size , double** totalCovMatrix , Edge *result );
void ZarFOTT( double** corMat, int size, double** ZarCovMat);

void svm_trueFalse_predict(const svm_model *model, const svm_node *x, double* est_value );

#ifdef __cplusplus
extern "C" {
#endif

	struct svm_node
	{
		int index;
		double value;
		double quality;
	};

	struct svm_problem
	{
		int l;
		double *y;
		struct svm_node **x;
	};

	enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
	enum { LINEAR, POLY, RBF, SIGMOID };	/* kernel_type */

	struct svm_parameter
	{
		int svm_type;
		int kernel_type;
		double degree;	/* for poly */
		double gamma;	/* for poly/rbf/sigmoid */
		double coef0;	/* for poly/sigmoid */

		/* these are for training only */
		double cache_size; /* in MB */
		double eps;	/* stopping criteria */
		double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
		int nr_weight;		/* for C_SVC */
		int *weight_label;	/* for C_SVC */
		double* weight;		/* for C_SVC */
		double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
		double p;	/* for EPSILON_SVR */
		int shrinking;	/* use the shrinking heuristics */
	};

	struct svm_model
	{
		svm_parameter param;	// parameter
		int nr_class;		// number of classes, = 2 in regression/one class svm
		int l;			// total #SV
		svm_node **SV;		// SVs (SV[l])
		double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[n-1][l])
		double *rho;		// constants in decision functions (rho[n*(n-1)/2])

		// for classification only

		int *label;		// label of each class (label[n])
		int *nSV;		// number of SVs for each class (nSV[n])
		// nSV[0] + nSV[1] + ... + nSV[n-1] = l
		// XXX
		int free_sv;		// 1 if svm_model is created by svm_load_model
		// 0 if svm_model is created by svm_train
	};

	struct svm_model *svm_train(const struct svm_problem *prob,
			const struct svm_parameter *param);

	int svm_save_model(const char*, const struct svm_model*);

	int get_sv(const struct svm_model*);

	struct svm_model *svm_load_model(const char*);

	double svm_predict(const struct svm_model*, const struct svm_node*);

	void svm_destroy_model(struct svm_model*);

	const char *svm_check_parameter(const struct svm_problem*, const struct svm_parameter*);

	void info(char*,...);

	void info_flush();

#ifdef __cplusplus
}
#endif


#endif /* CLASSIFIERS_H_ */
