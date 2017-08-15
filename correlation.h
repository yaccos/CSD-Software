#include <list>
#include <vector>

using namespace std;

typedef double* DblPtr;

typedef float* FltPtr;


void meanAndVar(double ***corrMatrix, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& numOfGenes, const int& numOfSubsamples, const int& threadCount);


void meanAndVar(double ***corrMatrix, double **avgCorr, double **var, const int& numOfGenes, const int& numOfSubsamples); 
// The estimates the mean and variance of the correlation parameter from the
// generated subsamples.


void pearsonCorrelation(double **rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& subsampleSize);
// The function calculates the pearson correlation coefficient for all the
// input subsamples.


void pearsonCorrelation(vector< vector<float> >& rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& subsampleSize, const int& threadCount);


void spearmanCorrelation(double **rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& numOfMeas, const int& subsampleSize);
// The function calculates the spearman correlation coefficient for all the
// input subsamples.


void spearmanCorrelation(vector< vector<float> >& rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& numOfMeas, const int& subsampleSize, const int& threadCount);

  
void pearson(double *x, double *y, const int& n, double& r);
// The function calculates the pearson correlation between
// the two input arrays x and y. The correlation coefficient
// is returned in the parameter "r".
//
// Precondition: Both arrays are of same type and have length n.


void spearman(double *x, double *y, const int& n, double& rho);
// The function calculates the spearman correlation between
// the two input arrays x and y. The correlation coefficient
// is returned in the parameter "rho".
//
// Precondition: Both arrays are of same type and have length n.


void biWeightMidcorr(double *x, double *y, int n, double& bc);
// This function calculates the bi-weight midcorrelation
// between the two input arrays, which is returned in the
// "bc" parameter.
//
// Precondition: Both arrays are of same type and have length n.




// UTILITY FUNCTIONS

void pQuicksort(DblPtr *A, int p, int r);
// The function uses quicksort to sort an array of pointers
// based on the value they point to.


int pPartition(DblPtr *A, int p, int r);
// This is the partition function used by the pQuicksort
// function.


void qsRank(double *A, const int& n);
// QuickSortRank: The function replaces the values within
// the input array by their relative rank. The function uses
// the QuickSort algorithm.
//
// Postcondition: The array consists of integer ranks, and
//                should be considered typecasted to an
//                int-array.


void msRank(double *A, const int& n);
// MergeSortRank: The function replaces the values within
// the input array by their relative rank. The function uses
// the MergeSort algorithm.
//
// Postcondition: The array consists of integer ranks, and
//                should be considered typecasted to an
//                int array.


void msRank(float *A, const int& n);

void msRank(vector<double> A, const int& n);

void msRank(vector<float> A, const int& n);


void pSwap(DblPtr& a, DblPtr& b);
// Swap the memory location of two pointers.

void pSwap(FltPtr& a, FltPtr& b);


void pMerge(DblPtr *A, int p, int q, int r);
// This is the merge function used be the pMergeSort function.


void pMerge(FltPtr *A, int p, int q, int r);


void pMergeSort(DblPtr *A, int p, int r);
// The function sorts an array of pointers based on the double
// value pointed to by each of the pointers. The function uses
// the MergeSort algorithm to do so.


void pMergeSort(FltPtr *A, int p, int r);
