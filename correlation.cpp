#include <cmath>
#include <cfloat>
#include <list>
#include <vector>
#include <omp.h>

#include "correlation.h"

using namespace std;


typedef double* DblPtr;

typedef float* FltPtr;


void meanAndVar(double ***corrMatrix, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& numOfGenes, const int& numOfSubsamples, const int& threadCount)
{
# pragma omp parallel for num_threads(threadCount) schedule(static,1) 
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i+1; j<numOfGenes; ++j)
	{
	  avgCorr[i][j] = 0;
	  for (int k=0; k<numOfSubsamples; ++k)
	    {
	      avgCorr[i][j] += corrMatrix[i][j][k];
	    }
	  avgCorr[i][j] = avgCorr[i][j]/double(numOfSubsamples);
	}
    }
# pragma omp parallel for num_threads(threadCount) schedule(static,1)
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i+1; j<numOfGenes; ++j)
	{
	  var[i][j] = 0;
	  for (int k=0; k<numOfSubsamples; ++k)
	    {
	      var[i][j] += pow((corrMatrix[i][j][k] - avgCorr[i][j]),2);
	    }
	  var[i][j] = var[i][j]/double(numOfSubsamples-1);
	}
    }  
}



void meanAndVar(double ***corrMatrix, double **avgCorr, double **var, const int& numOfGenes, const int& numOfSubsamples)
{
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i; j<numOfGenes; ++j)
	{
	  avgCorr[i][j] = 0;
	  for (int k=0; k<numOfSubsamples; ++k)
	    {
	      avgCorr[i][j] += corrMatrix[i][j][k];
	    }
	  avgCorr[i][j] = avgCorr[i][j]/double(numOfSubsamples);
	}
    }
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i; j<numOfGenes; ++j)
	{
	  var[i][j] = 0;
	  for (int k=0; k<numOfSubsamples; ++k)
	    {
	      var[i][j] += pow((corrMatrix[i][j][k] - avgCorr[i][j]),2);
	    }
	  var[i][j] = avgCorr[i][j]/double(numOfSubsamples-1);
	}
    }  
}




void pearsonCorrelation(double **rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& subsampleSize)
{
  int l = 0;
  double x[subsampleSize];
  double y[subsampleSize];
  for (list< vector<int> >::const_iterator it = subsamples.begin(); it != subsamples.end(); ++it)
    {
      for (int i=0; i<numOfGenes; ++i)
	{
	  for (int j=i+1; j<numOfGenes; ++j)
	    {
	      for (int k=0; k<subsampleSize; ++k)
		{
		  x[k] = rawData[(*it)[k]][i];
		  y[k] = rawData[(*it)[k]][j];
		}
	      pearson(x,y,subsampleSize,corrMatrix[i][j][l]);
	    }
	}
      ++l;
    }
}


void pearsonCorrelation(vector< vector<float> >& rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& subsampleSize, const int& threadCount)
{
  int l = 0;
  double x[subsampleSize];
  double y[subsampleSize];
  for (list< vector<int> >::const_iterator it = subsamples.begin(); it != subsamples.end(); ++it)
    {
#     pragma omp parallel for num_threads(threadCount) private(x,y) schedule(static,1)
      for (int i=0; i<numOfGenes; ++i)
	{
	  for (int j=i+1; j<numOfGenes; ++j)
	    {
	      for (int k=0; k<subsampleSize; ++k)
		{
		  x[k] = rawData[(*it)[k]][i];
		  y[k] = rawData[(*it)[k]][j];
		}
	      pearson(x,y,subsampleSize,corrMatrix[i][j][l]);
	    }
	}
      ++l;
    }
}


void spearmanCorrelation(double **rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& numOfMeas, const int& subsampleSize)
{
  for (int i=0; i<numOfMeas; ++i)
    {
      msRank(rawData[i], numOfGenes);
    }
   
  int l = 0;
  double x[subsampleSize];
  double y[subsampleSize];
  for (list< vector<int> >::const_iterator it = subsamples.begin(); it != subsamples.end(); ++it)
    {
      for (int i=0; i<numOfGenes; ++i)
	{
	  for (int j=i+1; j<numOfGenes; ++j)
	    {
	      for (int k=0; k<subsampleSize; ++k)
		{
		  x[k] = rawData[(*it)[k]][i];
		  y[k] = rawData[(*it)[k]][j];
		}
	      pearson(x,y,subsampleSize,corrMatrix[i][j][l]);
	    }
	}
      ++l;
    }
}


void spearmanCorrelation(vector< vector<float> >& rawData, double ***corrMatrix, const list< vector<int> >& subsamples, const int& numOfGenes, const int& numOfMeas, const int& subsampleSize, const int& threadCount)
{
# pragma omp parallel for num_threads(threadCount) 
  for (int i=0; i<numOfMeas; ++i)
    {
      msRank(rawData[i], numOfGenes);
    }
   
  int l = 0;
  double x[subsampleSize];
  double y[subsampleSize];
  for (list< vector<int> >::const_iterator it = subsamples.begin(); it != subsamples.end(); ++it)
    {
#     pragma omp parallel for num_threads(threadCount) private(x,y) schedule(static,1)      
      for (int i=0; i<numOfGenes; ++i)
	{
	  for (int j=i+1; j<numOfGenes; ++j)
	    {
	      for (int k=0; k<subsampleSize; ++k)
		{
		  x[k] = rawData[(*it)[k]][i];
		  y[k] = rawData[(*it)[k]][j];
		}
	      pearson(x,y,subsampleSize,corrMatrix[i][j][l]);
	    }
	}
      ++l;
    }
}




void pearson(double *x, double *y, const int& n, double& r)
{
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;

  for (int i=0; i<n; ++i)
  {
    sx += x[i];
    sy += y[i];
    sxx += (x[i]*x[i]);
    syy += (y[i]*y[i]);
    sxy += (x[i]*y[i]);
  }

  sxx -= ((sx*sx)/n);
  syy -= ((sy*sy)/n);
  sxy -= ((sx*sy)/n);

  r = sxy/sqrt(sxx*syy);
}


void spearman(double *x, double *y, const int& n, double& rho)
{
  msRank(x, n);
  msRank(y, n);
  pearson(x, y, n, rho);
}

/*
void biWeightMidcorr(double *x, double *y, int n, double& bc)
{
  #TODO
}
*/

void pQuicksort(DblPtr *A, int p, int r)
{
  if (p<r)
  {
    int q = pPartition(A, p, r);
    pQuicksort(A, p, q-1);
    pQuicksort(A, q+1, r);
  }
}


int pPartition(DblPtr *A, int p, int r)
{
  int x = *A[r];
  int i = p-1;
  for (int j=p; j<r; ++j)
  {
    if (*A[j] <= x)
    {
      i = i+1;
      pSwap(A[i], A[j]);
    }
  }
  pSwap(A[i+1], A[r]);
  return i+1;
}


void qsRank(double *A, const int& n)
{
  DblPtr B[n];
  for (int i=0; i<n; ++i)
  {
    B[i] = &A[i];
  }

  pQuicksort(B, 0, n-1);

  for (int i=0; i<n; ++i)
  {
    *B[i] = i;
  }
}


void msRank(double *A, const int& n)
{
  DblPtr B[n];
  for (int i=0; i<n; ++i)
  {
    B[i] = &A[i];
  }

  pMergeSort(B, 0, n-1);
  
  for (int i=0; i<n; ++i)
  {
    *B[i] = i;
  }
}


void msRank(float *A, const int& n)
{
  FltPtr B[n];
  for (int i=0; i<n; ++i)
  {
    B[i] = &A[i];
  }

  pMergeSort(B, 0, n-1);
  
  for (int i=0; i<n; ++i)
  {
    *B[i] = i;
  }
}


void msRank(vector<double> A, const int& n)
{
  DblPtr B[n];
  for (int i=0; i<n; ++i)
  {
    B[i] = &A[i];
  }

  pMergeSort(B, 0, n-1);
  
  for (int i=0; i<n; ++i)
  {
    *B[i] = i;
  }
}


void msRank(vector<float> A, const int& n)
{
  FltPtr B[n];
  for (int i=0; i<n; ++i)
  {
    B[i] = &A[i];
  }

  pMergeSort(B, 0, n-1);
  
  for (int i=0; i<n; ++i)
  {
    *B[i] = i;
  }
}




void pSwap(DblPtr& a, DblPtr& b)
{
  double *temp = a;
  a = b;
  b = temp;
}

void pSwap(FltPtr& a, FltPtr& b)
{
  float *temp = a;
  a = b;
  b = temp;
}



void pMerge(DblPtr *A, int p, int q, int r)
{
  int n1 = q-p+1;
  int n2 = r-q;
  DblPtr L[n1+1], R[n2+1];

  for (int i=0; i<n1; ++i)
  {
    L[i] = A[p+i];
  }
  for (int j=0; j<n2; ++j)
  {
    R[j] = A[q+j+1];
  }
  double inf = DBL_MAX;
  L[n1] = &inf;
  R[n2] = &inf;
  int i = 0;
  int j = 0;
  for (int k = p; k<=r; ++k)
  {
    if (*L[i] <= *R[j])
    {
      A[k] = L[i];
      ++i;
    }
    else
    {
      A[k] = R[j];
      ++j;
    }
  }
}


void pMerge(FltPtr *A, int p, int q, int r)
{
  int n1 = q-p+1;
  int n2 = r-q;
  FltPtr L[n1+1], R[n2+1];

  for (int i=0; i<n1; ++i)
  {
    L[i] = A[p+i];
  }
  for (int j=0; j<n2; ++j)
  {
    R[j] = A[q+j+1];
  }
  float inf = FLT_MAX;
  L[n1] = &inf;
  R[n2] = &inf;
  int i = 0;
  int j = 0;
  for (int k = p; k<=r; ++k)
  {
    if (*L[i] <= *R[j])
    {
      A[k] = L[i];
      ++i;
    }
    else
    {
      A[k] = R[j];
      ++j;
    }
  }
}




void pMergeSort(DblPtr *A, int p, int r)
{
  if (p<r)
  {
    int q = (p+r)/2;
    pMergeSort(A, p, q);
    pMergeSort(A, q+1, r);
    pMerge(A, p, q, r);
  }
}


void pMergeSort(FltPtr *A, int p, int r)
{
  if (p<r)
  {
    int q = (p+r)/2;
    pMergeSort(A, p, q);
    pMergeSort(A, q+1, r);
    pMerge(A, p, q, r);
  }
}

