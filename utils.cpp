#include <algorithm>

#include "utils.h"
#include "myStructs.h"

#include <iostream>  //remove after debug

using namespace std;

typedef double* DblArr;


void allocateMemory(double ***matrix, const int& dim1, const int& dim2, const int& dim3)
{
  matrix = new DblArr*[dim1];
  for (int i=0; i<dim2; ++i)
    {
      matrix[i] = new DblArr[dim2];
      for (int j=0; j<dim3; ++j)
	{
	  matrix[i][j] = new double[dim3];
	}
    }
}


void allocateMemory(double **matrix, const int& dim1, const int& dim2)
{
  matrix = new DblArr[dim1];
  for (int i=0; i<dim1; ++i) matrix[i] = new double[dim2];
}


void freeMemory(double ***matrix, const int& dim1, const int& dim2)
{
  for (int i=0; i<dim1; ++i)
    {
      for (int j=0; j<dim2; ++j)
	{
	  delete[] matrix[i][j];
	}
      delete[] matrix[i];
    }
  delete[] matrix;
}


void freeMemory(double **matrix, const int& dim1)
{
  for (int i=0; i<dim1; ++i) delete[] matrix[i];
  delete[] matrix;
}


bool geneLessThan(const Gene& a, const Gene& b)
{
  return a.geneID < b.geneID;
}


bool geneCompare(const Gene& a, const Gene& b)
{
  return a.geneID <= b.geneID;
}


bool indexCompare(const Gene& g1, const Gene& g2)
{
  return g1.index <= g2.index;
}


bool sortedPosCompare(const Gene& a, const Gene& b)
{
  return a.sortedPos <= b.sortedPos;
}


int numberOfCombinations(const int& numberOfFiles)
{
  int combos = 0;
  for (int i = 0; i<numberOfFiles; ++i) combos += i;
  return combos;
}



int invNumberOfCombinations(int m)
{
  int n=0;
  while (m>0)
    {
      ++n;
      m-=n;
    }

  if (m!=0)
    {
      cout << "Error: Invalid input to function invNumberOfCombinations()."
	   << "Input number must represent the number of combinations between an integer"
	   << "number of objects.\n";
      exit(1);
    }
  
  return n+1;
}
      



void quickSortGenes(vector<Gene>& A, int p, int r, bool (*compare)(const Gene&, const Gene&))
{
  if (p<r)
    {
      int q = partitionGenes(A, p, r, compare);
      quickSortGenes(A, p, q-1, compare);
      quickSortGenes(A, q+1, r, compare);
    }
}


int partitionGenes(vector<Gene>& A, int p, int r, bool (*compare)(const Gene&, const Gene&))
{
  int i = p-1;
  for (int j=p; j<r; ++j)
    {
      if (compare(A[j],A[r]))
	{
	  ++i;
	  //iter_swap(A.begin()+i, A.begin()+j);
	  geneSwap(A[i], A[j]);
	}
    }
  geneSwap(A[i+1], A[r]);
  //iter_swap(A[i+1], A[r]);
  return i+1;
}


void geneSwap(Gene& a, Gene& b)
{
  Gene temp;
  temp.geneID = a.geneID;
  temp.index = a.index;
  temp.sortedPos = a.sortedPos;
  temp.origin = a.origin;

  a.geneID = b.geneID;
  a.index = b.index;
  a.sortedPos = b.sortedPos;
  a.origin = b.origin;

  b.geneID = temp.geneID;
  b.index = temp.index;
  b.sortedPos = temp.sortedPos;
  b.origin = temp.origin;
}
