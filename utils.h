#include "myStructs.h"

using namespace std;


void allocateMemory(double ***matrix, const int& dim1, const int& dim2, const int& dim3);
// The function allocates memory for a 3D-matrix


void allocateMemory(double **matrix, const int& dim1, const int& dim2);
// The function allocates memory for a 2D-matrix


void freeMemory(double ***matrix, const int& dim1, const int& dim2);
// The function frees the memory used by a 3D-matrix


void freeMemory(double **matrix, const int& dim1);
// The function frees the memory used by a 2D-matrix


bool geneLessThan(const Gene& a, const Gene& b);


bool geneCompare(const Gene& a, const Gene& b);
// The function is used by the built in sort-function in the "algorithms"-library to compare
// data-structures of type "Gene" by their "geneID" member variable.


bool indexCompare(const Gene& g1, const Gene& g2);
// The function is used by the built in sort-function in the "algorithms"-library to compare
// data-structures of type "Gene" by their "index" member variable.


bool sortedPosCompare(const Gene& a, const Gene& b);


int numberOfCombinations(const int& numberOfFiles);
// The function returns the number of different ways to combine a group of n objects into
// pairs, where draw order is not taken into account.


int invNumberOfCombinations(int m);
// The function takes in the number of different ways of pairing an integer number of objects,
// and returns the number of objects.


void quickSortGenes(vector<Gene>& A, int p, int r, bool (*compare)(const Gene&, const Gene&));


int partitionGenes(vector<Gene>& A, int p, int r, bool (*compare)(const Gene&, const Gene&));


void geneSwap(Gene& a, Gene& b);
