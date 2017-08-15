#include <string>
#include <list>
#include <vector>

using namespace std;


void printFromFile(string filename);
// The function reads a textfile line by line
// and writes it to the console.


void printFromData(double **dataMatrix, int numOfMeas, int numOfGenes, vector<string>& measIDs, vector<string>& geneIDs);
// The function writes the gene data read using the
// functions described in readFile.h out to the
// console.


void printStringVector(vector<string>& strList);
// The function prints all the strings contained in the input
// list to the console.
