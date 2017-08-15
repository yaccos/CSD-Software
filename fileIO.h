#include <string>
#include <vector>

#include "myStructs.h"

using namespace std;


void writeFiltered(vector<CSD_Pair>& csdNodes, vector< vector<Gene> >& filteredGenes, const int& noc, string outPrefix);


//void writeFiltered(vector<CSD_Pair>& csdNodes, const int& noc);


void readRawFilesNew(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector< vector< vector<float> > >& rawData, const int& threadCount);

/*
void readRawFiles(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector< vector< vector<double> > >& rawData);


void readRawFiles(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector<double**>& rawData);
// The function reads the contents from the list of files in the "inputFiles"-vector. The
// funciton uses the result from a preprocessing step to read the gene-data from all files
// so that they occur in the same gene sequence in their matrices.
*/


void writeAll(vector<vector<vector<float> > >& avgCorr, vector<vector<vector<float> > >& var, vector<vector<vector<vector<float> > > >& csd, const vector< vector<Gene> >& filteredGenes, const int& numOfGenes, string outPrefix);


void writeCSDtoFile(const vector< vector<Gene> >& filteredGenes, const vector<CSD_Pair>& csdNodes, const unsigned int csdSize);
// The function writes filtered CSD values to file

/*
void readFile(const string filename, const int& NUMOFMEAS, const int& NUMOFGENES, double **matrix, vector<string>& mIDs, vector<string>& gIDs);
// The function reads the content of a file on the format specified in the README-file.
// The gene-IDs are stored in the list gIDs, and the measurement-IDs in the list named mIDs.
// The measurement data are stored in the two dimensional input-matrix.
//
// Preconditions: - The data must be on the format specified in the README-file.
//                - The input matrix must be of appropriate size.
//                - There should be no empty space after the last column
//
// Input parameters:
//     - string: "filename", the name of the file containing the gene expression measurements
//     - const int: "NUMOFMEAS", the number of measurements/samples in the gene data set
//     - const int: "NUMOFGENES", the number of genes expression values in each sample
//     - double**: "matrix", a two dimensional matrix for storing the gene expression data
//     - list<string>&: "mIDs", a list for collecting the measurement-IDs in the order presented
//                      in the input file
//     - list<string>&: "gIDs", a list for collecting the gene-IDs in the order presented
//                      in the input file
*/


void readCorr(const string& inFile, vector<Gene>& filteredGenes, vector< vector<float> >& avgCorr, vector< vector<float> >& var, int& numOfGenes);


void writeCorr(const string& inFile, const vector<Gene>& filteredGenes, vector< vector<float> >& avgCorr, vector< vector<float> >& var, const int& NUMOFGENES);


void writeCorr(const string& inFile, const vector<Gene>& gIDs, double **avgCorr, double **var, const int& NUMOFGENES);
// The functions writes all the possible gene pairs, along with the mean correlation and
// the variance calculated from the generated subsamples.
//
// Precondition: - The input matrices only contains data above the diagonal. Hence iterating
//                 over the matrices only involves i=0->end , j=i->end.
//
// Input parameters:
//     - string: "inFile" is the filename of the related rawData-file.
//     - list<string>: "gIDs" is a list containing all the gene-IDs
//     - double**: "avgCorr", matrix containing the average correlation between all gene pairs
//     - double**: "var", matrix containing the variance in the correlation distribution of all
//                  the gene pairs.
//     - const int: "NUMOFGENES", the total number of genes



void readUnfiltered(const string& infile, vector<Gene>& filteredGenes, vector< vector< vector<float> > >& csd, int& numOfGenes);


void writeUnfiltered(const vector< vector< vector< vector<float> > > >& csd, const vector< vector<Gene> >& filteredGenes, const int& numOfGenes, const int& threadCount);



void readAllGeneIDs(const vector<string>& inputFiles, const bool& vertical, vector< vector<string> >& allIDs, const int& threadCount);
// The function reads all the gene-IDs from all the input files, and stores them in a two
// dimensional vector to used in preprocessing the input data.


int getNumOfCols(string filename);
// The function counts the number of columns in the input file, including the column with the IDs.
//
// Preconditions: - The data must be on the format specified in the README-file.
//
// Input parameters:
//     - string: "filename", name of the input file


int getNumOfRows(string filename);
// The function counts the number of rows in the input file, including the column with the IDs.
//
// Preconditions: - The data must be on the format specified in the README-file.
//
// Input parameters:
//     - string: "filename", name of the input file


void getColIDs(const string& filename, vector<string>& ids, const int& NUMOFCOLS);
// The function reads all the IDs at the top of each column in the input file, and
// stores them in a vector.


void getRowIDs(const string& filename, vector<string>& ids, const int& NUMOFROWS);
