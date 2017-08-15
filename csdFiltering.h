#include <string>
#include <list>
#include <vector>

#include "myStructs.h"

using namespace std;



void dependencyScores(vector< vector<float> >& avgCorr1, vector< vector<float> >& avgCorr2, vector< vector<float> >& var1, vector< vector<float> >& var2, vector< vector< vector<float> > >& csd, const int& numOfGenes, const int& threadCount);


void dependencyScores(double **avgCorr1, double **avgCorr2, double **var1, double **var2, double ***csd, const int& numOfGenes);
// The function calculates the 'C', 'S' and 'D' scores given the correlation
// parameters for two data sets.


/*
void thresholdValues(double ***csd, const int& numOfGenes, const double& significance, double& cThresh, double& sThresh, double& dThresh);
// The function determines the treshold values for the CSD-scores given a
// certain level of significance.
*/
void thresholdValuesWithReplacement(vector< vector< vector<float> > >& csd, const int& numOfGenes, const double& significance, const int& numOfSamples, float& cThresh, float& sThresh, float& dThresh, int& rSeed);

void thresholdValues(vector< vector< vector<float> > >& csd, const int& numOfGenes, const double& significance, const int& numOfSamples, float& cThresh, float& sThresh, float& dThresh);


void thresholdValues(double ***csd, const int& numOfGenes, const double& significance, const int& numOfSamples, double& cThresh, double& sThresh, double& dThresh);


/*
void filterCSDscores(vector<double***>& csd, const int& NUMOFGENES, const double& cThresh, const double& sThresh, const double& dThresh, vector<CSD_Pair>& csdNodes);
// The function returns a vector containing all co-expressed genes.  
*/

void filterCSDscores(vector< vector< vector< vector<float> > > >& csd, const int& NUMOFGENES, const vector<float>& cThresh, const vector<float>& sThresh, const vector<float>& dThresh, vector<CSD_Pair>& csdNodes, const int& threadCount);


void filterCSDscores(vector<double***>& csd, const int& NUMOFGENES, const vector<double>& cThresh, const vector<double>& sThresh, const vector<double>& dThresh, vector<CSD_Pair>& csdNodes);

/*  
void filterCSDscores(double ***csd, const int& numOfGenes, const vector<double>& cTresh, const vector<double>& sTresh, const vector<double>& dTresh, vector<CSD_Pair>& csdNodes);
// The function returns a vector with all the significantly expressed gene pairs in
// the form of gene indices.
*/
