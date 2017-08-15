#include <string>
#include <vector>

#include "myStructs.h"

using namespace std;




// FUNCTIONS EXECUTING PROGRAM
void runProgramRawInput(const vector<string>& inputFiles, const Output& o, const int& k, const double& s, const string& c, const bool& vertical, const int& iSamples, const int& threadCount, const string& pt_filename);


void runProgramRawInput(const vector<string>& inputFiles, const Output& o, const int& k, const double& s, const string& c, const bool& format, const int& iSamples, const int& threadCount, const string& pt_filename, int& rSeed, int& maxSubs, int& terminationLimit, const string& outPrefix);



void runProgramCorrInput(const vector<string>& inputFiles, const Output& o, const double& s, const int& iSamples, const int& threadCount, const string& pt_filename, const string& outPrefix);


void runProgramUnfilteredInput(const vector<string>& inputFiles, const double& s, const int& iSamples, const int& threadCount, const string& pt_filename, const string& outPrefix);


void rawToCorrMem_Eval(vector< vector<double> >& rawData, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& NUMOFMEAS, const int& NUMOFGENES, const string& corrMethod, const int& subsampleSize, const int& threadCount);


void rawToCorr(vector< vector<double> >& rawData, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& NUMOFMEAS, const int& NUMOFGENES, const string& corrMethod, const int& subsampleSize, const int& threadCount);







void runSpearmanCorrelation(vector<vector<float> >& rawData, vector<vector<float> >& avgCorr, const int& numOfMeas, const int& NUMOFGENES, const int& threadCount);


void calculateVariance(vector <vector<float> >& rawData, vector <vector<float> >& avgCorr, vector< vector<float> >& var, const int& numOfMeas, const int& NUMOFGENES, const int subsampleSize, const int& terminationLimit, const int& maxNumSubsamples, const int& randomSeed, const int& threadCount);
