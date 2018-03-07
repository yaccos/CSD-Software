#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <sys/time.h>

#include "queryValidation.h"
#include "runProgram.h"
#include "myStructs.h"

using namespace std;


struct timeval start, endi;


int main(int argc, char* argv[])
{
  // For determining the wall-time of the program:
  gettimeofday(&start, NULL);

  
  // Number of threads for parallellization
  int threadCount = 1;

  
  // Input data type
  InputType input;

  
  // Input parameters
  int k = 7;   //default
  double s = 0.001;   //default
  string c = "spearman";   //spearman-correlation as default
  bool vertical = 1; //genes arranged vertically as default
  int iSamples = 10000; //the number of samples for threshold-estimation
  
  // Output data requested (may be multiple types)
  Output o;

  // Parameters for subsampling
  int maxSubs = 100; //default
  int terminationLimit = 1000; //default
  
  // Seed for random number generators
  int rSeed = 0; //rSeed=0 will cause the program to call srand(time(NULL))

  // Performance test filename
  string pt_filename = "pt_CSD_network.txt";

  // Name prefix of out files
  string outPrefix = "CSD_network";
  
  // Reading input parameters 
  vector<string> inputFiles;
  for (int i=1; i<argc; ++i)
    {
      if (!string(argv[i]).compare("-i")) inputFiles.push_back(string(argv[i+1]));
      if (!string(argv[i]).compare("-k")) k = atoi(argv[i+1]);
      if (!string(argv[i]).compare("-p")) s = atof(argv[i+1]);
      if (!string(argv[i]).compare("--threads")) threadCount = atoi(argv[i+1]);
      if (!string(argv[i]).compare("--input")) setInputType(string(argv[i+1]), input);
      if (!string(argv[i]).compare("--output")) setOutput(string(argv[i+1]), o);
      if (!string(argv[i]).compare("--corrMethod")) setCorrMethod(string(argv[i+1]), c);
      if (!string(argv[i]).compare("--format")) setFormat(string(argv[i+1]), vertical);
      if (!string(argv[i]).compare("--iSamples")) iSamples = atoi(argv[i+1]);
      if (!string(argv[i]).compare("--performanceTest")) pt_filename = string(argv[i+1]);
      if (!string(argv[i]).compare("--randomSeed")) rSeed = atoi(argv[i+1]);
      if (!string(argv[i]).compare("--maxSubsamples")) maxSubs = atoi(argv[i+1]);
      if (!string(argv[i]).compare("--terminationLimit")) terminationLimit = atoi(argv[i+1]);
      if (!string(argv[i]).compare("--outfilePrefix")) outPrefix = string(argv[i+1]);
    }




  
  // Checks for illegal combinations of input- and output-types.
  validateJob(input, o);

  
  // Call function for running the program
  if (input.raw) runProgramRawInput(inputFiles, o, k, s, c, vertical, iSamples, threadCount, pt_filename, rSeed, maxSubs, terminationLimit, outPrefix);
  else if (input.corr) runProgramCorrInput(inputFiles, o, s, iSamples, threadCount, pt_filename, outPrefix);
  else if (input.unfiltered) runProgramUnfilteredInput(inputFiles, s, iSamples, threadCount, pt_filename, outPrefix);
  else
    {
      cout << "Error: No input data type set.\n";
      exit(1);
    }



  // Determining wall-time
  gettimeofday(&endi, NULL);
  double delta = ((endi.tv_sec - start.tv_sec)* 1000000u + endi.tv_usec - start.tv_usec)/1.e6;

  cout << "\n\nWall-time:\t" << delta << "\tseconds\n";
  
  return 0;	   
}

