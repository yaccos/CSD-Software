#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <omp.h>

#include "queryValidation.h"
#include "preprocessing.h"
#include "fileIO.h"
#include "subsampling2.h"
#include "correlation.h"
#include "csdFiltering.h"
#include "runProgram.h"
#include "myStructs.h"
#include "utils.h"


using namespace std;

typedef float* FltArr;

struct timeval start_wt, end_wt;

struct timeval start_wt_corr, end_wt_corr;

struct timeval start_wt_var, end_wt_var;



void runProgramRawInput(const vector<string>& inputFiles, const Output& o, const int& k, const double& s, const string& c, const bool& format, const int& iSamples, const int& threadCount, const string& pt_filename, int& rSeed, int& maxSubs, int& terminationLimit, const string& outPrefix)
{
  // Setting for performance test
  ofstream fout(pt_filename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed for performance test outfile.\n";
      exit(1);
    }
  
  fout << "----------------------------------------------\n"
       << "              PERFORMANCE TEST\n"
       << "----------------------------------------------\n\n\n"
       << "Running times for program using gene expression data as input.\n\n"
       << "FUNCTIONS/SECTION\t\tTIME (s)\tTOTAL TIME (s)\n"
       << "-------------------------------------------------------------------\n";

  clock_t begin = clock();
  


  
  cout << "\nStart of runProgram...\n";
  cout << "\tNumber of input files: " << inputFiles.size() << endl;
  if (format) cout << "\tInput file format: vertical\n";
  else cout << "\tInput file format: horizontal\n";
  
  if (o.corr) cout << "\nOutput: corr";
  if (o.unfiltered) cout << "\nOutput: unfiltered";
  if (o.filtered) cout << "\nOutput: filtered";

  cout << endl;


  
  // PREPROCESSING INPUT DATA
  vector< vector<string> > allIDs;
  // Parallellized reading inside function, one thread for each file
  readAllGeneIDs(inputFiles, format, allIDs, threadCount);
  vector< vector<Gene> > allGenes;
  vector< vector<Gene> > filteredGenes;


  clock_t end = clock();
  double section_time = double(end-begin)/CLOCKS_PER_SEC;
  double total_time = section_time;
  begin = end;
  fout << "readAllGeneIDs()\t\t" << section_time << "\t" << total_time << endl;

  
  cout << "Number of genes in file 1: " << (allIDs[0]).size() << endl;
  cout << "Number of genes in file 2: " << (allIDs[1]).size() << endl;

  

  // Parallellized soring and initialization
  filterAndSortGenes(allIDs, allGenes, filteredGenes, threadCount);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "filterAndSortGenes()\t\t" << section_time << "\t" << total_time << endl;
  
  
  //
  //  SECTION ABOVE DEBUGGED AND TESTED: PERFORMS AS EXPECTED
  //
  
  
  cout << "Preprocessing done...\n";
  
  // READING FROM FILE(-S)
  const int NUMOFGENES = (filteredGenes[0]).size(); 
  vector<int> numOfMeas;
  numOfMeas.resize(int(inputFiles.size()));

# pragma omp parallel for num_threads(threadCount)
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      if (format) numOfMeas[i] = getNumOfCols(inputFiles[i])-1;
      else numOfMeas[i] = getNumOfRows(inputFiles[i])-1;
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "determining numOfMeas[i]\t" << section_time << "\t" << total_time << endl;
  
  
  cout << "\tReading numOfMeas done...\n";

  vector< vector< vector<float> > > rawData(inputFiles.size());
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      (rawData[i]).resize(numOfMeas[i]);
      for (int j=0; j<numOfMeas[i]; ++j)
	{
	  (rawData[i][j]).resize(NUMOFGENES);
	}
    }
  
  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for rawData\t" << section_time << "\t" << total_time << endl;

  
  cout << "\tAllocating memory for rawData done...\n";

  // Parallellized outer loop, one thread per input file
  readRawFilesNew(inputFiles, allGenes, filteredGenes, numOfMeas, format, rawData, threadCount);  

  /*
  for (unsigned int fi=0; fi<rawData.size(); ++fi)
    {
      for (int i=0; i<numOfMeas[fi]; ++i)
	{
	  for (int j=0; j<NUMOFGENES; ++j)
	    {
	      cout << "\t" << rawData[fi][i][j];
	    }
	  cout << endl;
	}
    }
  */
  
  
  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "readRawFilesNew\t" << section_time << "\t" << total_time << endl;

  
  cout << "File reading done...\n";

  /*
  cout << "Sample from rawData to verify file reading:\n"
       << rawData[0][5][4] << endl
       << rawData[0][3][8] << endl
       << rawData[1][12][43] << endl
       << rawData[1][20][17] << endl;
  */
  
  //
  //
  // FILE READING APPEARS TO WORK AS EXPECTED
  //
  //

  
  // CALCULATING CORRELATIONS AND VARIANCES

  vector< vector< vector<float> > > avgCorr(inputFiles.size());
  vector< vector< vector<float> > > var(inputFiles.size());
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      (avgCorr[i]).resize(NUMOFGENES);
      (var[i]).resize(NUMOFGENES);
      for (int j=0; j<NUMOFGENES; ++j)
	{
	  (avgCorr[i][j]).resize(NUMOFGENES);
	  (var[i][j]).resize(NUMOFGENES);
	}
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for corr/var\t" << section_time << "\t" << total_time << endl;  
  

  cout << "Before correlations...\n";

  /*
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      // Parallellized inside correlation functions
      rawToCorrMem_Eval(rawData[i], avgCorr[i], var[i], numOfMeas[i], NUMOFGENES, c, k, threadCount);
    }
  */

  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      runSpearmanCorrelation(rawData[i], avgCorr[i], numOfMeas[i], NUMOFGENES, threadCount);
      calculateVariance(rawData[i], avgCorr[i], var[i], numOfMeas[i], NUMOFGENES, k, terminationLimit, maxSubs, rSeed, threadCount);
    }



  
  cout << "After correlations...\n";

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "rawToCorr\t\t\t" << section_time << "\t" << total_time << endl;  
  
  cout << endl;
  cout << "Example value from avgCorr:\t" << avgCorr[0][0][1] << endl;
  cout << "Example value from var:    \t" << var[0][0][1] << endl;
  cout << endl;
  

  
  // Swap trick for deallocating vector memory
  vector< vector< vector<float> > >().swap(rawData);


  cout << "Calculating correlations and variances done...\n";
  
  // WRITING CORR AND VAR VALUES TO FILE UPON REQUEST
  if (o.corr)
    {
#     pragma omp parallel for num_threads(threadCount) schedule(static,1)      
      for (unsigned int i=0; i<inputFiles.size(); ++i)
	{
	  writeCorr(inputFiles[i], filteredGenes[i], avgCorr[i], var[i], NUMOFGENES);	  
	}
        end = clock();
	section_time = double(end-begin)/CLOCKS_PER_SEC;
	total_time += section_time;
	begin = end;
	fout << "writeCorr()\t\t\t" << section_time << "\t" << total_time << endl;
    }

  cout << "Before corr as output check...\n";

  // CHECKING IF USER REQUEST IS FULFILLED
  if (!o.unfiltered && !o.filtered && !o.all)
    { 
      exit(0);
    }

  cout << "After free memory sequence...\n";
  

  // CALCULATING CSD-SCORES
  
  cout << "number of inputFiles: " << inputFiles.size() << endl;
  int noc = numberOfCombinations(int(inputFiles.size()));
  cout << "number of combinations: " << noc << endl;
  vector< vector< vector< vector<float> > > > csd(noc);
  cout << "csd.size(): " << csd.size() << endl;
  
  for (int i=0; i<noc; ++i)
    {
      (csd[i]).resize(NUMOFGENES);
      for (int j=0; j<NUMOFGENES; ++j)
	{
	(csd[i][j]).resize(NUMOFGENES);	  
	  for (int l=0; l<NUMOFGENES; ++l)
	    {
	      (csd[i][j][l]).resize(3);
	    }
	}
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for csd\t" << section_time << "\t" << total_time << endl;

  
  cout << "After memory allocation...\n";

  
  int it = 0;
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      for (unsigned int j=i+1; j<inputFiles.size(); ++j)
	{
	  // Parallellized loop inside function
	  dependencyScores(avgCorr[i], avgCorr[j], var[i], var[j], csd[it], NUMOFGENES, threadCount);
	  ++it;
	}
    }

  if (o.all)
    {
      // NB! does not support multiple files as input yet, only writes results
      //     the first network (from file1 and file2)
      writeAll(avgCorr, var, csd, filteredGenes, NUMOFGENES, outPrefix);
    }


  // Freeing memory
  //vector< vector< vector<float> > >().swap(avgCorr);
  //vector< vector< vector<float> > >().swap(var);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "dependencyScores()\t\t" << section_time << "\t" << total_time << endl;  

  
  
  cout << "After calculating dependency scores...\n";

  cout << "Dependency scores for G0-G1:\n"
       << "C:\t" << csd[0][0][1][0] << endl
       << "S:\t" << csd[0][0][1][1] << endl
       << "D:\t" << csd[0][0][1][2] << endl;
  cout << "Dependency scores for G1-G2:\n"
       << "C:\t" << csd[0][1][2][0] << endl
       << "S:\t" << csd[0][1][2][1] << endl
       << "D:\t" << csd[0][1][2][2] << endl;
  cout << "Dependency scores for G2-G5:\n"
       << "C:\t" << csd[0][2][5][0] << endl
       << "S:\t" << csd[0][2][5][1] << endl
       << "D:\t" << csd[0][2][5][2] << endl;  
  
  // CHECKING IF USER REQUESTED UNFILTERED-NETWORKS TO BE WRITTEN TO FILE
  if (o.unfiltered)
    {
      // Parallellized inside function; one thread per network
      writeUnfiltered(csd, filteredGenes, NUMOFGENES, threadCount);
      
      end = clock();
      section_time = double(end-begin)/CLOCKS_PER_SEC;
      total_time += section_time;
      begin = end;
      fout << "writeUnfiltered()\t\t" << section_time << "\t" << total_time << endl;
    }

  
  // CHECKING IF USER REQUEST IS FULFILLED
  if (!o.filtered)
    {
      exit(0);
    }


  cout << "Before calculating threshold values...\n";

  
  // CALCULATING CSD-THRESHOLDS FOR ALL NETWORKS
  vector<float> cThresh(csd.size());
  vector<float> sThresh(csd.size());
  vector<float> dThresh(csd.size());

# pragma omp parallel for num_threads(threadCount)
  for (unsigned int i=0; i<csd.size(); ++i)
    {
      cout << "\tCalculating threshold values iteration: " << i+1 << "...\n";
      thresholdValuesWithReplacement(csd[i], NUMOFGENES, s, iSamples, cThresh[i], sThresh[i], dThresh[i], rSeed);
      cout << "\tCalculating threshold values (draw with replacement) iteration: " << i+1 << "...DONE\n";
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "thresholdValues()\t\t" << section_time << "\t" << total_time << endl;  
  

  cout << "After calculating threshold values...\n";

  cout << "c-thresh[0]: " << cThresh[0] << endl;
  cout << "s-thresh[0]: " << sThresh[0] << endl;
  cout << "d-thresh[0]: " << dThresh[0] << endl;
  cout << endl;
  //  cout << "c-thresh[1]: " << cThresh[1] << endl;
  //  cout << "s-thresh[1]: " << sThresh[1] << endl;
  //  cout << "d-thresh[1]: " << dThresh[1] << endl;
  
  
  // FILTERING CSD-NODES
  vector<CSD_Pair> csdNodes;

  //////////////////////////
  // Can this be parallellized somehow??
  // Possible: Can parallellize, but the csdNodes-vector will not be sorted by gene-ID anymore.
  //           This is okay, because it is not required for neighter cytoscape or function
  //           writing filter to file.
  filterCSDscores(csd, NUMOFGENES, cThresh, sThresh, dThresh, csdNodes, threadCount);
  ///////////////////////


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "filterCSDscores()\t\t" << section_time << "\t" << total_time << endl;

  
  cout << "After filtering CSD-networks...\n";


  cout << "Number of important correlations: " << csdNodes.size() << endl;  
  
  
  // WRITING FILTERED NETWORK TO FILE
  writeFiltered(csdNodes, filteredGenes, noc, outPrefix);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "writeFiltered()\t\t\t" << section_time << "\t" << total_time << endl;  

  
  //  vector< vector< vector< vector<float> > > >().swap(csd);

  cout << "End of program!\n";
  
}




void runProgramCorrInput(const vector<string>& inputFiles, const Output& o, const double& s, const int& iSamples, const int& threadCount, const string& pt_filename, const string& outPrefix)
{  
  // Setting up output file for performance test
  ofstream fout(pt_filename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed for performance test outfile.\n";
      exit(1);
    }  

  fout << "----------------------------------------------\n"
       << "              PERFORMANCE TEST\n"
       << "----------------------------------------------\n\n"
       << "Running times for program using corr/var-data as input.\n\n\n"
       << "FUNCTIONS/SECTION\t\tTIME (s)\tTOTAL TIME (s)\n"
       << "-------------------------------------------------------------------\n";
  
  clock_t begin = clock();
  
  
  // PREPROCESSING INPUT
  int numOfInputFiles = inputFiles.size();

  vector<int> numOfGenes(numOfInputFiles);
  for (int i=0; i<numOfInputFiles; ++i)
    {
      numOfGenes[i] = invNumberOfCombinations(getNumOfRows(inputFiles[i]));
    }

  clock_t end = clock();
  double section_time = double(end-begin)/CLOCKS_PER_SEC;
  double total_time = section_time;
  begin = end;
  fout << "determining numOfGenes[i]\t" << section_time << "\t" << total_time << endl;

  
  cout << "numOfGenes[0]: " << numOfGenes[0] << endl;
  cout << "numOfGenes[1]: " << numOfGenes[1] << endl;
  
  vector< vector< vector<float> > > avgCorr(numOfInputFiles);
  vector< vector< vector<float> > > var(numOfInputFiles);;
  vector< vector<Gene> > filteredGenes(numOfInputFiles);  
  for (int fi=0; fi<numOfInputFiles; ++fi)
    {
      (filteredGenes[fi]).resize(numOfGenes[fi]);
      (avgCorr[fi]).resize(numOfGenes[fi]);
      (var[fi]).resize(numOfGenes[fi]);
      for (int i=0; i<numOfGenes[fi]; ++i)
	{
	  (avgCorr[fi][i]).resize(numOfGenes[fi]);
	  (var[fi][i]).resize(numOfGenes[fi]);
	}
    }


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for corr/var\t" << section_time << "\t" << total_time << endl;    
  

  cout << "Before file reading...\n";

# pragma omp parallel for num_threads(threadCount)  
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      readCorr(inputFiles[i], filteredGenes[i], avgCorr[i], var[i], numOfGenes[i]);
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "readCorr()\t\t\t" << section_time << "\t" << total_time << endl;  
  
  cout << "After file reading...\n";
  /*cout << "Writing all corr and var values to screen for debugging file reading:\n";
  for (int fi=0; fi<numOfInputFiles; ++fi)
    {
      for (int i=0; i<numOfGenes[0]; ++i)
	{
	  for (int j=0; j<numOfGenes[0]; ++j)
	    {
	      cout << "\t" << avgCorr[fi][i][j];
	    }
	  cout << endl;
	}
      cout << endl;
    }

  cout << endl;
  for (int fi=0; fi<numOfInputFiles; ++fi)
    {
      for (int i=0; i<numOfGenes[0]; ++i)
	{
	  for (int j=0; j<numOfGenes[0]; ++j)
	    {
	      cout << "\t" << var[fi][i][j];
	    }
	  cout << endl;
	}
      cout << endl;
    }  
  */
  
  // CALCULATING CSD-SCORES
  cout << "number of inputFiles: " << inputFiles.size() << endl;
  int noc = numberOfCombinations(int(inputFiles.size()));
  cout << "number of combinations: " << noc << endl;
  vector< vector< vector< vector<float> > > > csd(noc);
  cout << "csd.size(): " << csd.size() << endl;
  
  for (int i=0; i<noc; ++i)
    {
      (csd[i]).resize(numOfGenes[0]);
      for (int j=0; j<numOfGenes[0]; ++j)
	{
	(csd[i][j]).resize(numOfGenes[0]);	  
	  for (int l=0; l<numOfGenes[0]; ++l)
	    {
	      (csd[i][j][l]).resize(3);
	    }
	}
    }


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for csd\t" << section_time << "\t" << total_time << endl;  

  
  cout << "After csd memory allocation...\n";

  
  int it = 0;
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      for (unsigned int j=i+1; j<inputFiles.size(); ++j)
	{
	  // Parallellized loop inside function
	  dependencyScores(avgCorr[i], avgCorr[j], var[i], var[j], csd[it], numOfGenes[0], threadCount);
	  ++it;
	}
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "dependencyScores()\t\t" << section_time << "\t" << total_time << endl;  

  
  vector< vector< vector<float> > >().swap(avgCorr);
  vector< vector< vector<float> > >().swap(var);

  cout << "After calculating dependency scores...\n";

  cout << "Dependency scores for G0:\n"
       << "C:\t" << csd[0][0][1][0] << endl
       << "S:\t" << csd[0][0][1][1] << endl
       << "D:\t" << csd[0][0][1][2] << endl;
  

  // CHECKING IF USER REQUESTED UNFILTERED-NETWORKS TO BE WRITTEN TO FILE
  if (o.unfiltered)
    {
      cout << "Writing unfiltered to file...\n";
      // Parallellized inside function; one thread per network
      writeUnfiltered(csd, filteredGenes, numOfGenes[0], threadCount);

      end = clock();
      section_time = double(end-begin)/CLOCKS_PER_SEC;
      total_time += section_time;
      begin = end;
      fout << "writeUnfiltered()\t\t" << section_time << "\t" << total_time << endl;      
    }

  
  // CHECKING IF USER REQUEST IS FULFILLED
  if (!o.filtered)
    {
      exit(0);
    }


  cout << "Before calculating threshold values...\n";

  
  // CALCULATING CSD-THRESHOLDS FOR ALL NETWORKS
  vector<float> cThresh(csd.size());
  vector<float> sThresh(csd.size());
  vector<float> dThresh(csd.size());
  
# pragma omp parallel for num_threads(threadCount)
  for (unsigned int i=0; i<csd.size(); ++i)
    {
      cout << "\tCalculating threshold values iteration: " << i+1 << "...\n";
      thresholdValues(csd[i], numOfGenes[0], s, iSamples, cThresh[i], sThresh[i], dThresh[i]);
      cout << "\tCalculating threshold values iteration: " << i+1 << "...DONE\n";
    }

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "thresholdValues()\t\t" << section_time << "\t" << total_time << endl;  
  

  cout << "After calculating threshold values...\n";

  cout << "c-thresh[0]: " << cThresh[0] << endl;
  cout << "s-thresh[0]: " << sThresh[0] << endl;
  cout << "d-thresh[0]: " << dThresh[0] << endl;
  cout << endl;
  cout << "c-thresh[1]: " << cThresh[1] << endl;
  cout << "s-thresh[1]: " << sThresh[1] << endl;
  cout << "d-thresh[1]: " << dThresh[1] << endl;
  
  
  // FILTERING CSD-NODES
  vector<CSD_Pair> csdNodes;
  // Parallellized loop inside function
  filterCSDscores(csd, numOfGenes[0], cThresh, sThresh, dThresh, csdNodes, threadCount);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "filterCSDscores()\t\t" << section_time << "\t" << total_time << endl;  
  

  cout << "After filtering CSD-networks...\n";


  cout << "Number of significantly correlated genes: " << csdNodes.size() << endl;  
  
  
  // WRITING FILTERED NETWORK TO FILE
  // No good way of parallellizing function
  writeFiltered(csdNodes, filteredGenes, noc, outPrefix);

  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "writeFiltered()\t\t\t" << section_time << "\t" << total_time << endl;  

  
  vector< vector< vector< vector<float> > > >().swap(csd);

  cout << "End of program!\n";
}













void runProgramUnfilteredInput(const vector<string>& inputFiles, const double& s, const int& iSamples, const int& threadCount, const string& pt_filename, const string& outPrefix)
{
  // Setting up output file for performance test
  ofstream fout(pt_filename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed for performance test outfile.\n";
      exit(1);
    }

  fout << "----------------------------------------------\n"
       << "              PERFORMANCE TEST\n"
       << "----------------------------------------------\n\n"
       << "Running times for program using unfiltered network as input.\n\n\n"
       << "FUNCTIONS/SECTION\t\tTIME (s)\tTOTAL TIME (s)\n"
       << "-------------------------------------------------------------------\n";
  
  clock_t begin = clock();

  

  int numOfInputFiles = inputFiles.size();

  vector<int> numOfGenes(numOfInputFiles);
  for (int i=0; i<numOfInputFiles; ++i)
    {
      numOfGenes[i] = invNumberOfCombinations(getNumOfRows(inputFiles[i]));
    }

  
  clock_t end = clock();
  double section_time = double(end-begin)/CLOCKS_PER_SEC;
  double total_time = section_time;
  begin = end;
  fout << "determining numOfGenes[i]\t" << section_time << "\t" << total_time << endl;
  
  
  cout << "numOfGenes[0]: " << numOfGenes[0] << endl;
  cout << "numOfGenes[1]: " << numOfGenes[1] << endl;
  

  vector< vector<Gene> > filteredGenes(numOfInputFiles);  
  vector< vector< vector< vector<float> > > > csd(numOfInputFiles);
  for (int fi=0; fi<numOfInputFiles; ++fi)
    {
      (filteredGenes[fi]).resize(numOfGenes[fi]);
      (csd[fi]).resize(numOfGenes[0]);
      for (int i=0; i<numOfGenes[0]; ++i)
	{
	  (csd[fi][i]).resize(numOfGenes[0]);
	  for (int j=0; j<numOfGenes[0]; ++j)
	    {
	      (csd[fi][i][j]).resize(3);
	    }
	}
    }


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "allocating memory for csd\t" << section_time << "\t" << total_time << endl;

  

  cout << "Before reading from file...\n";
  
# pragma omp parallel for num_threads(threadCount)
  for (unsigned int i=0; i<inputFiles.size(); ++i)
    {
      readUnfiltered(inputFiles[i], filteredGenes[i], csd[i], numOfGenes[0]);
    }


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "readUnfiltered()\t\t" << section_time << "\t" << total_time << endl;  

  
  cout << "After reading from file...\n";

  cout << "After reading dependency scores...\n";

  cout << "Dependency scores for G0:\n"
       << "C:\t" << csd[0][0][1][0] << endl
       << "S:\t" << csd[0][0][1][1] << endl
       << "D:\t" << csd[0][0][1][2] << endl;
  
  
  cout << "Before calculating threshold values...\n";

  
  // CALCULATING CSD-THRESHOLDS FOR ALL NETWORKS
  vector<float> cThresh(csd.size());
  vector<float> sThresh(csd.size());
  vector<float> dThresh(csd.size());

  for (unsigned int i=0; i<csd.size(); ++i)
    {
      cout << "\tCalculating threshold values iteration: " << i+1 << "...\n";
      thresholdValues(csd[i], numOfGenes[0], s, iSamples, cThresh[i], sThresh[i], dThresh[i]);
      cout << "\tCalculating threshold values iteration: " << i+1 << "...DONE\n";
    }


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "thresholdValues()\t\t" << section_time << "\t" << total_time << endl;  
  

  cout << "After calculating threshold values...\n";

  cout << "c-thresh[0]: " << cThresh[0] << endl;
  cout << "s-thresh[0]: " << sThresh[0] << endl;
  cout << "d-thresh[0]: " << dThresh[0] << endl;
  cout << endl;
  cout << "c-thresh[1]: " << cThresh[1] << endl;
  cout << "s-thresh[1]: " << sThresh[1] << endl;
  cout << "d-thresh[1]: " << dThresh[1] << endl;
  
  
  // FILTERING CSD-NODES
  vector<CSD_Pair> csdNodes;
  // Parallellized for loop inside function
  filterCSDscores(csd, numOfGenes[0], cThresh, sThresh, dThresh, csdNodes, threadCount);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "filterCSDscores()\t\t" << section_time << "\t" << total_time << endl;  
  

  cout << "After filtering CSD-networks...\n";


  cout << "Number of significantly correlated genes: " << csdNodes.size() << endl;  
  
  
  // WRITING FILTERED NETWORK TO FILE
  int noc = numberOfCombinations(int(inputFiles.size()));
  writeFiltered(csdNodes, filteredGenes, noc, outPrefix);


  end = clock();
  section_time = double(end-begin)/CLOCKS_PER_SEC;
  total_time += section_time;
  begin = end;
  fout << "writeFiltered()\t\t\t" << section_time << "\t" << total_time << endl;    

  
  vector< vector< vector< vector<float> > > >().swap(csd);

  cout << "End of program!\n";
  
}



/*
void rawToCorr(vector< vector<double> >& rawData, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& NUMOFMEAS, const int& NUMOFGENES, const string& corrMethod, const int& subsampleSize, const int& threadCount)
{
  
  // Reading input file
  const int NUMOFGENES = getNumOfCols(inFile)-1;
  const int NUMOFMEAS = getNumOfRows(inFile)-1;
  vector<string> mIDs;
  vector<string> gIDs;
  DblArr *rawData = new DblArr[NUMOFMEAS];
  for (int i=0; i<NUMOFMEAS; ++i) rawData[i] = new double[NUMOFGENES];
  readFile(inFile, NUMOFMEAS, NUMOFGENES, rawData, mIDs, gIDs);
  

  
  // Generating subsamples
  list< vector<int> > subsamples;
  //generateSubsamples(subsamples, NUMOFMEAS, subsampleSize);
  randomSubsampling(subsamples, NUMOFMEAS, subsampleSize, 1000);
  
  cout << "Number of subsamples: " << subsamples.size() << endl;
  
  // Allocating memory for correlation calculations
  const int NUMOFSUBSAMPLES = subsamples.size();
  FltArr **corrMatrix = new FltArr*[NUMOFGENES];
  for (int i=0; i<NUMOFGENES; ++i)
    {
      corrMatrix[i] = new FltArr[NUMOFGENES];
      for (int j=0; j<NUMOFGENES; ++j)
	{
	  corrMatrix[i][j] = new float[NUMOFSUBSAMPLES];
	}
    }


  // Running correlations based on corrMethod-parameter
  if (!corrMethod.compare("pearson"))
    {
      pearsonCorrelation(rawData, corrMatrix, subsamples, NUMOFGENES, subsampleSize, threadCount);
    }
  else if (!corrMethod.compare("spearman"))
    {
      spearmanCorrelation(rawData, corrMatrix, subsamples, NUMOFGENES, NUMOFMEAS, subsampleSize, threadCount);
    }
  else  if (!corrMethod.compare("bwm"))
    {
      cout << "Method not implemented yet\n";
      exit(0);
      // biWeightMidcorrelation();
    }
  else
    {
      cout << "Error: No valid correlation method entered\n";
      exit(1);
    }

  
  for (int l=0; l<NUMOFSUBSAMPLES; ++l)
    {
      for (int i=0; i<NUMOFGENES; ++i)
	{
	  for (int j=0; j<NUMOFGENES; ++j)
	    {
	      cout << "\t" << corrMatrix[i][j][l];
	    }
	  cout << endl;
	}
      cout << endl;
    }



  
  cout << "Test values after correlation:\n";
  cout << "Value from corrMatrix over diagonal: \t" << corrMatrix[0][1][0] << endl;
  cout << "Value from corrMatrix under diagonal:\t" << corrMatrix[1][0][0] << endl;
  
  
  // Calculating mean and variances
  meanAndVar(corrMatrix, avgCorr, var, NUMOFGENES, NUMOFSUBSAMPLES, threadCount);

  
  // Freeing memory used by corrMatrix as it is not needed anymore 
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=0; j<NUMOFGENES; ++j)
	{
	  delete[] corrMatrix[i][j];
	}
      delete[] corrMatrix[i];
    }
  delete[] corrMatrix;

}
*/
/*
void rawToCorrMem_Eval(vector< vector<double> >& rawData, vector< vector<double> >& avgCorr, vector< vector<double> >& var, const int& NUMOFMEAS, const int& NUMOFGENES, const string& corrMethod, const int& subsampleSize, const int& threadCount)
{
  
  // Reading input file
  const int NUMOFGENES = getNumOfCols(inFile)-1;
  const int NUMOFMEAS = getNumOfRows(inFile)-1;
  vector<string> mIDs;
  vector<string> gIDs;
  DblArr *rawData = new DblArr[NUMOFMEAS];
  for (int i=0; i<NUMOFMEAS; ++i) rawData[i] = new double[NUMOFGENES];
  readFile(inFile, NUMOFMEAS, NUMOFGENES, rawData, mIDs, gIDs);
  

  
  // Generating subsamples
  vector< vector<int> > subsamples;
  //generateSubsamples(subsamples, NUMOFMEAS, subsampleSize);
  randomSubsampling(subsamples, NUMOFMEAS, subsampleSize, 1000);
  
  cout << "Number of subsamples: " << subsamples.size() << endl;

  unsigned int numOfSubsamples = subsamples.size();

  cout << "NumOfSubsamples: " << subsamples.size() << endl;
  
  //vector<double> ssCorrs(subsamples.size());

  //cout << "ssCorrs.size():  " << ssCorrs.size() << endl;
  

  // Allocating memory for correlation calculations
  const int NUMOFSUBSAMPLES = subsamples.size();
  DblArr **corrMatrix = new DblArr*[NUMOFGENES];
  for (int i=0; i<NUMOFGENES; ++i)
    {
      corrMatrix[i] = new DblArr[NUMOFGENES];
      for (int j=0; j<NUMOFGENES; ++j)
	{
	  corrMatrix[i][j] = new double[NUMOFSUBSAMPLES];
	}
    }




  
  ///////////////////////////////////////////
  
  // This section calculates the spearman correlation. The section should be
  // rewritten after debugging and testing to accomodate the software's
  // support for other correlation measures.

  // The following parameters are needed as input for the "function"
  int numOfMeas = NUMOFMEAS;
  int numOfGenes = NUMOFGENES;
  //////

  
  // For determining the wall-time of the program:
  gettimeofday(&start_wt, NULL);
  
  
  if (!corrMethod.compare("spearman"))
    {
      cout << "Running Spearman correlations\n";
#     pragma omp parallel for num_threads(threadCount) 
      for (int i=0; i<numOfMeas; ++i)
	{
	  msRank(rawData[i], numOfGenes);
	}
   
      int l = 0;

      cout << "Done sorting rawData...\n";


      vector<vector<vector<int> > > subsampleList(threadCount);
      for (int tc=0; tc<threadCount; ++tc)
	{
	  (subsampleList[tc]).resize(numOfSubsamples);
	  for (unsigned int si=0; si<numOfSubsamples; ++si)
	    {
	      (subsampleList[tc][si]).resize(subsampleSize);
	      for (int ss=0; ss<subsampleSize; ++ss)
		{
		  subsampleList[tc][si][ss] = subsamples[si][ss];
		}
	    }
	}

      cout << "After allocating subsampleList\n";
      
#     pragma omp parallel num_threads(threadCount) private(l)
      {
	//vector<vector<int> > subsampleList(numOfSubsamples);
	int my_rank = omp_get_thread_num();

	//cout << "My_rank = " << my_rank << endl;
	//int numOfMeas = NUMOFMEAS;
	int numOfGenes = NUMOFGENES;
	unsigned int numOfSubsamples = subsamples.size();
	//ssCorrs.resize(subsampleSize);
	double* ssCorrs = new double[numOfSubsamples];

	double* x = new double[subsampleSize];
	double* y = new double[subsampleSize];
	
	//cout << "After ssCorrs.resize()\n";

	for (unsigned int si=0; si<numOfSubsamples; ++si)
	  {
	    for (int ss=0; ss<subsampleSize; ++ss)
	      {
		subsampleList[si][ss] = subsamples[si][ss];
	      }
	  }

	    
	//cout << "Before pragama omp for\n";
#       pragma omp for schedule(static,1)
	for (int i=0; i<numOfGenes; ++i)
	  {
	    //cout << "fredrik\n";
	    for (int j=i+1; j<numOfGenes; ++j)
	      {
		//cout << "magnus\n";
		l=0;
		
		//cout << "subsamples.size(): "<< subsamples.size() << endl;
	      
		avgCorr[i][j] = 0.0;
		var[i][j] = 0.0;
		//cout << "Gene pair: " << i << " " << j << endl;
		for (unsigned int it = 0; it < numOfSubsamples; ++it)
		  {
		    //cout << ".....\n";
		    
		    for (int k=0; k<subsampleSize; ++k)
		      {
			//cout << "...\n";
			x[k] = rawData[(subsampleList[my_rank][it][k])][i];
			y[k] = rawData[(subsampleList[my_rank][it][k])][j];
		      }
		    //cout << "Before pearson\n";
		    pearson(x,y,subsampleSize,ssCorrs[it]);
		    //cout << "ssCorrs[l] = " << ssCorrs[l] << endl;
		    avgCorr[i][j] += ssCorrs[it];
		    ++l;
		    //cout << "BLALALALA: this is for subsample " << l-1 << endl;
		  }
		//cout << "After calculating corr for all subsamples\n";
		avgCorr[i][j] = avgCorr[i][j]/double(numOfSubsamples);
		
		for (unsigned int s=0; s<numOfSubsamples; ++s)
		  {
		    var[i][j] += pow((ssCorrs[s] - avgCorr[i][j]),2);
		  }
		var[i][j] = var[i][j]/double(numOfSubsamples -1);
	      }
	  }
	//cout << "Before freeing memory\n";

	delete[] x;
	delete[] y;
	delete[] ssCorrs;	

	//cout << "After freeing memory\n";	
      }
    }
  else
    {
      cout << "Invalid corrlation method/Correlation method not supported.\n";
      exit(1);
    }


  // Determining wall-time
  gettimeofday(&end_wt, NULL);
  double delta = ((end_wt.tv_sec - start_wt.tv_sec)* 1000000u + end_wt.tv_usec - start_wt.tv_usec)/1.e6;

  cout << endl << "Wall-time of spearman-correlations:\t" << delta << endl << endl;
  
  cout << "Done calculating avgCorr and var\n";
  /////////////////////////////////////////////////





  // Running correlations based on corrMethod-parameter
  if (!corrMethod.compare("pearson"))
    {
      pearsonCorrelation(rawData, corrMatrix, subsamples, NUMOFGENES, subsampleSize, threadCount);
    }
  else if (!corrMethod.compare("spearman"))
    {
      spearmanCorrelation(rawData, corrMatrix, subsamples, NUMOFGENES, NUMOFMEAS, subsampleSize, threadCount);
    }
  else  if (!corrMethod.compare("bwm"))
    {
      cout << "Method not implemented yet\n";
      exit(0);
      // biWeightMidcorrelation();
    }
  else
    {
      cout << "Error: No valid correlation method entered\n";
      exit(1);
    }


  for (int l=0; l<NUMOFSUBSAMPLES; ++l)
    {
      for (int i=0; i<NUMOFGENES; ++i)
	{
	  for (int j=0; j<NUMOFGENES; ++j)
	    {
	      cout << "\t" << corrMatrix[i][j][l];
	    }
	  cout << endl;
	}
      cout << endl;
    }



  
  cout << "Test values after correlation:\n";
  cout << "Value from corrMatrix over diagonal: \t" << corrMatrix[0][1][0] << endl;
  cout << "Value from corrMatrix under diagonal:\t" << corrMatrix[1][0][0] << endl;

  
  // Calculating mean and variances
  //  meanAndVar(corrMatrix, avgCorr, var, NUMOFGENES, NUMOFSUBSAMPLES, threadCount);


  // Freeing memory used by corrMatrix as it is not needed anymore 
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=0; j<NUMOFGENES; ++j)
	{
	  delete[] corrMatrix[i][j];
	}
      delete[] corrMatrix[i];
    }
  delete[] corrMatrix;

}
*/

void runSpearmanCorrelation(vector<vector<float> >& rawData, vector<vector<float> >& avgCorr, const int& numOfMeas, const int& NUMOFGENES, const int& threadCount)
{
  // For determining the wall-time of the program:
  gettimeofday(&start_wt_corr, NULL);

  
# pragma omp parallel num_threads(threadCount)
  {
    int nom = numOfMeas;
    int nog = NUMOFGENES;
    
    double* x = new double[numOfMeas];
    double* y = new double[numOfMeas];

    double rho;
    
#   pragma omp for schedule(static,1)
    for (int g1=0; g1<nog; ++g1)
      {
	for (int g2=g1+1; g2<nog; ++g2)
	  {
	    for (int i=0; i<nom; ++i)
	      {
		x[i] = rawData[i][g1];
		y[i] = rawData[i][g2];
	      }

	    msRank(x, nom);
	    msRank(y, nom);
	    pearson(x,y,nom,rho);
	    avgCorr[g1][g2] = float(rho);
	  }
      }
    delete[] x;
    delete[] y;
  }

  // Determining wall-time
  gettimeofday(&end_wt_corr, NULL);
  double delta = ((end_wt_corr.tv_sec - start_wt_corr.tv_sec)* 1000000u + end_wt_corr.tv_usec - start_wt_corr.tv_usec)/1.e6;

  cout << "Wall-time of correlation calculations: " << delta << endl;
}



void calculateVariance(vector <vector<float> >& rawData, vector <vector<float> >& avgCorr, vector< vector<float> >& var, const int& numOfMeas, const int& NUMOFGENES, const int subsampleSize, const int& terminationLimit, const int& maxNumSubsamples, const int& randomSeed, const int& threadCount)
{

  
  vector< vector<int > > subsamples;

  randomSubsamplingWithLimit(subsamples, numOfMeas, subsampleSize, terminationLimit, maxNumSubsamples, randomSeed);

  unsigned int numOfSubsamples = subsamples.size();
  
  vector<vector<vector<int> > > subsampleList(threadCount);
  for (int tc=0; tc<threadCount; ++tc)
    {
      (subsampleList[tc]).resize(numOfSubsamples);
      for (unsigned int si=0; si<numOfSubsamples; ++si)
	{
	  (subsampleList[tc][si]).resize(subsampleSize);
	  for (int ss=0; ss<subsampleSize; ++ss)
	    {
	      subsampleList[tc][si][ss] = subsamples[si][ss];
	    }
	}
    }

  // For determining the wall-time of the program:
  gettimeofday(&start_wt_var, NULL);

  
# pragma omp parallel num_threads(threadCount)
  {
    int my_rank = omp_get_thread_num();

    int numOfGenes = NUMOFGENES;
    unsigned int numOfSubsamples = subsamples.size();
    double* ssCorrs = new double[numOfSubsamples];

    double* x = new double[subsampleSize];
    double* y = new double[subsampleSize];


#   pragma omp for schedule(static,1)
    for (int i=0; i<numOfGenes; ++i)
      {
	for (int j=i+1; j<numOfGenes; ++j)
	  {
	    var[i][j] = 0.0;
	    for (unsigned int it = 0; it < numOfSubsamples; ++it)
	      {
		for (int k=0; k<subsampleSize; ++k)
		  {
		    x[k] = rawData[(subsampleList[my_rank][it][k])][i];
		    y[k] = rawData[(subsampleList[my_rank][it][k])][j];
		  }
		spearman(x,y,subsampleSize,ssCorrs[it]);
	      }
	    for (unsigned int s=0; s<numOfSubsamples; ++s)
	      {
		var[i][j] += pow((ssCorrs[s] - avgCorr[i][j]),2);
	      }
	    var[i][j] = var[i][j]/float(numOfSubsamples -1);
	  }
      }
    delete[] x;
    delete[] y;
    delete[] ssCorrs;
  }

  // Determining wall-time
  gettimeofday(&end_wt_var, NULL);
  double delta = ((end_wt_var.tv_sec - start_wt_var.tv_sec)* 1000000u + end_wt_var.tv_usec - start_wt_var.tv_usec)/1.e6;

  cout << "Wall-time of variance calculations: " << delta << endl;
}

