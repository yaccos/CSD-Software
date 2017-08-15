#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "preprocessing.h"
#include "fileIO.h"
#include "readFileTestUtils.h"
#include "utils.h"
#include "myStructs.h"


using namespace std;




void writeFiltered(vector<CSD_Pair>& csdNodes, vector< vector<Gene> >& filteredGenes, const int& noc, string outPrefix)
{
  
  quickSortGenes(filteredGenes[0], 0, int(filteredGenes.size()-1), indexCompare);
  
  string filename = outPrefix + "_filtered_network.txt";
  ofstream fout(filename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed in function writeFiltered().\n";
      exit(1);
    }

  for (unsigned int i=0; i<csdNodes.size(); ++i)
    {
      fout << (filteredGenes[0][(csdNodes[i]).gene_index1]).geneID << "\t"
	   << (filteredGenes[0][(csdNodes[i]).gene_index2]).geneID;
      for (int j=0; j<noc; ++j)
	{
	  fout << "\t" << (csdNodes[i]).type[j];
	}
      fout << endl;
    }

  if (fout.fail())
    {
      cout << "Error: fail-bit set before closing ofStream in function writeFiltered().\n";
      fout.close();
      exit(1);
    }

  fout.close();
}






void readRawFilesNew(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector< vector< vector<float> > >& rawData, const int& threadCount)
{
  //cout << "\t\tBeginning of readRawFiles()...\n";
  


  //cout << "\t\tOpening of mismatch-file done...\n";

  // Looping over all inputFiles and reading content
# pragma omp parallel for num_threads(threadCount) schedule(static,1)
  for (unsigned int fi=0; fi<inputFiles.size(); ++fi)
    {
      //cout << "\t\tReading from file number " << fi+1 << endl;
      
      int rank = omp_get_thread_num();

      // Connecting stream to output file for writing mismatches to.
      stringstream ss;
      ss << rank;
      string rank_str = ss.str();
      string outFilename = "mismatchedGenes" + rank_str + ".txt";
      ofstream fout(outFilename.c_str());
      if (fout.fail())
	{
	  cout << "Error: Output file opening failed in function readRawFiles().\n";
	  exit(1);
	}
      
      
      // Sorting filteredGenes by index
      // sort((filteredGenes[fi]).begin(), (filteredGenes[fi]).end(),indexCompare);
      //quickSortGenes(filteredGenes[fi], 0, (filteredGenes[fi]).size()-1, indexCompare);
      quickSortGenes(allGenes[fi], 0, (allGenes[fi]).size()-1, indexCompare);

      //cout << "Number of genes in allGenes[0]: " << (allGenes[0]).size() << endl;
      //cout << "Number of genes in allGenes[1]: " << (allGenes[1]).size() << endl;
      
      //cout << "\t\t\tSoring Genes by index done...\n";
      
      // Connecting input file to stream
      ifstream fin((inputFiles[fi]).c_str());
      if (fin.fail())
	{
	  cout << "Error: Input file opening failed in function readRawFiles().\n";
	  exit(1);
	}
      
      /*
      cout << "File " << fi+1 << " after sorting by index in fileIO:\n";
      for (unsigned int j=0; j<(allGenes[fi]).size(); ++j)
	{
	  cout << (allGenes[fi][j]).geneID << "\t" << (allGenes[fi][j]).index << "\t" << (allGenes[fi][j]).sortedPos << endl;
	}
      */
      
      //cout << "\t\t\tOpening of inputFile done...\n";
            
      // Reading from file, throwing away all measurements of genes not in "filteredGenes"
      // and write their IDs to the mismatch-file.
      if (vertical)
	{
	  //cout << "\t\t\t\tVertical input format flag set\n";
	  
	  string line, garbage, mismatch;
	  //double garb;
	  stringstream ss;
	  //int i = 0;
	  int gi = 0;
	  int j;
	  //int numOfGenes = (filteredGenes[fi]).size();

	  // Writing mismatches to file
	  int i=0;
	  gi=0;
	  while (i< int((allGenes[fi]).size()))
	    {
	      if ((allGenes[fi][i]).sortedPos == -1)
		{
		  fout << (allGenes[fi][i]).geneID << "\t" << inputFiles[fi] << endl;
		  ++i;
		}
	      else
		{
		  ++i;
		}
	    }
	  
	  // Removing first line
	  //for (int i=0; i<numOfMeas[fi]+1; ++i)
	  //  {
	  //    fin >> garbage;
	  //  }
	  unsigned int it=0;
	  getline(fin, line); // throwing away the top line containing measurement-IDs
	  while (!fin.eof() && it<(allGenes[fi]).size()  && !fin.fail())
	    {
	      //if (i == (filteredGenes[fi][gi]).index)
	      if ( !((allGenes[fi][gi]).sortedPos == -1) )
		{
		  //getline(fin,line);
		  //ss << line;
		  //ss >> garbage; // throwing away gene-ID at beginning of line
		  fin >> garbage;
		  j = 0;
		  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
		    {
		      //ss >> rawData[fi][j][(allGenes[fi][gi]).sortedPos];
		      fin >> rawData[fi][j][(allGenes[fi][gi]).sortedPos];		      
		      ++j;
		    }
		  //ss.str(""); //emptying stringstream
		  //++i;
		  ++gi;
		  ++it;
		}
	      else
		{
		  //fin >> mismatch;
		  //fout << mismatch << "\t" << inputFiles[fi] << endl;
		  //for (int i=0; i<numOfMeas[fi]; ++i)
		  //  {
		  //    fin >> garb;
		  //  }
		  getline(fin, line);
		  //ss << line;
		  //ss >> mismatch;
		  //fout << mismatch << "\t" << inputFiles[fi] << endl;
		  //ss.str(""); //emptying stringstream
		  //++i;
		  ++gi;
		  ++it;
		}
	    }
	}
      else
	{
	  //cout << "\t\t\t\tHorizontal input format flag set...\n";
	  
	  string line, garbage;
	  stringstream ss;
	  double dummy;
	  int i;
	  int gi;
	  int j;
	  //int j = 0;
	  //bool mismatchesRead = 0;
	  //int numOfGenes = (filteredGenes[fi]).size();

	  //Writing mismatches to file
	  //getline(fin, line);	  // throwing away the top line containing gene-IDs	  
	  //ss << line;
	  //ss >> garbage; // first element not a gene ID
	  fin >> garbage;
	  i=0;
	  gi=0;
	  while (i< int((allGenes[fi]).size()))
	    {
	      //ss >> dummy;
	      fin >> garbage;
	      //if (((allGenes[fi][i]).index == (filteredGenes[fi][gi]).index))
	      if ((allGenes[fi][i]).sortedPos == -1)
		{
		  //cout << "\t\t\t\tWriting to mismatch file...\n";
		  fout << (allGenes[fi][i]).geneID << "\t" << inputFiles[fi] << endl;
		  ++i;
		}
	      else
		{
		  //ss >> dummy;
		  ++i;
		}
	    }
	  j=0;
	  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
	    {
	      //getline(fin, line);
	      //ss << line;
	      //ss >> garbage; // throwing away measurement-ID at beginning of line
	      //cout << "\tReading to garbage\n";
	      fin >> garbage; // throwing away measurement-ID at beginning of line
	      i = 0;
	      gi = 0;
	      while (gi < int((allGenes[fi]).size()) && !fin.fail())
		{
		  fin >> dummy;
		  //if (i == (filteredGenes[fi][gi]).index)
		  if ( !((allGenes[fi][gi]).sortedPos == -1) )
		    {
		      //cout << "\t\t\t\t\tWriting to rawData...\n";
		      //ss >> rawData[fi][i][gi];
		      //++i
		      rawData[fi][j][(allGenes[fi][gi]).sortedPos] = dummy;
		      ++gi;
		    }
		  else
		    {
		      //    cout << "\t\t\t\t\tReading to mismatch-file...\n";
		      //    ss >> dummy;
		      //    fout << (allGenes[fi][i]).geneID << inputFiles[fi] << endl;
		      //++i;
		      //cout << "\tReading to rawData\n";
      		      //rawData[fi][j][(*((allGenes[fi][gi]).sortedPos)).index] = dummy;
		      //cout << "\tAfter reading to rawData\n";
		      ++gi;
		    }
		}
	      ++j;
	    }
	}

      if (fout.fail())
	{
	  cout << "Error: fail-bit set before closing fout in function readRawFiles()\n";
	  fout.close();
	  exit(1);
	}
      fout.close();
  
      if (fin.fail() && !fin.eof())
	{
	  cout << "Error: fail-bit set before closing fin in iteration number "
	       << fi << " in function readRawFiles()\n";
	  fin.close();
	  exit(1);
	}
      fin.close();
    }
}


















void readRawFiles(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector< vector< vector<double> > >& rawData)
{
  //cout << "\t\tBeginning of readRawFiles()...\n";
  
  // Connecting stream to output file for writing mismatches to.
  string outFilename = "mismatchedGenes.txt";
  ofstream fout(outFilename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed in function readFiles().\n";
      exit(1);
    }

  //cout << "\t\tOpening of mismatch-file done...\n";

  // Looping over all inputFiles and reading content
  for (unsigned int fi=0; fi<inputFiles.size(); ++fi)
    {
      //cout << "\t\tReading from file number " << fi+1 << endl;

      // Sorting filteredGenes by index
      // sort((filteredGenes[fi]).begin(), (filteredGenes[fi]).end(),indexCompare);
      quickSortGenes(filteredGenes[fi], 0, (filteredGenes[fi]).size()-1, indexCompare);
      quickSortGenes(allGenes[fi], 0, (allGenes[fi]).size()-1, indexCompare);

      //cout << "\t\t\tSoring Genes by index done...\n";
      
      // Connecting input file to stream
      ifstream fin((inputFiles[fi]).c_str());
      if (fin.fail())
	{
	  cout << "Error: Input file opening failed in function readFiles().\n";
	  exit(1);
	}

      //cout << "\t\t\tOpening of inputFile done...\n";
            
      // Reading from file, throwing away all measurements of genes not in "filteredGenes"
      // and write their IDs to the mismatch-file.
      if (vertical)
	{
	  //cout << "\t\t\t\tVertical input format flag set\n";
	  
	  string line, garbage, mismatch;
	  stringstream ss;
	  int i = 0;
	  int gi = 0;
	  int j;
	  int numOfGenes = (filteredGenes[fi]).size();
	  getline(fin, line); // throwing away the top line containing measurement-IDs
	  while (!fin.eof() && gi<numOfGenes && !fin.fail())
	    {
	      if (i == (filteredGenes[fi][gi]).index)
		{
		  getline(fin,line);
		  ss << line;
		  ss >> garbage; // throwing away gene-ID at beginning of line
		  j = 0;
		  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
		    {
		      ss >> rawData[fi][j][gi];
		      ++j;
		    }
		  ++i;
		  ++gi;
		}
	      else
		{
		  getline(fin, line);
		  ss << line;
		  ss >> mismatch;
		  fout << mismatch << "\t" << inputFiles[fi] << endl;
		  ++i;
		}
	    }
	}
      else
	{
	  //cout << "\t\t\t\tHorizontal input format flag set...\n";
	  
	  string line, garbage;
	  stringstream ss;
	  double dummy;
	  int i;
	  int gi;
	  int j;
	  //int j = 0;
	  //bool mismatchesRead = 0;
	  int numOfGenes = (filteredGenes[fi]).size();

	  //Writing mismatches to file
	  //getline(fin, line);	  // throwing away the top line containing gene-IDs	  
	  //ss << line;
	  //ss >> garbage; // first element not a gene ID
	  fin >> garbage;
	  i=0;
	  gi=0;
	  while (gi<numOfGenes && i< int((allGenes[fi]).size()))
	    {
	      //ss >> dummy;
	      fin >> garbage;
	      if (((allGenes[fi][i]).index == (filteredGenes[fi][gi]).index))
		{
		  ++i;
		  ++gi;
		}
	      else
		{
		  //cout << "\t\t\t\tWriting to mismatch file...\n";
		  //ss >> dummy;
		  fout << (allGenes[fi][i]).geneID << "\t" << inputFiles[fi] << endl;
		  ++i;
		}
	    }
	  j=0;
	  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
	    {
	      //getline(fin, line);
	      //ss << line;
	      //ss >> garbage; // throwing away measurement-ID at beginning of line
	      //cout << "\tReading to garbage\n";
	      fin >> garbage; // throwing away measurement-ID at beginning of line
	      i = 0;
	      gi = 0;
	      while (gi<numOfGenes && i < int((allGenes[fi]).size()) && !fin.fail())
		{
		  fin >> dummy;
		  if (i == (filteredGenes[fi][gi]).index)
		    {
		      //cout << "\t\t\t\t\tWriting to rawData...\n";
		      //ss >> rawData[fi][i][gi];
		      //cout << "\tReading to rawData\n";
		      rawData[fi][j][gi] = dummy;		      
		      ++i;
		      ++gi;
		    }
		  else
		    {
		      //    cout << "\t\t\t\t\tReading to mismatch-file...\n";
		      //    ss >> dummy;
		      //    fout << (allGenes[fi][i]).geneID << inputFiles[fi] << endl;
		      ++i;
		    }
		}
	      ++j;
	    }
	}

      if (fin.fail() && !fin.eof())
	{
	  cout << "Error: fail-bit set before closing fin in iteration number "
	       << fi << " in function readFile()\n";
	  fin.close();
	  exit(1);
	}
      fin.close();
    }

  if (fout.fail())
    {
      cout << "Error: fail-bit set before closing fout in function readFile()\n";
      fout.close();
      exit(1);
    }
  fout.close();
}      
  





void readRawFiles(const vector<string>& inputFiles, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const vector<int>& numOfMeas, const bool& vertical, vector<double**>& rawData)
{
  cout << "\t\tBeginning of readRawFiles()...\n";
  
  // Connecting stream to output file for writing mismatches to.
  string outFilename = "mismatchedGenes.txt";
  ofstream fout(outFilename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed in function readFiles().\n";
      exit(1);
    }

  cout << "\t\tOpening of mismatch-file done...\n";
  
  // Looping over all inputFiles and reading content
  for (unsigned int fi=0; fi<inputFiles.size(); ++fi)
    {
      cout << "\t\tReading from file number " << fi+1 << endl;
      
      // Sorting filteredGenes by index
      //sort((filteredGenes[fi]).begin(), (filteredGenes[fi]).end(),indexCompare);
      quickSortGenes(filteredGenes[fi], 0, (filteredGenes[fi]).size(), indexCompare);
      
      cout << "\t\t\tSoring Genes by index done...\n";
      
      // Connecting input file to stream
      ifstream fin((inputFiles[fi]).c_str());
      if (fin.fail())
	{
	  cout << "Error: Input file opening failed in function readFiles().\n";
	  exit(1);
	}

      cout << "\t\t\tOpening of inputFile done...\n";
            
      // Reading from file, throwing away all measurements of genes not in "filteredGenes"
      // and write their IDs to the mismatch-file.
      if (vertical)
	{
	  cout << "\t\t\t\tVertical input format flag set\n";
	  
	  string line, garbage, mismatch;
	  stringstream ss;
	  int i = 0;
	  int gi = 0;
	  int j;
	  int numOfGenes = (filteredGenes[fi]).size();
	  getline(fin, line); // throwing away the top line containing measurement-IDs
	  while (!fin.eof() && gi<numOfGenes && !fin.fail())
	    {
	      if (i == (filteredGenes[fi][gi]).index)
		{
		  getline(fin,line);
		  ss << line;
		  ss >> garbage; // throwing away gene-ID at beginning of line
		  j = 0;
		  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
		    {
		      ss >> (rawData[fi])[j][gi];
		      ++j;
		    }
		  ++i;
		  ++gi;
		}
	      else
		{
		  getline(fin, line);
		  ss << line;
		  ss >> mismatch;
		  fout << mismatch << "\t" << inputFiles[fi] << endl;
		  ++i;
		}
	    }
	}
      else
	{
	  cout << "\t\t\t\tHorizontal input format flag set...\n";
	  
	  string line, garbage;
	  stringstream ss;
	  double dummy;
	  int i;
	  int gi;
	  int j = 0;
	  int numOfGenes = (filteredGenes[fi]).size();
	  getline(fin, line); // throwing away the top line containing gene-IDs	  
	  while (!fin.eof() && j<numOfMeas[fi] && !fin.fail())
	    {
	      getline(fin, line);
	      ss << line;
	      ss >> garbage; // throwing away measurement-ID at beginning of line
	      i = 0;
	      gi = 0;
	      while (gi<numOfGenes && !fin.fail())
		{
		  if (i == (filteredGenes[fi][gi]).index)
		    {
		      cout << "\t\t\t\t\tReading to rawData...\n";
		      ss >> (rawData[fi])[j][gi];
		      ++i;
		      ++gi;
		    }
		  else
		    {
		      cout << "\t\t\t\t\tReading to mismatch-file...\n";
		      ss >> dummy;
		      fout << (allGenes[fi][i]).geneID << inputFiles[fi] << endl;
		      ++i;
		    }
		}
	    }
	}

      if (fin.fail())
	{
	  cout << "Error: fail-bit set before closing fin in iteration number "
	       << fi << " in function readFile()\n";
	  fin.close();
	  exit(1);
	}
      fin.close();
    }

  if (fout.fail())
    {
      cout << "Error: fail-bit set before closing fout in function readFile()\n";
      fout.close();
      exit(1);
    }
  fout.close();
}



void writeCSDtoFile(const vector< vector<Gene> >& filteredGenes, const vector<CSD_Pair>& csdNodes, const unsigned int csdSize)
{
  string outFilename = "filteredCSD.txt";
  ofstream fout(outFilename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed in function \"writeCSDtoFile()\".\n";
      exit(1);
    }

  for (unsigned int i=0; i<csdNodes.size(); ++i)
    {	
      fout << (filteredGenes[0][(csdNodes[i]).gene_index1]).geneID << "\t"
	   << (filteredGenes[0][(csdNodes[i]).gene_index2]).geneID;
      for (unsigned int j=0; j<csdSize; ++j)
	{
	  fout << "\t" << (csdNodes[i]).type[j];
	}
      fout << endl;
    }

  if (fout.fail())
    {
      cout << "Error: fail-bit set before closing fout in function \"writeCSDtoFile()\"\n";
      fout.close();
      exit(1);
    }
  fout.close();
}

/*
void readFile(const string filename, const int& NUMOFMEAS, const int& NUMOFGENES, double **matrix, vector<string>& mIDs, vector<string>& gIDs)
{
  // Connecting file to stream
  ifstream fin(filename.c_str());
  if (fin.fail())
    {
      cout << "Input file opening failed.\n";
      exit(1);
    }

  // Initializing variables
  string line, ID, garbage;
  stringstream ss;

  // Reading gene-IDs
  getline(fin, line);
  ss << line;
  ss >> garbage;
  while (! ss.eof())
    {
      ss >> ID;
      gIDs.push_back(ID);
    }

  // Reading measurement-IDs and measurements
  ss.str("");
  int i = 0;
  int j;
  while (!fin.eof() && i < NUMOFMEAS && !fin.fail())
    {
      fin >> ID;
      mIDs.push_back(ID);
      j = 0;
      while (!fin.eof() && j < NUMOFGENES && !fin.fail())
	{
	  fin >> matrix[i][j];
	  j++;
	}
      i++;
    }

    if (fin.fail())
      {
        cout << "Error: fail-bit set in function readFile()\n";
	fin.close();
	exit(1);
      }

  fin.close();
}
*/



void writeAll(vector<vector<vector<float> > >& avgCorr, vector<vector<vector<float> > >& var, vector<vector<vector<vector<float> > > >& csd, const vector< vector<Gene> >& filteredGenes, const int& numOfGenes, string outPrefix)
{
  // Initializing out-stream
  string outFilename = outPrefix + "_all_data.txt";
  ofstream fout(outFilename.c_str());
  if (fout.fail())
    {
      cout << "Error: Output file opening failed in function \"writeCSDtoFile()\".\n";
      exit(1);
    }
  
  
  // Writing first two lines of file
  fout << "Gene 1:\tGene 2:\tavg_corr 1:\tavg_corr 2:\tvariance 1:\t variance 2:\tC-score:\t"
       << "S-score:\tD-score:\n"
       << "---------------------------------------------------------------------------------"
       << "------------------------\n";
  
  // Iterating over all combinations of genes and writing data to file
  for (int g1=0; g1<numOfGenes; ++g1)
    {
      for (int g2=g1+1; g2<numOfGenes; ++g2)
	{
	  fout << (filteredGenes[0][g1]).geneID << "\t" << (filteredGenes[0][g2]).geneID << "\t"
	       << avgCorr[0][g1][g2] << "\t" << avgCorr[1][g1][g2] << "\t"
	       << var[0][g1][g2] << "\t" << var[1][g1][g2] << "\t"
	       << csd[0][g1][g2][0] << "\t"
	       << csd[0][g1][g2][1] << "\t"
	       << csd[0][g1][g2][2] << "\n";
	}
    }
      
  // Checking for fail-bit flag
  if (fout.fail())
    {
      cout << "Error: fail-bit set before closing fout in function \"writeCSDtoFile()\"\n";
      fout.close();
      exit(1);
    }
  fout.close();
}








void readCorr(const string& infile, vector<Gene>& filteredGenes, vector< vector<float> >& avgCorr, vector< vector<float> >& var, int& numOfGenes)
{
  ifstream fin(infile.c_str());
  if (fin.fail())
    {
      cout << "Error: Input file opening failed in function readCorr().\n";
      exit(1);
    }

  //int numOfRows = getNumOfRows(infile);
  //numOfGenes = invNumberOfCombinations(numOfRows);
  
  string garbage, ID;
  Gene g;
  
  //Reading first line
  fin >> ID;
  g.geneID = ID;
  g.index = 0;
  filteredGenes[0] = g;
  fin >> ID;
  g.geneID = ID;
  g.index = 1;
  filteredGenes[1] = g;
  fin >> avgCorr[0][1];
  fin >> var[0][1];

  //Reading lines until all genes are included in filteredGenes
  for (int i=2; i<numOfGenes; ++i)
    {
      fin >> garbage;
      fin >> ID;
      g.geneID = ID;
      g.index = i;
      filteredGenes[i] = g;
      fin >> avgCorr[0][i];
      fin >> var[0][i];
    }

  // Reading rest of input file measurements
  int n = numOfGenes-2;
  int i = 1;
  while (n>0 && !fin.fail() && !fin.eof())
    {
      for (int j=numOfGenes-n; j<numOfGenes; ++j)
	{
	  fin >> garbage;
	  fin >> garbage;
	  fin >> avgCorr[i][j];
	  fin >> var[i][j];
	}
      ++i;
      --n;
    }

  if (fin.fail())
    {
      cout << "Error: fail-bit set before closing fin in function readCorr().\n";
      fin.close();
      exit(1);
    }
  
  fin.close();
}



  
  
  

void writeCorr(const string& inFile, const vector<Gene>& filteredGenes, vector< vector<float> >& avgCorr, vector< vector<float> >& var, const int& NUMOFGENES)
{
  string outFile = "Corr_Var_" + inFile;
  ofstream fout(outFile.c_str());
  if (fout.fail())
    {
      cout << "Output file opening failed.\n";
      exit(1);
    }
  
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=i+1; j<NUMOFGENES; ++j)
	{
	  fout << (filteredGenes[i]).geneID << "\t" << (filteredGenes[j]).geneID << "\t"
	       << avgCorr[i][j] << "\t" << var[i][j] << endl;
	}
    }
  if (fout.fail())
    {
      cout << "Error: fail-bit set in function writeCorr()\n";
      fout.close();
      exit(1);
    }

  fout.close();	  
}









void writeCorr(const string& inFile, const vector<Gene>& filteredGenes, double **avgCorr, double **var, const int& NUMOFGENES)
{
  string outFile = "Corr_Var_" + inFile;
  ofstream fout(outFile.c_str());
  if (fout.fail())
    {
      cout << "Output file opening failed.\n";
      exit(1);
    }
  
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=i+1; j<NUMOFGENES; ++j)
	{
	  fout << (filteredGenes[i]).geneID << "\t" << (filteredGenes[j]).geneID << "\t"
	       << avgCorr[i][j] << "\t" << var[i][j] << endl;
	}
    }
  if (fout.fail())
    {
      cout << "Error: fail-bit set in function writeCorr()\n";
      fout.close();
      exit(1);
    }

  fout.close();	  
}


void readUnfiltered(const string& infile, vector<Gene>& filteredGenes, vector< vector< vector<float> > >& csd, int& numOfGenes)
{
  ifstream fin(infile.c_str());
  if (fin.fail())
    {
      cout << "Error: Input file opening failed in function readUnfiltered().\n";
      exit(1);
    }

  int numOfRows = getNumOfRows(infile);
  numOfGenes = invNumberOfCombinations(numOfRows);

  csd.resize(numOfGenes);
  for (int i=0; i<numOfGenes; ++i)
    {
      (csd[i]).resize(numOfGenes);
      for (int j=0; j<numOfGenes; ++j)
	{
	  (csd[i][j]).resize(3);
	}
    }

  string garbage, ID;
  Gene g;

  //Reading first line
  fin >> ID;
  g.geneID = ID;
  g.index = 0;
  filteredGenes[0] = g;
  fin >> ID;
  g.geneID = ID;
  g.index = 1;
  filteredGenes[1] = g;
  fin >> csd[0][1][0];
  fin >> csd[0][1][1];
  fin >> csd[0][1][2];

  //Reading lines until all genes are included in filteredGenes
  for (int i=2; i<numOfGenes; ++i)
    {
      fin >> garbage;
      fin >> ID;
      g.geneID = ID;
      g.index = i;
      filteredGenes[i] = g;
      fin >> csd[0][i][0];
      fin >> csd[0][i][1];
      fin >> csd[0][i][2];
    }

  //Reading rest of input file measruements
  int n = numOfGenes-2;
  int i = 1;
  while (n>0 && !fin.fail() && fin.eof())
    {
      for (int j=numOfGenes-n; j<numOfGenes; ++j)
	{
	  fin >> garbage;
	  fin >> garbage;
	  fin >> csd[i][j][0];
	  fin >> csd[i][j][1];
	  fin >> csd[i][j][2];
	}
      ++i;
      --n;
    }

  if (fin.fail())
    {
      cout << "Error: fail-bit set before closing fin in function readUnfiltered().\n";
      fin.close();
      exit(1);
    }

  fin.close();
}



void writeUnfiltered(const vector< vector< vector< vector<float> > > >& csd, const vector< vector<Gene> >& filteredGenes, const int& numOfGenes, const int& threadCount)
{
# pragma omp parallel for num_threads(threadCount)  
  for (unsigned int fi=0; fi<csd.size(); ++fi)
    {
      ostringstream oss;
      oss << fi;
      string index = oss.str();
      string outfile = "unfiltered_network_"+index+".txt";
      ofstream fout(outfile.c_str());

      cout << "Before for loop\n";
      for (int i=0; i<numOfGenes; ++i)
	{
	  for (int j=i+1; j<numOfGenes; ++j)
	    {
	      fout << (filteredGenes[fi][i]).geneID << "\t"
		   << (filteredGenes[fi][j]).geneID << "\t"
		   << csd[fi][i][j][0] << "\t"
		   << csd[fi][i][j][1] << "\t"
		   << csd[fi][i][j][2] << "\n";
	    }
	}
      if (fout.fail())
	{
	  cout << "Error: fail-bit set before closing fout in for file " << fi
	       << " function writeUnfiltered().\n";
	  fout.close();
	  exit(1);
	}

      fout.close();
    }
}










void readAllGeneIDs(const vector<string>& inputFiles, const bool& vertical, vector< vector<string> >& allIDs, const int& threadCount)
{
  allIDs.resize(int(inputFiles.size()));
  if (vertical)
    {
#     pragma omp parallel for num_threads(threadCount)      
      for (unsigned int i=0; i<inputFiles.size(); ++i)
	{
	  int n = getNumOfRows(inputFiles[i]);
	  vector<string> dummy(n-1);
	  getRowIDs(inputFiles[i], dummy, n);
	  allIDs[i] = dummy;
	}
    }
  else
    {
#     pragma omp parallel for num_threads(threadCount)      
      for (unsigned int i=0; i<inputFiles.size(); ++i)
	{
	  int n = getNumOfCols(inputFiles[i]);
	  vector<string> dummy(n-1);
	  getColIDs(inputFiles[i], dummy, n);
	  allIDs[i] = dummy;
	}
    }
}
  

int getNumOfCols(string filename)
{
  // Connecting file to stream
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail())
  {
    cout << "Input file opening failed.";
    exit(1);
  }

  string line, dummy;
  stringstream ss;
  int numOfCols = 0;
  getline(fin, line);
  ss << line;
  
  while (ss >> dummy)
    {
      numOfCols += 1;
    }

  fin.close();

  return numOfCols;
}


int getNumOfRows(string filename)
{
  // Connecting file to stream
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail())
  {
    cout << "Input file opening failed.";
    exit(1);
  }

  string line;
  double numOfRows = 0;

  while (!fin.eof() && !fin.fail())
    {
      getline(fin, line);
      if (line.empty()) continue;
      numOfRows += 1;
    }

  fin.close();

  return numOfRows;
}



void getColIDs(const string& filename, vector<string>& ids, const int& NUMOFCOLS)
{
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail())
    {
      cout << "Error: Input file opening failed in function getColIDs.\n";
      exit(1);
    }

  // Initializing variables
  string line, ID, garbage;
  stringstream ss;

  // Reading gene-IDs
  getline(fin, line);
  ss << line;
  ss >> garbage; // remove first element, as this is expected to be "///"
  int i = 0;
  while (! ss.eof() && i < (NUMOFCOLS-1))
    {
      ss >> ID;
      ids[i] = ID;
      ++i;
    }

  if (fin.fail())
    {
      cout << "Error: fail-bit set while writing from file in function getColIDs.\n";
      exit(1);
    }
  fin.close();
}



void getRowIDs(const string& filename, vector<string>& ids, const int& NUMOFROWS)
{
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail())
    {
      cout << "Error: Input file opening failed in function getColIDs.\n";
      exit(1);
    }

  string line, ID;
  stringstream ss;
  getline(fin, line); // Throwing away first line
  
  // Reading measurement-IDs
  int i = 0;
  while (!fin.eof() && i < (NUMOFROWS-1) && !fin.fail())
    {
      getline(fin, line);
      ss << line;
      ss >> ID;
      ids[i] = ID;
      ss.str("");
      ++i;
    }

  if (fin.fail())
    {
      cout << "Error: fail-bit set while reading from file in function getRowIDs.\n";
      exit(1);
    }
  fin.close();
}


/*
void filterAndSortGenes(vector<string>& gIDs1, vector<string>& gIDs2, vector<string>& sortedGeneSequence, vector<int> geneIndices1, vector<int> geneIndices2)
{
  // Matching genes from the two datasets
  const int NUMOFGENES1 = gIDs1.size();
  const int NUMOFGENES2 = gIDs2.size();

  vector<int> mismatchIndices;

  int j;
  
  for (int i=0; i<NUMOFGENES1; ++i)
    {
      j=0;
      while (j<NUMOFGENES2)
	{
	  if (gIDs1[i] == gIDs2[j])
	    {
	      sortedGeneSequence.push_back(gIDs1[i]);
	      geneIndices1.push_back(i);
	      geneIndices2.push_back(j);
	      gIDs2[j] = "---";
	      break;
	    }
	  else
	    {
	      ++j;
	    }	  
	}
      if (j == NUMOFGENES2) mismatchIndices.push_back(i);
      j=0;
    }

  //Writing mismatches to file
  string filename = "mismatched_genes.txt";
  
  ofstream fout(filename.c_str());
  if (fout.fail())
    {
      cout << "Error: Input file opening failed.\n";
      exit(1);
    }

  fout << "Mismatched genes from file 1:\n";
  for (unsigned int i=0; i<mismatchIndices.size(); ++i)
    {
      fout << gIDs1[mismatchIndices[i]] << "\n";
    }

  fout << "\n\nMismatched genes from file 2:\n";
  for (int i=0; i<NUMOFGENES2; ++i)
    {
      if (gIDs2[i].compare("---"))
	{
	  fout << gIDs2[i] << "\n";
	}
    }

  if (fout.fail())
    {
      cout << "Error: fail-bit set while writing mismathces to file.\n";
      exit(1);
    }
  fout.close();    
}
*/
  
