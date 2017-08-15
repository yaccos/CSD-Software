#include <iostream>
#include <string>
#include <fstream>
#include <list>
#include <cstdlib>
#include <vector>

#include "readFileTestUtils.h"

using namespace std;



void printFromFile(string filename)
{
  // Connecting file to stream
  ifstream fin;
  fin.open(filename.c_str());
  if (fin.fail())
  {
  cout << "Input file opening failed.";
  exit(1);
  }

  // Loop reading to screen
  string line;
  while (! fin.eof())
  {
    getline(fin, line);
    cout << line << endl;
  }
}


void printFromData(double **dataMatrix, int numOfMeas, int numOfGenes, vector<string>& measIDs, vector<string>& geneIDs)
{
  cout << "///\t";
  // Writes all gene-IDs
  vector<string>::const_iterator it;
  for (it = geneIDs.begin(); it != geneIDs.end(); ++it)
    {
      cout << *it  << "\t";
    }
  cout << endl;

  
  // Writes all measurment-IDs and data
  it = measIDs.begin();
  for (int i=0; i < numOfMeas; i++)
    {
      cout << *it << "\t"; 
      for (int j=0; j < numOfGenes; j++)
	{
	  cout << dataMatrix[i][j] << "\t";
	}
      ++it;
      cout << endl;
    }
}
  

void printStringVector(vector<string>& strList)
{
  vector<string>::const_iterator iterator;
  for (iterator = strList.begin(); iterator != strList.end(); ++iterator)
  {
    cout << *iterator << endl;
  }
}
