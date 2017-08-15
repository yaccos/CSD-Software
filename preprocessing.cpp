#include <string>
#include <vector>
#include <algorithm>
#include <iostream>  // remove after debugging
#include <omp.h>

#include "preprocessing.h"
#include "utils.h"
#include "myStructs.h"


using namespace std;




void initializeGeneVectors(vector< vector<string> >& allIDs, vector< vector<Gene> >& allGenes, const int& threadCount)
{
  for (unsigned int i=0; i<allIDs.size(); ++i)
    {
      vector<Gene> dummy;
      allGenes.push_back(dummy);
      (allGenes[i]).resize((allIDs[i]).size());
    }
  
  //Gene dummy;
  //vector<Gene> vecDummy;
# pragma omp parallel for num_threads(threadCount)
  for (unsigned int fi=0; fi<allIDs.size(); ++fi)
    {
      //allGenes.push_back(vecDummy);

      //cout << endl << "\t\t\tInitializing allGenes vector " << i << endl;
      
      for (unsigned int gi=0; gi<(allIDs[fi]).size(); ++gi)
	{
	  //(allGenes[i]).push_back(dummy);
	  (allGenes[fi][gi]).geneID = allIDs[fi][gi];
	  (allGenes[fi][gi]).index = gi;
	  (allGenes[fi][gi]).sortedPos = -1;
	  (allGenes[fi][gi]).origin = fi;
	  //cout << endl << "\t\t\t\tgeneID = " << (allGenes[i][j]).geneID << endl;
	  //cout << "\t\t\t\tindex = " << (allGenes[i][j]).index << endl;
	  //cout << "\t\t\t\torigin = " << (allGenes[i][j]).origin << endl;
	}
    }
}


bool indicesInRange(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices)
{
  //cout << "\t\tCall to indices in range...\n";
  for (unsigned int i=0; i<allGenes.size(); ++i)
    {
      if (indices[i] < (allGenes[i]).size())
	{
	  continue;
	}
      else
	{
	  //cout << "\t\t\tReturn: FALSE\n";
	  return 0;
	}
    }
  //cout << "\t\t\tReturn TRUE\n";
  return 1;
}


bool allAreEqual(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices)
{
  //cout << "\t\tCall to allAreEqual..\n";
  for (unsigned int i=1; i<allGenes.size(); ++i)
    {
      if ((allGenes[i-1][indices[i-1]]).geneID == (allGenes[i][indices[i]]).geneID)
	{
	  continue;
	}
      else
	{
	  //cout << "\t\t\tReturn: FALSE\n";
	  return 0;
	}
    }
  //cout << "\t\t\tReturn: TRUE\n";
  return 1;
}


void incrementAllIndices(vector<unsigned int>& indices)
{
  for (unsigned int i=0; i<indices.size(); ++i)
    {
      indices[i] += 1;
    }
}


void filterAndSortGenes(vector< vector<string> >& allIDs, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const int& threadCount)
{
  //cout << "\t\tBeginning of filterAndSortGenes...\n";
  
  // Special case: only one input file
  //        - add all genes to filtered genes and return
  if (allIDs.size() == 1)
    {
      initializeGeneVectors(allIDs, filteredGenes, threadCount);
      return;
    }

  
  // Initializing gene vectors
  initializeGeneVectors(allIDs, allGenes, threadCount);
  //  cout << "After initialization of allGenes...\n";
  
  vector<Gene> dummy;
  for (unsigned int i=0; i<allIDs.size(); ++i) filteredGenes.push_back(dummy);
      
  
  //cout << "\t\tInitializing gene vectors done...\n";
  //cout << "\t\t\tSize of allIDs vector 1: " << (allIDs[0]).size() << endl;
  //cout << "\t\t\tSize of allIDs vector 2: " << (allIDs[1]).size() << endl;
  //cout << "\t\t\tSize of gene vector 1: " << (allGenes[0]).size() << endl;
  //cout << "\t\t\tSize of gene vector 2: " << (allGenes[1]).size() << endl;
  
  // Sorting genes

# pragma omp parallel for num_threads(threadCount)
  for (unsigned int i=0; i<allGenes.size(); ++i)
    {
      /*
      cout << "File " << i+1 << " before sorting by geneID:\n";
      for (unsigned int j=0; j<(allGenes[i]).size(); ++j)
	{
	  cout << (allGenes[i][j]).geneID << "\t" << (allGenes[i][j]).index << endl;
	}
      */
      //cout << "\t\tSorting genes iteration " << i+1 << endl;
      //sort((allGenes[i]).begin(), (allGenes[i]).end(), geneCompare); NB! geneCompare have been modified
      quickSortGenes(allGenes[i], 0, (allGenes[i]).size()-1, geneCompare);
      //cout << "\t\tSorting genes iteration " << i+1 << " done...\n";
      /*
      cout << "File " << i+1 << " after sorting by geneID:\n";
      for (unsigned int j=0; j<(allGenes[i]).size(); ++j)
	{
	  cout << (allGenes[i][j]).geneID << "\t" << (allGenes[i][j]).index << "\t" << (allGenes[i][j]).sortedPos << endl;
	}
      */
      //for (unsigned int j=0; j<(allGenes[i]).size(); ++j)
      //{
      //   cout << "\n\t\t\t\t" << (allGenes[i][j]).geneID;
      //}
    }







  
  // SO FAR ALL GOOD


  
  //cout << "\t\tSorting gene vectors done...\n";
  
  // Filtering Genes
  vector<unsigned int> indices(allGenes.size());
  for (unsigned int i=0; i<allGenes.size(); ++i) indices[i] = 0;

  //cout << "\tAfter indices initialization...\n";
  
  unsigned int i=0;
  Gene biggest;

  //cout << "\t\tPreparing extraction of genes done...\n";
  int j=0;
  while (indicesInRange(allGenes, indices))
    {
      if (allAreEqual(allGenes, indices))
	{
	  //cout << "\t\t\tPushing similar genes to vector...\n";
	  pushToFilteredGenes(allGenes, indices, filteredGenes);

	  for (unsigned int l=0; l<allIDs.size(); ++l)
	    {
	      (filteredGenes[l][j]).index = j;
	      (allGenes[l][indices[l]]).sortedPos = j;
	    }
	  ++j;
	  //cout << "\t\t\t\tPushing done...\n";
	  //cout << "\t\t\tIncrementing all indices...\n";
	  incrementAllIndices(indices);
	  //cout << "\t\t\t\tIncrementation done...\n";
	}
      else
	{
	  i = 0;
	  //cout << "\t\t\t Call to findBiggest()\n";
	  findBiggest(allGenes, indices, biggest);
	  while (!allAreEqual(allGenes, indices) && indicesInRange(allGenes,indices) && i<allGenes.size())
	    {
	      while (geneLessThan(allGenes[i][indices[i]], biggest) && indices[i]<(allGenes[i]).size())
		{
		  indices[i] += 1;
		}
	      if ((allGenes[i][indices[i]]).geneID == biggest.geneID)
		{
		  i++;
		  continue;
		}
	      else
		{
		  biggest.geneID = (allGenes[i][indices[i]]).geneID;
		  i = 0;
		  continue;
		}
	    }
	}
    }

  /*
  // Important step to maintain the wanted order of genes upon sorting them later
  for (unsigned int i=0; i<filteredGenes.size(); ++i)
    {
      for (unsigned int j=0; j<(filteredGenes[i]).size(); ++j)
	{
	  (filteredGenes[i][j]).index = j;
	}
    }
  */
  
}


void findBiggest(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices, Gene& biggest)
{
  biggest.geneID = "";
  for (unsigned int i=0; i<allGenes.size(); ++i)
    {
      if (geneLessThan(biggest, allGenes[i][indices[i]]))
	{
	  biggest.geneID = (allGenes[i][indices[i]]).geneID;
	}
    }
}


void pushToFilteredGenes(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices, vector< vector<Gene> >& filteredGenes)
{
  for (unsigned int i=0; i<allGenes.size(); ++i)
    {
      (filteredGenes[i]).push_back(allGenes[i][indices[i]]);
      //(allGenes[i][indices[i]]).sortedPos = &((filteredGenes[i]).back());
    }
}
