#include <iostream> //remove after debugging
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <omp.h>

#include "csdFiltering.h"
#include "myStructs.h"
#include "subsampling2.h"

using namespace std;

typedef double* DblPtr;
typedef bool* BoolPtr;


void dependencyScores(vector< vector<float> >& avgCorr1, vector< vector<float> >& avgCorr2, vector< vector<float> >& var1, vector< vector<float> >& var2, vector< vector< vector<float> > >& csd, const int& numOfGenes, const int& threadCount)
{
  double v,c;
# pragma omp parallel for num_threads(threadCount) private(c,v) schedule(static,1)
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i+1; j<numOfGenes; ++j)
	{
	  //cout << var1[i][j] << "\t" << var2[i][j] << endl;
	  //double bla = var1[i][j]+var2[i][j];
	  //cout << "bla: " << bla << endl;
	  v  = sqrt(abs(var1[i][j]+var2[i][j]));
	  //cout << "v: " << v << endl;
	  c = abs(avgCorr1[i][j]+avgCorr2[i][j]);
	  //cout << "c: " << c << endl;
	  //C-score
	  csd[i][j][0] = c/v;
	  //S-score
	  csd[i][j][1] = abs(abs(avgCorr1[i][j])-abs(avgCorr2[i][j]))/v;
	  //D-score
	  csd[i][j][2] = (abs(avgCorr1[i][j])+abs(avgCorr2[i][j])-c)/v;
	}
    }
}


void dependencyScores(double **avgCorr1, double **avgCorr2, double **var1, double **var2, double ***csd, const int& numOfGenes)
{
  double v,c;
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i; j<numOfGenes; ++j)
	{
	  v  = sqrt(var1[i][j]+var2[i][j]);
	  c = abs(avgCorr1[i][j]+avgCorr2[i][j]);
	  //C-score
	  csd[i][j][0] = c/v;
	  //S-score
	  csd[i][j][1] = abs(abs(avgCorr1[i][j])-abs(avgCorr2[i][j]))/v;
	  //D-score
	  csd[i][j][2] = (abs(avgCorr1[i][j])+abs(avgCorr2[i][j])-c)/v;
	}
    }
}




void thresholdValuesWithReplacement(vector< vector< vector<float> > >& csd, const int& numOfGenes, const double& significance, const int& numOfSamples, float& cThresh, float& sThresh, float& dThresh, int& rSeed)
{
  int subsampleSize = (1/significance)+1;

  if (rSeed == 0) srand(time(NULL));
  else srand(rSeed*3 +1);
  
  int r1, r2;

  int numDrawn = 0;
  int numOfGeneratedSamples = 0;
  
  double maxC, maxS, maxD;
  double cumMaxC = 0, cumMaxS = 0, cumMaxD = 0;

  while (numOfGeneratedSamples<numOfSamples)
    {
      maxC = 0;
      maxS = 0;
      maxD = 0;
      numDrawn = 0;
      while (numDrawn < subsampleSize)
	{
	  r1 = rand()%numOfGenes;
	  r2 = rand()%numOfGenes;
	  if (r1==r2) continue;
	  else if (r1<r2)
	    {
	      if (csd[r1][r2][0] > maxC) maxC = csd[r1][r2][0];
	      if (csd[r1][r2][1] > maxS) maxS = csd[r1][r2][1];
	      if (csd[r1][r2][2] > maxD) maxD = csd[r1][r2][2];
	      numDrawn++;
	    }
	  else if (r1>r2)
	    {
	      if (csd[r2][r1][0] > maxC) maxC = csd[r2][r1][0];
	      if (csd[r2][r1][1] > maxS) maxS = csd[r2][r1][1];
	      if (csd[r2][r1][2] > maxD) maxD = csd[r2][r1][2];
	      numDrawn++;
	    }
	}
      cumMaxC += maxC;
      cumMaxS += maxS;
      cumMaxD += maxD;
      numOfGeneratedSamples += 1;
    }

  cThresh = cumMaxC/numOfSamples;
  sThresh = cumMaxS/numOfSamples;
  dThresh = cumMaxD/numOfSamples;
}




void thresholdValues(vector< vector< vector<float> > >& csd, const int& numOfGenes, const double& significance, const int& numOfSamples, float& cThresh, float& sThresh, float& dThresh)
{
  //cout << "\t\tStart of thresholValues()\n";
  
  int subsampleSize = (1/significance)+1;
  int m = numOfSamples;  
  DblPtr *subsampleMax = new DblPtr[m];
  for (int i=0; i<m; ++i) subsampleMax[i] = new double[3];
  BoolPtr *isTaken = new BoolPtr[numOfGenes];
  for (int i=0; i<numOfGenes; ++i)
    {
      isTaken[i] = new bool[numOfGenes];
      for (int j=0; j<numOfGenes; ++j)
	{
	  if (i == j) isTaken[i][j] = 1;
	  else isTaken[i][j] = 0;
	}
    }

  //cout << "\t\tAfter allocating memory...\n";
  
  srand(time(NULL));
  int r1;
  int r2;
  bool taken;
  int index;
  for (int i=0; i<m; ++i)
    {
      index = 0;
      while (index < subsampleSize)
	{
	  taken = 1;
	  while (taken)
	    {
	      r1 = rand()%numOfGenes;
	      r2 = rand()%numOfGenes;
	      //cout << "\t\tDrawing random indices...\n";
	      if (isTaken[r1][r2] || isTaken[r2][r1] || (r1==r2)) continue;
	      else
		{
		  //cout << "\t\tEntering else statement...\n";
		  taken = 0;
		  isTaken[r1][r2] = 1;
		  isTaken[r2][r1] = 1;
		  if (r1 < r2)
		    {
		      //cout << "\t\tBefore r1 < r2...\n";
		      if (csd[r1][r2][0] > subsampleMax[i][0]) subsampleMax[i][0] = csd[r1][r2][0];
		      if (csd[r1][r2][1] > subsampleMax[i][1]) subsampleMax[i][1] = csd[r1][r2][1];
		      if (csd[r1][r2][2] > subsampleMax[i][2]) subsampleMax[i][2] = csd[r1][r2][2];
		      //cout << "\t\tAfter r1 < r2...\n";
		    }
		  else
		    {
		      //cout << "\t\tBefore r2 > r1...\n";
		      if (csd[r2][r1][0] > subsampleMax[i][0]) subsampleMax[i][0] = csd[r2][r1][0];
		      if (csd[r2][r1][1] > subsampleMax[i][1]) subsampleMax[i][1] = csd[r2][r1][1];
		      if (csd[r2][r1][2] > subsampleMax[i][2]) subsampleMax[i][2] = csd[r2][r1][2];
		      //cout << "\t\tBefore r2 > r1...\n";
		    }
		}
	    }
	  ++index;
	}
    }

  cThresh = 0;
  sThresh = 0;
  dThresh = 0;
  for (int i=0; i<m; ++i)
    {
      //cout << "\t\tBefore subsampleMax counting...\n";
      cThresh += subsampleMax[i][0];
      sThresh += subsampleMax[i][1];
      dThresh += subsampleMax[i][2];
      //cout << "\t\tAfter subsampleMax counting...\n";
    }
  cThresh = cThresh/m;
  sThresh = sThresh/m;
  dThresh = dThresh/m;

  //cout << "\t\tAfter calculating average max...\n";
  
  for (int i=0; i<m; ++i) delete[] subsampleMax[i];
  delete[] subsampleMax;
  for (int i=0; i<numOfGenes; ++i) delete[] isTaken[i];
  delete[] isTaken;

  //cout << "\t\tEnd of thresholdValues()\n";
}





void thresholdValues(double ***csd, const int& numOfGenes, const double& significance, const int& numOfSamples, double& cThresh, double& sThresh, double& dThresh)
{
  int subsampleSize = (1/significance)+1;
  int m = numOfSamples;  
  DblPtr *subsampleMax = new DblPtr[m];
  for (int i=0; i<m; ++i) subsampleMax[i] = new double[3];
  BoolPtr *isTaken = new BoolPtr[numOfGenes];
  for (int i=0; i<numOfGenes; ++i)
    {
      isTaken[i] = new bool[numOfGenes];
      for (int j=0; j<numOfGenes; ++j)
	{
	  if (i == j) isTaken[i][j] = 1;
	  else isTaken[i][j] = 0;
	}
    }
  srand(time(NULL));
  int r1;
  int r2;
  bool taken;
  int index;
  for (int i=0; i<m; ++i)
    {
      index = 0;
      while (index < subsampleSize)
	{
	  taken = 1;
	  while (taken)
	    {
	      r1 = rand()%numOfSamples;
	      r2 = rand()%numOfSamples;
	      if (isTaken[r1][r2] || isTaken[r2][r1]) continue;
	      else
		{
		  taken = 0;
		  isTaken[r1][r2] = 1;
		  isTaken[r2][r1] = 1;
		  if (r1 < r2)
		    {
		      if (csd[r1][r2][0] > subsampleMax[i][0]) subsampleMax[i][0] = csd[r1][r2][0];
		      if (csd[r1][r2][1] > subsampleMax[i][1]) subsampleMax[i][1] = csd[r1][r2][1];
		      if (csd[r1][r2][2] > subsampleMax[i][2]) subsampleMax[i][2] = csd[r1][r2][2];
		    }
		  else
		    {
		      if (csd[r2][r1][0] > subsampleMax[i][0]) subsampleMax[i][0] = csd[r2][r1][0];
		      if (csd[r2][r1][1] > subsampleMax[i][1]) subsampleMax[i][1] = csd[r2][r1][1];
		      if (csd[r2][r1][2] > subsampleMax[i][2]) subsampleMax[i][2] = csd[r2][r1][2];
		    }
		}
	    }
	  ++index;
	}
    }

  cThresh = 0;
  sThresh = 0;
  dThresh = 0;
  for (int i=0; i<m; ++i)
    {
      cThresh += subsampleMax[i][0];
      sThresh += subsampleMax[i][1];
      dThresh += subsampleMax[i][2];
    }
  cThresh = cThresh/m;
  sThresh = sThresh/m;
  dThresh = dThresh/m;

  for (int i=0; i<m; ++i) delete[] subsampleMax[i];
  delete[] subsampleMax;
  for (int i=0; i<numOfGenes; ++i) delete[] isTaken[i];
  delete[] isTaken;
}



void filterCSDscores(vector< vector< vector< vector<float> > > >& csd, const int& NUMOFGENES, const vector<float>& cThresh, const vector<float>& sThresh, const vector<float>& dThresh, vector<CSD_Pair>& csdNodes, const int& threadCount)
{
  CSD_Pair node(csd.size());

  //# pragma omp parallel for num_threads(threadCount) schedule(static, 1) private(node)
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=i+1; j<NUMOFGENES; ++j)
	{
	  for (unsigned int l=0; l<csd.size(); ++l)
	    {
	      if ( csd[l][i][j][0] > cThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 'c';
		}
	      else if ( csd[l][i][j][1] > sThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 's';
		}
	      else if ( csd[l][i][j][2] > dThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 'd';
		}
	      else
		{
		  node.type[l] = '-';
		}
	    }
	  if (node.activated)
	    {
	      csdNodes.push_back(node);
	      node.activated = 0;
	    }
	}
    }
}




void filterCSDscores(vector<double***>& csd, const int& NUMOFGENES, const vector<double>& cThresh, const vector<double>& sThresh, const vector<double>& dThresh, vector<CSD_Pair>& csdNodes)
{
  CSD_Pair node(csd.size());
  for (int i=0; i<NUMOFGENES; ++i)
    {
      for (int j=i+1; j<NUMOFGENES; ++j)
	{
	  for (unsigned int l=0; l<csd.size(); ++l)
	    {
	      if ( (csd[l])[i][j][0] > cThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 'c';
		}
	      else if ( (csd[l])[i][j][1] > sThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 's';
		}
	      else if ( (csd[l])[i][j][2] > dThresh[l])
		{
		  if (!node.activated)
		    {
		      node.gene_index1 = i;
		      node.gene_index2 = j;
		      node.activated = 1;
		    }
		  node.type[l] = 'd';
		}
	      else
		{
		  node.type[l] = '-';
		}
	    }
	  if (node.activated)
	    {
	      csdNodes.push_back(node);
	    }
	}
    }
}

void filterCSDscores(double ***csd, const int& numOfGenes, const double& cThresh, const double& sThresh, const double& dThresh, list<CSD_Pair>& csdNodes)
{
  CSD_Pair node;
  for (int i=0; i<numOfGenes; ++i)
    {
      for (int j=i; j<numOfGenes; ++j)
	{
	  if (csd[i][j][0]>cThresh)
	    {
	      node.gene_index1 = i;
	      node.gene_index2 = j;
	      node.type[0] = 'c';
	      csdNodes.push_back(node);
	    }
	  if (csd[i][j][1]>sThresh)
	    {
	      node.gene_index1 = i;	      
	      node.gene_index2 = j;
	      node.type[0] = 's';
	      csdNodes.push_back(node);
	    }
	  if (csd[i][j][2]>dThresh)
	    {
	      node.gene_index1 = i;	      
	      node.gene_index2 = j;
	      node.type[0] = 'd';
	      csdNodes.push_back(node);
	    }
	}
    }
}

