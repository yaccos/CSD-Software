#include <vector>
#include <list>
#include <cstdlib>
#include <ctime>
#include <iostream>  //remove after debugging
#include <algorithm>

#include "subsampling2.h"
//#include "testUtils.h"   //remove after debugging

using namespace std;




void generateSubsamples(list< vector<int> >& subsamples, const int& n, const int& k)
{
  vector< vector<int> > domains;
  initializeDomains(domains, n);

  generateFirstSamples(subsamples, domains, n, k);

  int var = 0;
  bool success;
  vector<int> currentDomain(n);
  for (int i=0; i<n; ++i) currentDomain[i] = i;

  while (var < n)
    {
      vector<int> subsample;

      success = backtrackingSearch(domains, currentDomain, subsample, n, k, var);

      if (success)
	{
	  updateDomains(domains, subsample);
	  subsamples.push_back(subsample);
	  continue;
	}
      else
	{
	  currentDomain.erase(currentDomain.begin());
	  ++var;
	}
    }
}


void generateSubsamplesNoInitialSampling(list< vector<int> >& subsamples, const int& n, const int& k)
{
  vector< vector<int> > domains;
  initializeDomains(domains, n);
  
  int var = 0;
  bool success;
  vector<int> currentDomain(n);
  for (int i=0; i<n; ++i) currentDomain[i] = i;

  while (var < n)
    {
      vector<int> subsample;

      success = backtrackingSearch(domains, currentDomain, subsample, n, k, var);

      if (success)
	{
	  updateDomains(domains, subsample);
	  subsamples.push_back(subsample);
	  continue;
	}
      else
	{
	  currentDomain.erase(currentDomain.begin());
	  ++var;	  
	}
    }
}


void randomSubsampling(list< vector<int> >& subsamples, const int& n, const int& k, const int& terminationLimit)
{
  srand(time(NULL));
  vector< vector<int> > domains;
  initializeDomains(domains, n);

  vector<int> currentDomain(n);
  for (int i=0; i<n; ++i) currentDomain[i] = i;
  int fail = 0;
  bool success;

  while (fail<terminationLimit)
    {
      vector<int> subsample;

      success = randomDraw(domains, currentDomain, subsample, n, k);

      if(success)
	{
	  updateDomains(domains, subsample);
	  sort(subsample.begin(),subsample.end());  //subsamples must be in sorted order 
	  subsamples.push_back(subsample);
	  continue;
	}
      else
	{
	  ++fail;
	}
    }
}











void randomSubsampling(vector< vector<int> >& subsamples, const int& n, const int& k, const int& terminationLimit)
{
  srand(time(NULL));
  vector< vector<int> > domains;
  initializeDomains(domains, n);

  vector<int> currentDomain(n);
  for (int i=0; i<n; ++i) currentDomain[i] = i;
  int fail = 0;
  bool success;

  while (fail<terminationLimit)
    {
      vector<int> subsample;

      success = randomDraw(domains, currentDomain, subsample, n, k);

      if(success)
	{
	  updateDomains(domains, subsample);
	  subsamples.push_back(subsample);
	  continue;
	}
      else
	{
	  ++fail;
	}
    }
}



void randomSubsamplingWithLimit(vector< vector<int> >& subsamples, const int& n, const int& k, const int& failLimit, const int& maxNumSubsamples, const int& randomSeed)
{
  if (randomSeed>0) srand(randomSeed);
  else srand(time(NULL));


  vector< vector<int> > domains;
  initializeDomains(domains, n);

  
  vector<int> currentDomain(n);
  for (int i=0; i<n; ++i) currentDomain[i] = i;
  int fails = 0;
  int successes = 0;
  bool success;

  
  while (fails<failLimit && successes<maxNumSubsamples)
    {
      vector<int> subsample;

      success = randomDraw(domains, currentDomain, subsample, n, k);

      if(success)
	{
	  updateDomains(domains, subsample);
	  sort(subsample.begin(),subsample.end());  //subsamples must be in sorted order
	  subsamples.push_back(subsample);
	  ++successes;
	}
      else
	{
	  ++fails;
	}
    }
}




bool backtrackingSearch(const vector< vector<int> >& domains, const vector<int> currentDomain, vector<int>& subsample, const int& n, const int& k, int var)
{
  if (int(subsample.size()) == k) return 1;
  if (var >= n) return 0;
  int newVar = selectNewVar(currentDomain);
  bool result;
  if (newVar != -1)
    {
      vector<int> varDomain = domains[newVar];
      vector<int> newDomain;

      //cout << "\t\t new variable selected:\t" << newVar << endl;

      subsample.push_back(newVar);
      
      intersect(currentDomain, varDomain, newDomain);
      /*
      cout << "\t\t\t current domain:\t ";
      printDomain(currentDomain);
      cout << "\n\t\t\t variable domain:\t ";
      printDomain(varDomain);
      cout << "\n\t\t\t new domain:\t\t ";
      printDomain(newDomain);
      cout << endl;
      */
      if (validDomain(newDomain, subsample, k))
	{
	  //cout << "\t\t call to backtracking search\n";
	  result = backtrackingSearch(domains, newDomain, subsample, n, k, newVar);
	  if (result) return 1;
	}
    }
  subsample.pop_back();
  //var = subsample.back();
  return 0;
}


bool randomDraw(const vector< vector<int> >& domains, vector<int> currentDomain, vector<int>& subsample, const int& n, const int& k)
{
  bool success = 1;
    
  while ((int(subsample.size())<k) && success)
    {

      int var = selectRandomVar(currentDomain);
      if (var == -1) return 0;
      subsample.push_back(var);
      vector<int> varDomain = domains[var];
      vector<int> newDomain;
      intersect(currentDomain, varDomain, newDomain);
      currentDomain = newDomain;
      if (!validDomain(currentDomain,subsample,k)) success=0;
    }
  return success;
}


int selectNewVar(const vector<int>& currentDomain)
{
  if (currentDomain.empty()) return -1;
  else return currentDomain[0];
}


int selectRandomVar(const vector<int>& currentDomain)
{
  if (currentDomain.empty()) return -1;
  else return currentDomain[rand()%currentDomain.size()];
}


void generateFirstSamples(list< vector<int> >& subsamples, vector< vector<int> >& domains, const int& n, const int& k)
{
  vector<int> subsample(k);
  for (int i=0; i<n; ++i)
    {
      subsample[i%k] = i;
      if ((i%k) == (k-1))
	{
	  subsamples.push_back(subsample);
	  updateDomains(domains, subsample);
	}
    }
}


bool validDomain(const vector<int>& currentDomain, const vector<int>& subsample, const int& k)
{
  if (currentDomain.size() < (k-subsample.size())) return 0;
  else return 1;
}


void updateDomains(vector< vector<int> >& domains, const vector<int>& subsample)
{
  for (vector<int>::const_iterator it = subsample.begin(); it != subsample.end(); ++it)
    {
      for (vector<int>::const_iterator jt = it+1; jt != subsample.end(); ++jt)
	{
	  removeValuePair(domains, *it, *jt);
	}
    }
}
  
void removeValuePair(vector< vector<int> >& domains, const int& val1, const int& val2)
{
  for (vector<int>::iterator it = domains[val1].begin(); it != domains[val1].end(); ++it)
    {
      if (*it == val2)
	{
	  domains[val1].erase(it);
	  break;
	}
    }
  for (vector<int>::iterator it = domains[val2].begin(); it != domains[val2].end(); ++it)
    {
      if (*it == val1)
	{
	  domains[val2].erase(it);
	  break;
	}
    }
}


void initializeDomains(vector< vector<int> >& domains, const int& n)
{
  vector<int> temp(n);
  int i=0;
  for (vector<int>::iterator it = temp.begin(); it != temp.end(); ++it)
    {
      *it = i;
      ++i;
    }

  domains.resize(n,temp);

  for (int i=0; i<n; ++i)
    {
      domains[i].erase(domains[i].begin()+i);
    }
}


void intersect(const vector<int>& domain1, const vector<int>& domain2, vector<int>& intersection)
{
  for (vector<int>::const_iterator it = domain1.begin(); it != domain1.end(); ++it)
    {
      for (vector<int>::const_iterator jt = domain2.begin(); jt != domain2.end(); ++jt)
	{
	  if (*it == *jt) intersection.push_back(*it);
	}
    }
}
