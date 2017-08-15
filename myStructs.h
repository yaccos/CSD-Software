#ifndef MYSTRUCTS_H
#define MYSTRUCTS_H

#include <string>
#include <vector>

using namespace std;



struct InputType
{
  bool raw;
  bool corr;
  bool unfiltered;

  InputType()
  {
    raw = 1;
    corr = 0;
    unfiltered = 0;
  }
};


struct Output
{
  bool corr;
  bool unfiltered;
  bool filtered;
  bool all;
   
  Output()
  {
    corr = 0;
    unfiltered = 0;
    filtered = 0;
    all = 0;
  }
};


struct Gene
{
  string geneID;
  int index;
  int sortedPos; 
  int origin;

  /*  
  Gene& operator=(const Gene& g)
  {
    //Gene a;
    geneID = g.geneID;
    index = g.index;
    origin = g.origin;
    return *this;
  }
  */
/* 
  void swap(Gene& g)
  {
    geneID = g.geneID;
    index = g.index;
    origin = g.origin;
  }
  */
};


struct CSD_Pair
{
  int gene_index1;
  int gene_index2;
  bool activated;
  vector<char> type;

  CSD_Pair()
  {
    gene_index1 = -1;
    gene_index2 = -1;
    activated = 0;
    type.resize(1);
  }
  CSD_Pair(unsigned int n)
  {
    gene_index1 = -1;
    gene_index2 = -1;
    activated = 0;
    type.resize(n);
  } 
};

#endif

