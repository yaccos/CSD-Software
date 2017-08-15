#include <vector>
#include <list>

using namespace std;



// FUNCTIONS FOR GENERATING COLLECTIONS OF SUBSAMPLES
void generateSubsamples(list< vector<int> >& subsamples, const int& n, const int& k);


void generateSubsamplesNoInitialSampling(list< vector<int> >& subsamples, const int& n, const int& k);


void randomSubsampling(list< vector<int> >& subsamples, const int& n, const int& k, const int& terminationLimit);


void randomSubsampling(vector< vector<int> >& subsamples, const int& n, const int& k, const int& terminationLimit);


void randomSubsamplingWithLimit(vector< vector<int> >& subsamples, const int& n, const int& k, const int& failLimit, const int& maxNumSubsamples, const int& randomSeed);



// FUNCTIONS FOR DRAWING INDIVIDUAL SUBSAMPLES
bool backtrackingSearch(const vector< vector<int> >& domains, const vector<int> currentDomain, vector<int>& subsample, const int& n, const int& k, int var);


bool randomDraw(const vector< vector<int> >& domains, vector<int> currentDomain, vector<int>& subsample, const int& n, const int& k);




// FUNCTIONS FOR SELECTING VARIABLE FROM A DOMAIN 
int selectNewVar(const vector<int>& currentDomain);


int selectRandomVar(const vector<int>& currentDomain);




// UTILITY FUNCTIONS
void generateFirstSamples(list< vector<int> >& subsamples, vector< vector<int> >& domains, const int& n, const int& k);


bool validDomain(const vector<int>& currentDomain, const vector<int>& subsample, const int& k);


void updateDomains(vector< vector<int> >& domains, const vector<int>& subsample);


void removeValuePair(vector< vector<int> >& domains, const int& val1, const int& val2);


void initializeDomains(vector< vector<int> >& domains, const int& n);


void intersect(const vector<int>& domain1, const vector<int>& domain2, vector<int>& intersection);
