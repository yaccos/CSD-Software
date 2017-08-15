#include <string>
#include <vector>

#include "myStructs.h"

using namespace std;


// FUNCTION CALLED FOR PREPROCESSING
void filterAndSortGenes(vector< vector<string> >& allIDs, vector< vector<Gene> >& allGenes, vector< vector<Gene> >& filteredGenes, const int& threadCount);
// The function sorts all genes by their ID and then add the genes that occur in all the 
// input-files to the "filteredGenes"-container. Each "Gene"-structure  contains a gene ID,
// the index at which it occurs in its respective file, and an index indicating the file in
// mention. These indices are used in the file-reading step that follows and allows for all the
// matrices containing the measurements to store the gene-measurements in the same order.
//
// Input parameters:
//           - allIDs: a vector of vectors containing the gene-IDs from each of the input
//                     files.
//           - allGenes: an empty vector of vectors for storing Gene-stucts, which serves
//                       as a temporary storage for the Gene-structs during sorting and
//                       filtering.
//           - filteredGenes: an empty vector of vectors for storing Gene-structs, in which
//                            the sorted and filtered Gene-structs are returned.



// UTILITY FUNCTIONS
void initializeGeneVectors(vector< vector<string> >& allIDs, vector< vector<Gene> >& allGenes, const int& threadCount);
// The function reads the gene-IDs from allIDs, converts them to Gene-structs and return
// them in allGenes.


bool indicesInRange(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices);
// The function iterates through the "indices"-vector checking that all indices are in
// range of their related Gene-vector in "allGenes".


bool allAreEqual(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices);
// The function checks whether the gene-ID of all Gene-structs pointed to by an index
// in "indices" are idetical.


void incrementAllIndices(vector<unsigned int>& indices);
// The function increments all indices in "indices" by +1.


void findBiggest(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices, Gene& biggest);
// The function compares the gene-IDs pointed to by the indices in "indices" and returns the
// the "biggest" by comparison using the funtion geneCompare().


void pushToFilteredGenes(vector< vector<Gene> >& allGenes, const vector<unsigned int>& indices, vector< vector<Gene> >& filteredGenes);
// The function push all Gene-structs currently pointed to by the indices in "indices" to their
// respective vectors in "filteredGenes".

