#include <string>

#include "myStructs.h"

using namespace std;




void setInputType(const string& parameter, InputType& input);
// The function validates and sets the input parameter. 


void setOutput(const string& parameter, Output& o);
// The function validates and sets the output parameter.


void setCorrMethod(const string& parameter, string& c);
// The function validates and sets the correlation method.


void setFormat(const string& parameter, bool& vertical);
// The function validates the parameter and sets the the data format to
// vertical or not (not=horizontal).


void validateJob(const InputType& input, const Output& o);
// The function checks for illegal combinations of input and output parameters.
