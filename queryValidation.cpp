#include <string>
#include <cstdlib>
#include <iostream>

#include "queryValidation.h"
#include "myStructs.h"

using namespace std;


void setInputType(const string& parameter, InputType& input)
{
  if (!parameter.compare("ged"))
    {
      input.raw = 1;
      input.corr = 0;
      input.unfiltered = 0;
    }
  else if (!parameter.compare("corr"))
    {
      input.raw = 0;
      input.corr = 1;
      input.unfiltered = 0;
    }
  else if (!parameter.compare("csd"))
    {
      input.raw = 0;
      input.corr = 0;
      input.unfiltered = 1;
    }
  else 
    {
      cout << "Error: Unknown argument for input.\n";
      exit(1);
    }
}

void setOutput(const string& parameter, Output& o)
{
  if (!parameter.compare("corr"))
    {
      o.corr = 1;
    }
  else if (!parameter.compare("csd"))
    {
     o.unfiltered = 1;
    }
  else if (!parameter.compare("network"))
    {
      o.filtered = 1;
    }
  else if (!parameter.compare("all"))
    {
      o.all = 1;
    }
  else
    {
      cout << "Error: Unknown argument for output.\n";
      exit(1);
    }
}


void setCorrMethod(const string& parameter, string& c)
{
  if (!parameter.compare("pearson"))
    {
      c = parameter;
    }
  else if (!parameter.compare("spearman"))
    {
      c = parameter;
    }
  else
    {
      cout << "Error: Unknown argument for correlation method.\n";
      exit(1);
    }
}


void setFormat(const string& parameter, bool& vertical)
{
  if (!parameter.compare("horizontal")) vertical = 0;
  else if (!parameter.compare("vertical")) vertical = 1;
  else
    {
      cout << "Error: Unknown argument for data format given. Expected either\n"
	   << "\"horizontal\" or \"vertical\"\n";
      exit(1);
    }
}


void validateJob(const InputType& input, const Output& o)
{
  if ((input.corr && o.corr) || (input.unfiltered && o.unfiltered))
    {
      cout << "Error: Input type = Output type.\n";
      exit(1);
    }
  if (input.unfiltered && o.corr)
    {
      cout << "Error: Cannot calculate correlation coefficients with unfiltered\n"
	   << "CSD-scores as input.";
      exit(1);
    }
}

