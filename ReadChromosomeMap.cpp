#include "ReadChromosomeMap.h"
using namespace std;

/*------------------------------------------------------------------------------- 
This file processes a file containing the chromosome map file which is formatted 
as follows:

Chromosome name from the annotation \t Chromosome name from the alignment
--------------------------------------------------------------------------------*/

//Constructor:  assumes a lower case sequence including a, c, g, t, u only
ReadChromosomeMap::ReadChromosomeMap(string inFileName, map<string, string>& chrMap)
{
	Handy h(0);


	ifstream inFile(inFileName.c_str(), ios::in);	 
	string line;
	

    // ensure that the file can be opened
	if ( !inFile )
	{
		cerr << "ERROR!! Chromosome Map: " << inFileName << " could not be opened" << endl;
		exit (1);
	}
  
	int linect=0;


	while (getline(inFile, line , '\n' )) 
	{
		vector <string> record=h.getSplitString(line, "\t");

		if(record.size()>=2)
		{
			chrMap.insert(map< string, string >::value_type(record[0], record[1]));
		}
		else
		{
			cout<<"Line "<<linect<<" of "<<inFileName<<" does not have the right number of records.  Line: "
				<<line<<endl;//maybe here
		}

		++linect;
	}
}

