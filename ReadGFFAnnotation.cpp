#include "ReadGFFAnnotation.h"

using namespace std;

/*------------------------------------------------------------------------------- 
This file processes a FASTA formatted file of a series of sequences and returns
e.g. a map or vectors
--------------------------------------------------------------------------------*/

//Constructor:  assumes a lower case sequence including a, c, g, t, u only
ReadGFFAnnotation::ReadGFFAnnotation(string fileName, map < string, string >& chrMap, bool u)
{
    cout<<"Reading Annotation file "<<fileName<<endl;
		
	uniqueOnly=u;
	
	readFile(fileName, chrMap);
	cout<<"Annotation file read"<<endl;
	
}
/***********************************************************************************
	Getters
***********************************************************************************/
const vector <AnnotationRecord> & ReadGFFAnnotation::getAnnotations()
{
	return annotations;
}

const vector <AnnotationRecord> & ReadGFFAnnotation::getGenes()
{
	return genes;
}

const vector <GeneRelationships> & ReadGFFAnnotation::getGeneRelationships()
{
	unsigned int currPosition=0;
	bool found=0;

	for(unsigned int annCtr=0; annCtr<annotations.size(); ++annCtr)
	{
		found=0;

		for(unsigned int geneCtr=0; geneCtr<genes.size(); ++geneCtr)
		{
			if(found==1 && annotations[annCtr].chromosome!=genes[geneCtr].chromosome)
			{
				geneCtr=genes.size();
			}

			else if((annotations[annCtr].startPos >=genes[geneCtr].startPos && annotations[annCtr].startPos <=genes[geneCtr].endPos))
			{
				GeneRelationships gr;
				gr.annotation_id=annotations[annCtr].id;
				gr.gene_id=genes[geneCtr].id;
				geneRelationships.push_back(gr);
			}		
		}	
	}

	return geneRelationships;
}
/***********************************************************************************
	Read File
***********************************************************************************/

//Reads in the Annotation file
void ReadGFFAnnotation::readFile(string fileName, map < string, string >& chrMap)
{
	Handy h(0);

	ifstream inAnnotation(fileName.c_str(), ios::in);	 
	string line;
	

    // ensure that the file can be opened
	if ( !inAnnotation )
	{
		cerr << "ERROR!! Annotation file: " << fileName << " could not be opened" << endl;
		exit (1);
	}
  
	cout<<"Annotation file opened, about to begin reading in."<<endl;

	int lineCt=0;

	bool fastaSeqs=0;
	while (getline(inAnnotation, line , '\n' ) && fastaSeqs==0) 
	{	
	

		if(lineCt%100000==0)
		{
			cout<<"ReadGFFAnnotation: Reading line "<<lineCt<<endl;	
				
		}

		//Some gff annotations put fasta sequences at end.  We want to ignore them
		if(line.length()>0 && line.substr(0,1)==">")
		{
			fastaSeqs=1;		
		}

		if (line.length()>0 && line.substr(0,1)!="#" && fastaSeqs==0)
		{	
			vector <string> record=h.getSplitString(line, "\t");			
			
			

			if(record.size()>=9 )
			{
						
				string uniqueInfo=record[0]+record[3]+record[4]+record[6];
				
				if(uniqueOnly==0 || uniqueAnns.find(uniqueInfo)==uniqueAnns.end())//If it wants the non unique or if this one hasn't been found yet
				{
					if(uniqueOnly==1)
					{
						uniqueAnns[uniqueInfo]=1;
					}
					
					AnnotationRecord currAnn;

					//This expects chromosome in gff file to be e.g. chr5
					
					currAnn.id=h.getStringFromInt(lineCt);
				
					if(chrMap.find(record[0])==chrMap.end())
					{
						cout<<"WARNING: Chromosome "<< record[0] <<" not found in chromosome map file." <<endl;
						currAnn.chromosome=record[0];
					}
					else
					{
						currAnn.chromosome=chrMap[record[0]];
					}
					

					if(record[6].length()>0)
					{
						currAnn.strand=record[6][0];
					}

					currAnn.startPos=h.getIntFromString(record[3])-1; //assumes annotation will start at 1 but alignment starts at 0
					currAnn.endPos=h.getIntFromString(record[4])-1; //same
					currAnn.annotationType=record[2];

					//Cycling through the last, long gff record to get the record name
					vector <string> lastRecord=h.getSplitString(record[8], ";");	
					currAnn.annName=record[8];	

						
					if(currAnn.annotationType=="gene")
					{
						genes.push_back(currAnn);				
					}

					else 
					{
						annotations.push_back(currAnn);				
					}
					

				}
			}//close cds
		}	

		++lineCt;
	}//close while loop


	map <string, int> names;

	for(unsigned int i=0; i<annotations.size(); ++i)
	{
		++names[annotations[i].id];
	}

}

