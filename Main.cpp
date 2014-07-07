#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "api/BamReader.h"
#include "Handy.h"
#include "ReadGFFAnnotation.h"
#include "ReadChromosomeMap.h"


unsigned int checkErrors();
void displayHelp();
void somethingsGoneWrong(string);
void processBamFile();
void getUnique();

//Variables taken from the inputs
string bamFileName="";
string bamIndexFileName="";
string gtfFileName="";
string outFileName="";
string chrMapFileName="";
map<string, string> chrMap;
vector <AnnotationRecord> annotations;
int totalReads=0;
int uniqueReads=0;
int uniquePairs=0;
int duplicateReads=0;
int duplicateUniqueReads=0;
int duplicateUniquePairs=0;
int tandemUniquePairs=0;
int tandemPairs=0;
int alignedPairs=0;
int uniquePairEitherSide=0;
int uniqueReadsSameStrand=0;
int totalReadsSameStrand=0;
int uniquePairEitherSideSameStrand=0;
int uniquePairsAllIn=0;
int overlappingBasesSameStrand=0;
int overlappingBasesOppositeStrand=0;
		
map<string, doubleCt> ctReads;


map <string, tripleCt> unmappedReadCts;
double mapQualCut=0;

using namespace std;
using namespace BamTools;

int main(int argc, char* argv[]) 
{
		
	//Intialize library of useful things
	Handy h(0);

	bool booltest=true;
	
	int optind=1;


	while ((optind < argc) && (argv[optind][0]=='-')) 
	{	
        string sw = argv[optind];
		
		if (sw=="-h") 
		{	
            optind++;
			displayHelp();
			return 1;
        }
		
		if(optind >=  argc-1)
		{
			cerr<<"Your final parameter, "<<sw<<" is missing a value."<<endl;
			return 1;
		}

		if (sw=="-bam")
		{	
            optind++;
			bamFileName = argv[optind];		
			optind++;
        }
		
		else if (sw=="-bamIndex")
		{	
            optind++;
			bamIndexFileName = argv[optind];			
			optind++;
        }
		
		
		else if (sw=="-gtf")
		{	
            optind++;
			gtfFileName = argv[optind];			
			optind++;
        }
		
		else if (sw=="-chrMap")
		{	
            optind++;
			chrMapFileName = argv[optind];			
			optind++;
        }
		
		else if (sw=="-out")
		{	
            optind++;
			outFileName = argv[optind];		
			optind++;
        }
					
		else
		{
			cerr<<"Main: Unknown parameter:"<<sw<<endl;
			return 1;
		}
	}	
	
	unsigned int errCheck = checkErrors();

	if(errCheck==0)
	{
		
	
	}
	
	ReadChromosomeMap(chrMapFileName, chrMap);
	cout<<"Chromosome Map Read"<<endl;
	
	ReadGFFAnnotation ran(gtfFileName, chrMap, false);
	annotations=ran.getAnnotations();	
	
	cout<<"Annotations read. Size is "<<annotations.size()<<endl;
	
	getUnique();
		
	cout.flush();
	
	processBamFile();
		
}


/*==============================================================================
Process the bam file 

==============================================================================*/
void processBamFile()
{

	BamReader reader;
	BamAlignment al;
	string readSeq;
	tripleCt counts;
	Handy h(0);
	RefVector refVector;
	map< string, bool> missingChromosomes;
	
	remove(outFileName.c_str());
	ofstream outputStream;
	outputStream.open(outFileName.c_str());
			
	outputStream<<"AnnotationName\tAnnotationType\t";
	outputStream<<"Chromosome\t";
	outputStream<<"StartPosition\tEndPos\t";
	outputStream<<"totalReads\tuniquePairs\t";
	outputStream<<"duplicateUniquePairs\tandemUniquePairs\tuniquePairEitherSide\tuniqueReads\t";
	outputStream<<"totalReadsSameStrand\tuniqueReadsSameStrand\tuniquePairEitherSideSameStrand\tOverlappingBasesSameStrand\tOverlappingBasesOppositeStrand"<<endl;
			
	
	
	if ( !reader.Open(bamFileName) ) {
		cerr << "Could not open input BAM file." << endl;
		return;
	}
	
	//reader.CreateIndex();
	
	refVector=reader.GetReferenceData();
	
	for(int i=0; i<refVector.size(); ++i)
	{
		cout<<refVector[i].RefName<<'|'<<endl;
	}
	
	if ( !reader.OpenIndex(bamIndexFileName) ) {
		cerr << "Could not open input BAM index file."<<bamIndexFileName << endl;
		return;
	}
	
	if ( reader.HasIndex() ) {
		cout<<"Reader has indexes."<<endl;
	}
	else
	{
		cout<<"No Indexes found."<<endl;
	}	
		
	
	for(int i=0; i<annotations.size(); ++i)
	{	
	
		if(i%10000==0)
		{
			cout<<"Reading annotation line "<<i<<endl;
		}
		//reset values to zero
		totalReads=0;
		uniquePairs=0;
		duplicateUniquePairs=0;
		tandemUniquePairs=0;
		duplicateReads=0;		
		uniquePairEitherSide=0;
		uniqueReads=0;
		uniqueReadsSameStrand=0;
		totalReadsSameStrand=0;
		uniquePairEitherSideSameStrand=0;
		uniquePairsAllIn=0;
		overlappingBasesSameStrand=0;
		overlappingBasesOppositeStrand=0;
			
		AnnotationRecord thisAnn=annotations[i];
		
		int refID= reader.GetReferenceID(thisAnn.chromosome);
		
		BamRegion region(refID, (int) thisAnn.startPos, refID, (int) thisAnn.endPos );	
				
		if(refID==-1)
		{
			if(missingChromosomes.find(thisAnn.chromosome)==missingChromosomes.end())
			{
				missingChromosomes[thisAnn.chromosome]=1;
				cout<<"Warning! This chromosome not mapped in bam file: "<<thisAnn.chromosome<<endl;
			}
		}				
		else if ( !reader.SetRegion(region.LeftRefID, region.LeftPosition, region.RightRefID, region.RightPosition) )
		{
			
			//If the region is not found in the bam file write to output
			cout<<"Warning! Region not found in bam file: "<<thisAnn.annName<<" at "<<thisAnn.chromosome<<" "<<thisAnn.startPos<<"-"<<thisAnn.endPos<<endl;
		}
		else
		{
			reader.SetRegion(region.LeftRefID, region.LeftPosition, region.RightRefID, region.RightPosition) ;
			while ( reader.GetNextAlignment(al) )
			{				
				++totalReads;
				
				if(al.IsDuplicate()==true)
				{
					++duplicateReads;
				}
				
				int s=al.Position; //Start
				int e=al.MatePosition;//End
					
				if(s>e)
				{
					swap(s,e);
				}
				
				if( ( al.IsPaired()==false || (al.IsFirstMate()==true && al.IsMateMapped()==true )) 
					&& (al.MateRefID==al.RefID) &&  (ctReads[al.Name].ct1==1 || ctReads[al.Name].ct2==1) ) 
				{		
					++uniquePairs;	
					
					if(al.IsDuplicate()==true)
					{
						++duplicateUniquePairs;
					}				
					
					if(al.IsMateReverseStrand() == al.IsReverseStrand()) 
					{
						++tandemUniquePairs;				
					}		
					
					if( al.MatePosition >(int) thisAnn.startPos && al.MatePosition <(int) thisAnn.endPos)
					{
						++uniquePairsAllIn;
					}			
					
					
				}
				
				//If it's the first mate or if it's the secend mate but the first mate is not mapped
				if( (al.IsPaired()==false || al.IsFirstMate()==true || al.IsMateMapped()==false)  &&  
					(ctReads[al.Name].ct1==1 || ctReads[al.Name].ct2==1) ) 
				{		
					++uniquePairEitherSide;
								
				}
				//If it's the first mate or if it's the secend mate but the first mate is not mapped
				if( (al.IsFirstMate()==true & ctReads[al.Name].ct1==1)
					|| (al.IsFirstMate()==false  &&  ctReads[al.Name].ct2==1) ) 
				{		
					++uniqueReads;
								
				}
				//if it is the same strand as the annotation
				if( (al.IsFirstMate()==true & ( (al.IsReverseStrand()==false & thisAnn.strand=='+') || (al.IsReverseStrand()==true & thisAnn.strand=='-') ) )
					|| (al.IsFirstMate()==false & ( (al.IsReverseStrand()==true & thisAnn.strand=='+') || (al.IsReverseStrand()==false & thisAnn.strand=='-') ) )
					)
				{
					//If it's the first mate or if it's the secend mate but the first mate is not mapped
					++totalReadsSameStrand;
					
					int thisOverlappingBases=al.Length;
					
					if (al.Position<thisAnn.startPos)
					{
						thisOverlappingBases=thisOverlappingBases-(thisAnn.startPos - al.Position);
					}
					
					if (al.Position+al.Length>thisAnn.endPos)
					{
						thisOverlappingBases=thisOverlappingBases-( (al.Position+al.Length)-thisAnn.endPos);
					}
					
					overlappingBasesSameStrand=overlappingBasesSameStrand+thisOverlappingBases;
					
					if( (al.IsFirstMate()==true & ctReads[al.Name].ct1==1)
						|| (al.IsFirstMate()==false  &&  ctReads[al.Name].ct2==1) ) 
					{		
						++uniqueReadsSameStrand;
									
					}
					
					if( (al.IsPaired()==false || al.IsFirstMate()==true || al.IsMateMapped()==false)  &&  
					(ctReads[al.Name].ct1==1 || ctReads[al.Name].ct2==1) ) 
					{		
						++uniquePairEitherSideSameStrand;
									
					}

				}
				else
				{
					int thisOverlappingBases=al.Length;
					
					if (al.Position<thisAnn.startPos)
					{
						thisOverlappingBases=thisOverlappingBases-(thisAnn.startPos - al.Position);
					}
					
					if (al.Position+al.Length>thisAnn.endPos)
					{
						thisOverlappingBases=thisOverlappingBases-( (al.Position+al.Length)-thisAnn.endPos);
					}
					
					overlappingBasesOppositeStrand=overlappingBasesOppositeStrand+thisOverlappingBases;
				
				}
			}				
		}
			
		outputStream<<thisAnn.annName<<"\t"<<thisAnn.annotationType<<"\t";
		outputStream<<thisAnn.chromosome<<"\t"<<thisAnn.startPos<<"\t"<<thisAnn.endPos<<"\t";
		outputStream<<totalReads<<"\t"<<uniquePairs<<"\t";
		outputStream<<duplicateUniquePairs<<"\t"<<tandemUniquePairs<<"\t"<<uniquePairEitherSide<<"\t"<<uniqueReads<<"\t";
		outputStream<<totalReadsSameStrand<<"\t"<<uniqueReadsSameStrand<<"\t"<<uniquePairEitherSideSameStrand<<"\t"<<overlappingBasesSameStrand<<"\t"<<overlappingBasesOppositeStrand<<endl;
		
	
			
	}		

	reader.Close();
	outputStream.close();
	
}
/*===========================================================================================
Figure out which reads are uniquely aligned
This step could be avoided if aligners used the same output to identify uniquely
aligned reads, but the coding is aligner-specific
===========================================================================================*/
void getUnique()
{
	Handy h(0);
	BamReader reader;
	BamAlignment al;
	RefVector refVector;
	int forwardAligned=0;
	int reverseAligned=0;

	
	if ( !reader.Open(bamFileName) ) 
	{
		cerr << "Could not open input BAM file." << endl;
		return;
	}
	

	int i=0;
		
	while(reader.GetNextAlignment(al))
	{	
		++i;		
		if(i%1000000==0)
		{
			cout<<"Reading line "<<i<<" of alignment for unique alignment map."<<endl;
		}
		
		if(ctReads.find(al.Name)==ctReads.end())
		{
			{
				doubleCt dc;
				dc.ct1=0;
				dc.ct2=0;
				ctReads.insert(map< string, doubleCt >::value_type(al.Name, dc));										
			}
		}				
	
	
		if(al.IsFirstMate()==true || al.IsPaired()==false)
		{						
			++ctReads[al.Name].ct1;
		}				
		else
		{
			++ctReads[al.Name].ct2;				
		}		

		if(al.IsMapped()==true)
		{
			if( (al.IsFirstMate()==true && al.IsReverseStrand()==false) || (al.IsFirstMate()==false && al.IsReverseStrand()==true) )
			{
				++forwardAligned;
			}
			else
			{
				++reverseAligned;
			}
		}
	}
	
	reader.Close();
	
	string outFileNameCt=outFileName+".cts";
	remove(outFileNameCt.c_str());
	ofstream outputStream;
	outputStream.open(outFileNameCt.c_str());		
	
	outputStream<<"Forward aligned: "<<forwardAligned<<endl;
	outputStream<<"Reverse aligned: "<<reverseAligned<<endl;
	
	outputStream.close();
}


/*===========================================================================================
Check that all of the necessary fields exist.
===========================================================================================*/

unsigned int checkErrors()
{
	//Errors 

	int err=0;
	Handy h(0);
	string problems="";
		
	if(bamFileName.length()==0)
	{
		problems.append("A bam file containing the reads is needed (-bam).\n");
		++err;
	}	
	
	if(bamIndexFileName.length()==0)
	{
		problems.append("An index file for the bam file is required(-bamIndex).\n");
		++err;
	}
	
	if(outFileName.length()==0)
	{
		problems.append("An file name is needed to write the output (-out).\n");	
		++err;
	}	
	
	if(chrMapFileName.length()==0)
	{
		problems.append("A chromosome map (-chrMap).\n");
		++err;
	}	
	
	/*============================================================
	Check that all the relevant files can be read/written to
	==============================================================*/

	err=err+h.checkRead(bamFileName);
	
	if (problems.length()>0)
	{
		somethingsGoneWrong(problems);
	}

	return err;
	
}


void displayHelp()
{
	
	
	cout<<"\nFor example, type something like this: \n";
	cout<<"./CountByGenesStrandedness -bam alignmentName.bam -bamIndex alignmentName.bam.bai -gtf annotation.gtf -chrMap chrMap.txt -out outputFile.txt"<<endl;
	cout<<"\n\nRequired options:\n";	
	cout<<"-bam Name of the bam file\n";
	cout<<"-bamIndex Name of the bam file index\n";
	cout<<"-gtf Annotation of genes\n";
	cout<<"-chrMap This is a map that maps the chromosome name as it is defined in the bam file with the chromosome name as it is defined in the GTF.\n";
	cout<<"-out Output file\n";
}


void somethingsGoneWrong(string whatsGoneWrong)
{

	cout<<"ERROR: Oh no! Something has gone horribly wrong.\n";	
	cout<<whatsGoneWrong;
	cout<<"\n";	
	cerr<<"\nPlease try again or display help using the -h option (CountByBarcode -h).\n";

	exit (EXIT_FAILURE);
}
