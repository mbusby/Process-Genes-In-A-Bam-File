/*------------------------------------------------------------------------------- 
This class reads the SGD GFF annotation format
--------------------------------------------------------------------------------*/
#pragma once
#include <string>#include <vector>#include <iostream>#include <fstream>#include <algorithm>#include "Handy.h"struct AnnotationRecord{	string id;//unique identifier	string annName;//not necessarily unique idenifier because this is how gff does it	string chromosome;	char strand; //strand is +, -, or . for things like ARS (replication origin) which occurs on both strands	unsigned int startPos;	unsigned int endPos;	unsigned int readsAligned;	unsigned int uniqueReads;	double normalizedReads;	unsigned int uniqueLength;	string annotationType;		AnnotationRecord(void)		: id("-")		, annName("-")		, chromosome("-")		, strand('.')		, startPos(0)		, endPos(0)		, readsAligned(0)		, uniqueReads(0)		, normalizedReads(0)		, uniqueLength(0)		, annotationType("-")	{}	//sorting function for multiple things in structure in vector		bool operator<(const AnnotationRecord& as) const 	{ 						if(chromosome == as.chromosome && startPos==as.startPos && endPos == as.endPos && annotationType == as.annotationType) return annName<as.annName;		if(chromosome == as.chromosome && startPos==as.startPos && endPos == as.endPos) return annotationType<as.annotationType;		if(chromosome == as.chromosome && startPos==as.startPos) return endPos < as.endPos; //chromosome start pos same then end		if(chromosome == as.chromosome) return startPos < as.startPos; // chromosome then start		return chromosome < as.chromosome;  //else just chromosome	}};struct GeneRelationships{	string gene_id;	string annotation_id;};using namespace std;
class ReadGFFAnnotation
{
	public:
		ReadGFFAnnotation(string , map < string, string >& , bool );
		const vector <AnnotationRecord> & getAnnotations();
		const vector <AnnotationRecord> & getGenes();
		const vector <GeneRelationships> & getGeneRelationships();		
	private:
		vector <AnnotationRecord> annotations;   
		vector <AnnotationRecord> genes;   
		vector <GeneRelationships> geneRelationships; //From Parent annotations in GFF file		map < string, int > uniqueAnns; //For keeping track of unique annotations		bool uniqueOnly;  //Include unique regions only
		void readFile(string, map < string, string >&);	
};

