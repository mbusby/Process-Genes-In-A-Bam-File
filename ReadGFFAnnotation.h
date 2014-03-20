/*------------------------------------------------------------------------------- 
This class reads the SGD GFF annotation format
--------------------------------------------------------------------------------*/
#pragma once
#include <string>
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
		vector <GeneRelationships> geneRelationships; //From Parent annotations in GFF file
		void readFile(string, map < string, string >&);	
};
