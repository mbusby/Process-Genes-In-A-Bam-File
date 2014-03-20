Process-Genes-In-A-Bam-File
===========================

This is the custom script that was used to quantify the strand-specific reads mapped to genes so we could caluclate 
the RPKM in Smanski et al. 2014.  While we mainly uploaded it to ensure our data is reproducible it is a general purpose
script that uses Derek Barnett's Bamtools API to count many things that we find useful in the lab.


What does it do?
Produces a tab-delimitted text file with various attributes quantified for each gene that is in an input GTF file. 

What would I want to do that?
There are lots of different stuff you might need for different purposes.


Installation:

The binary is compiled for Unix.  That is the easiest way to go.

Compiling the data requires the bamtools API.  You will also need to edit your Make file to point to your install of 
bamtools.


--------------------------------------
How does it work? 

Type:

./CountByGenesStrandedness -bam alignmentName.bam -bamIndex alignmentName.bam.bai -gtf annotation.gtf -chrMap chrMap.txt -out outputFile.txt 


-bam - the bam alignment
-bamIndex - an index of the file 
-gtf - the annotation file
-chrMap - sometimes chromsome names are different in the alignment and the GTF.  
-out - the output file


Bam Index
The bam file needs an index.  Usually it has the same name as the alignment file with the addition of .bai.  If you cannot find one, 
it is easy to make using samtools.

From a unix command line type:

samtools index alignmentName.bam


GTF

The GTF is a an annotation file in a standard GTF format.  It is a tab-delimited file (make it in Excel and save it as Text - tab delimmited) using the format here:
http://www.ensembl.org/info/website/upload/gff.html

The third column (feature - feature type name, e.g. Gene, Variation, Similarity) should say "gene" or the program won't recognize it.  It ignores everything that is not a gene.

chrMap
This is a map that maps the chromosome name as it is defined in the bam file with the chromosome name as it is defined in the GTF.

This is required whether the names in the two files are the same or not.

The file is a tab-delimitted text file as follows:

Chr_name_from_the_annotation[tab]Chr_name_from_the_alignment

e.g.:

chr1	1
mito	MT
chrY	Y

To see the names of the chromosomes that are used in the .bam file type:

samtools view -H alignmentName.bam

For example, chromosome 1 is found after the sequence name (SN:)
@SQ     SN:1 


---------------------------------------------
What do the output fields mean?

AnnotationName - name of the annotation
AnnotationType - what type of annotation it is (always gene)
Chromosome - chromosome annotation is on
StartPosition - start position of annotation
EndPos - end position of annotation (add one for length!)
totalReads- how many reads align to the gene.  This includes all reads that overlap the gene at all, even by one base.  Each read pair counts as two reads.  Multiply mapped reads are counted whereever they fall.
uniquePairs - how many pairs align uniquely to the read.
duplicateUniquePairs - how many of the pairs that align uniquely to the gene are duplicates.
tandemUniquePairs - how many of the unique pairs are tandem reads
uniquePairEitherSide - how many unique reads where either read 1 or read 2 is inside the annotation region.
uniqueReads - how many read align uniquely to this region.
totalReadsSameStrand - how many reads are on the same strand as the annotation.  I define "same strand" as follows:
	First read on the same strand as the annotation
	 (e.g. + strand if the annotation strand is +)
	Second read on the opposite strand at the annotation  
	(e.g. - strand if the annotation strand is +)
uniqueReadsSameStrand - how many uniquely aligned reads are on the same strand as the annotation
uniquePairEitherSideSameStrand-how many unique reads where either read 1 or read 2 is inside the annotation region and on the same side as the annotation
UniquePairsBothIn-the number of unique pairs where both reads are completely inside the annotation region.



