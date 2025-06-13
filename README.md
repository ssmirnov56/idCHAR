# idCHAR
version D. Processing or separtion of prolines added.

version C. Two reporting bugs fixed: for residue position and for single-digit % values.

version B. Percent-score is added for the output: score/#aa

version A. Processes partitioning by charge (KR vs ED) and by aromatics (FYW)


Purpose:

Finds the point in aa sequence which conveys the greatest partitioning (difference) in target feature of aa composition.

Usage:

 ./idCHAR_a.pl [-h] [-v] [-d] -fasta {fasta_filename} [-charge] [-aroma] [-pro]

Params:

	-fasta:		Protein sequence fasta filename;

  -charge:	Sequence partitioning by charge: Lys+Arg vs Asp+Glu content;
	
  -aroma:	Sequence partitioning by aromatics content: Trp+Tyr+Phe content vs. their absence;
	
  -pro:	Sequence partitioning by proline content: Pro vs its absence;

(C):

 Serge Smirnov, Western Washington University, smirnos@wwu.edu, 2020
