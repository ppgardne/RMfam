# RMfam: multiple alignments and probabilistic models of RNA motifs. 

### Release 0.3

0. CONTENTS
1. INTRODUCTION
2. AVAILABILITY
3. FILES
4. HISTORY
5. HOW TO CITE RMFAM

## 1. INTRODUCTION

RMfam is a database of annotated multiple sequence alignments and
probabilistic models for a number of RNA motifs. For the purposes of
this work an RNA motif is a recurring RNA sequence and/or secondary
structure found within larger structures that can be modelled by
either a covariance model or a profile HMM. RMfam is modelled loosely
on the Rfam and Pfam databases and draws a lot of ideas from these
resources.  

RMfam is licensed under a Creative Commons Attribution 3.0 Unported License
(http://creativecommons.org/licenses/by/3.0/).

## 2. AVAILABILITY

RMfam is available on the web at the following URL:

`https://github.com/ppgardne`

## 3. FILES

<pre>
README                 - this file
USERMAN                - a description of the RMfam flatfile formats
LICENSE                - Legalese stating RMfam's copyright is CC BY 3.0
*/SEED                 - Annotated alignment files of motifs (named *)
RMfam.cm               - a concatenated set of RMfam covariance models in ascii INFERNAL format
Rfam.stk               - Rfam seed alignments in STOCKHOLM format that have been annotated with RMfam motifs
Rfam.tab               - A table of mappings between RMfam and Rfam. The format is:
		          --column 1: Rfam accessions                
		          --column 2: Rfam ID                
		          --column 3: RMfam ID  
		          --column 4: fraction of sequences in Rfam seed that contain the RMfam motif          
		          --column 5: sum of bit scores
			  --column 6: weighted sum of bit scores (using GSC tree weights from esl-weights, 
			              weight the sum of bits, normalise with the total weights). 
				      I.e. (\sum_i (w_i * bits_i))/( \sum_i w_i )
PDB.gff                - Annotations of nucleotide sequences derived from PDB entries
accessions.txt         - Tab-delimited file of RMfam accessions and corresponding IDs
rmfam_scan.pl          - A script for annotating sequences and alignments with RMfam families
</pre>

## 4. HISTORY

<pre>
Version     Date        # families    Rfam        PDB
-------     ----        ----------    -------     -----
  0.1       10/11       14            Rfam 10.1   03/11/11
  0.2       05/12       24            Rfam 10.1   03/11/11
  0.3	    06/14	34	      Rfam 11.0   01/06/14
</pre>

## 5. HOW TO CITE RMfam

If you make use of RMfam in your work we ask that you cite the
following publications:

`Gardner PP and Eldai H (2015) Annotating RNA motifs in sequences and alignments. 
Nucl. Acids Res. 43 (2): 691-698.`

Feb 2016
