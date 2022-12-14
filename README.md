# Intron Classification Pipeline

Scripts for the classification of introns into major, major-like, hybrid, minor-like and minor bins. The pipeline has been described in the following preprint:

***Taxonomy of introns, their evolution, and the role of minor introns in stress response***

*Anouk M. Olthof<sup>1,2</sup>, Charles F. Schwoerer<sup>1</sup>, Audrey L. Weber<sup>1</sup>, Iswarya Arokiadhas<sup>3</sup>, Karen Doggett<sup>4</sup>, Stephen Mieruszynski<sup>4</sup>, Avner Cnaani<sup>3</sup>, Joan K. Heath<sup>4</sup>, Jakob Biran<sup>3</sup> & Rahul N. Kanadia<sup>1,5</sup>*

*<sup>1</sup>Physiology and Neurobiology Department, University of Connecticut, Storrs, CT, 06269, USA*<br>*<sup>2</sup>Current address: Institutite for Cellular and Molecular Medicine, University of Copenhagen, Copenhagen, Denmark*<br>*<sup>3</sup>Department of Poultry and Aquaculture, Institute of Animal Science, Agricultural Research Organization, Rishon LeTsiyon, Israel*<br>*<sup>4</sup>Walter and Eliza Hall Institute of Medical Research, Parkville, VIC 3052, Australia*<br>*<sup>5</sup>Institute of Systems Genomics, University of Connecticut, Storrs, CT, 06269, USA*
___

[splice-site-generation.md](splice-site-generation.md)

This is a bash markdown document for the intron classification pipeline that calls all other required scripts.

[extract-U2BPS.md](extract-U2BPS.md)

This is a bash markdown document for the extraction of splice site sequences from intron FASTA files.

[getNucleotideCounts_1.md](getNucleotideCounts_1.md)

This is a bash markdown document for the initial binning of introns by terminal dinucleotide sequences.

[mergePWM_1.md](mergePWM_1.md)

This is a bash markdown document for the generation of initial position-weight matrices.

[scoreInitialPWM.md](scoreInitialPWM.md)

This is a bash markdown document for the scoring of introns against the initial position-weight matrices.

[getNucleotideCounts_2.md](getNucleotideCounts_2.md)

This is a bash markdown document for the binning of introns based on the initial PWM scores.

[mergePWM_2.md](mergePWM_2.md)

This is a bash markdown document for the generation of refined position-weight matrices.

[scoreRefinedPWM.md](scoreRefinedPWM.md)

This is a bash markdown document for the scoring of introns against the refined position-weight matrices.

[binIntrons_afterScoring.md](binIntrons_afterScoring.md)

This is a bash markdown document for the binning of introns based on the refined PWM scores.

[scripts](scripts)

This directory contains plain script versions of the markdown documents.	
