# Intron Classification Pipeline

Scripts for the classification of introns into major, major-like, hybrid, minor-like and minor bins. The pipeline has been described in the following preprint:

***Taxonomy of introns, their evolution, and the role of minor introns in stress response***

*Anouk M. Olthof<sup>1,2</sup>, Charles F. Schwoerer<sup>1</sup>, Audrey L. Weber<sup>1</sup>, Iswarya Arokiadhas<sup>3</sup>, Karen Doggett<sup>4</sup>, Stephen Mieruszynski<sup>4</sup>, Avner Cnaani<sup>3</sup>, Joan K. Heath<sup>4</sup>, Jakob Biran<sup>3</sup> & Rahul N. Kanadia<sup>1,5</sup>*

*<sup>1</sup>Physiology and Neurobiology Department, University of Connecticut, Storrs, CT, 06269, USA*<br>*<sup>2</sup>Current address: Institutite for Cellular and Molecular Medicine, University of Copenhagen, Copenhagen, Denmark*<br>*<sup>3</sup>Department of Poultry and Aquaculture, Institute of Animal Science, Agricultural Research Organization, Rishon LeTsiyon, Israel*<br>*<sup>4</sup>Walter and Eliza Hall Institute of Medical Research, Parkville, VIC 3052, Australia*<br>*<sup>5</sup>Institute of Systems Genomics, University of Connecticut, Storrs, CT, 06269, USA*
___

*INTRODUCTION*

The Intron Classification Pipeline is a suite of scripts for classifying introns into major, major-like, hybrid, minor-like, and minor bins. This pipeline is designed to facilitate the study of intron taxonomy, their evolution, and the role of minor introns in stress response, as detailed in Olthof and Schwoerer 2023. 

*DEPENDENCIES:*

1. GNU's awk
2. bedtools (to generate IntronFASTAs)
3. python3 (to parse genome annotation)

*INSTALLATION:*

`git clone <github URL>`

*USAGE:*

`bash bin/RunClassificationPipeline.sh`

*WORKFLOW:*

1. Extract intron FASTA sequences from all species in MIDB_v2.0-SpeciesList.tsv (generates IntronFASTAs)
2. Generate initial PWMs using introns from 263 species (generates PWMs/Initial)
3. Scoring: Score introns using the initial PWMs (generates output/${species}/Initial)
4. Generate 169 titrated PWMs using introns from 263 species (generates PWMs/$(t})
5. Calculate precision/recall based using responsive introns in PipelineCalibration
6. Select the optimum titration by maximizing the average F1 score across species
7. Classify introns based on the selected titration (generates output/${species}/${t})

Please note that each time `RunClassificationPipeline.sh` is called, it will go through each step, check to see if the output directory exists, and responds accordingly. Thus, to bypass computationally heavy steps (such as FASTA or PWM generation) one can extract the corresponding files from FigShare: https://figshare.com/projects/Taxonomy_of_introns_and_the_evolution_of_minor_introns/167411. In order to classify an additional species, add that species to MIDB_v2.0-SpeciesList.tsv.

*CONTACT*

For queries or assistance with the pipeline, please contact:

Anouk M. Olthof: anouk.olthof@uconn.edu
Charles F. Schwoerer: charles.schwoerer@columbia.edu
