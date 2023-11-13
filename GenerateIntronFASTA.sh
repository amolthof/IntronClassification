#!/bin/bash

# Set inputs: working directory

DIR=$1		# Working directory

# Set species name, assembly, source

species=$2
assembly=$3
source=$4

# Set source files

FileType=$5
AssemblyFile=${DIR}/output/${species}/${species}.${assembly}.${source}-genome.fa
AnnotationFile=${DIR}/output/${species}/${species}.${assembly}.${source}-genome.${FileType}

# Parse annotation file

if [ ${FileType} == "gtf" ]; then

	${DIR}/bin/parseGTF.sh ${DIR} ${species} ${assembly} ${source} ${AnnotationFile}

elif [ ${FileType} == "gff3" ]; then

	${DIR}/bin/parseGFF3.sh ${DIR} ${species} ${assembly} ${source} ${AnnotationFile}

else

	echo "ERROR: ${FileType} annotation unavailable"
	exit

fi

# Check if parsing worked

if [ ! -s ${DIR}/output/${species}/exon.info ]; then

	echo "ERROR: no exon.info for ${species}"
	exit

fi

# Extract unique transcript IDs

awk -F"\t" 'NR==1 {for (i=1; i<=NF; i++) if ($i=="transcript_id") col=i} NR>1 && col {print $col}' ${DIR}/output/${species}/exon.info | sort | uniq > ${DIR}/output/${species}/transcript_IDs.tmp

# Process transcript IDs

cat ${DIR}/output/${species}/transcript_IDs.tmp | parallel ${DIR}/bin/processTranscript.sh ${DIR} ${species} ${assembly} ${source} {}

# Generate final bed file 

find "${DIR}/output/${species}" -type f -name "*::exons.bed" -exec cat {} \; | sortBed > ${DIR}/output/${species}/exons.bed
find "${DIR}/output/${species}" -type f -name "*::introns.bed" -exec cat {} \; | sortBed > ${DIR}/output/${species}/introns.bed
find "${DIR}/output/${species}" -type f -name "*::transcript.bed" -exec cat {} \; | sortBed > ${DIR}/output/${species}/transcripts.bed

# Extract intron sequence

${DIR}/bin/extractIntronSequence.sh ${DIR} ${species} ${assembly} ${source} ${AssemblyFile}

# Clean-up files

find "${DIR}/output/${species}" -type f -name "*::*.bed" -exec rm -f {} \;
rm -f ${DIR}/output/${species}/*.tmp
