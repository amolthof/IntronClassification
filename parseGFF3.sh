#!/bin/bash

# Set base path

DIR=$1

# Set species info

species=$2
assembly=$3
source=$4
AnnotationFile=$5

# Extract exons

awk -F"\t" -v OFS="\t" '{if ($3=="exon") print $1,$4,$5,($5-$4+1),$7,$NF}' ${AnnotationFile} | sort | uniq > ${DIR}/output/${species}/exon.tmp

# Check if exons exist

if [ ! -s ${DIR}/output/${species}/exon.tmp ]; then

	exit

fi

# Get coordinates / parse info
# fields = chr, start, stop, len, strand, [parsed]

echo "chr" "start" "stop" "len" "strand" | sed 's/ /\t/g' > ${DIR}/output/${species}/exon_coords.tmp
cut -f1-5 ${DIR}/output/${species}/exon.tmp >> ${DIR}/output/${species}/exon_coords.tmp
cut -f6 ${DIR}/output/${species}/exon.tmp | awk '{print $0";"}' | sed 's/=/ "/g' | sed 's/;/"; /g' | sed 's/ID /feature_id /' | sed 's/Parent /transcript_id /' | python3 ${DIR}/bin/parseInfo.py > ${DIR}/output/${species}/exon_parsed.tmp

paste ${DIR}/output/${species}/exon_coords.tmp ${DIR}/output/${species}/exon_parsed.tmp > ${DIR}/output/${species}/exon.info

rm -f ${DIR}/output/${species}/*.tmp