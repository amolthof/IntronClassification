#!/bin/bash

# Set base path

DIR=$1

# Set species info

species=$2
assembly=$3
source=$4
AssemblyFile=$5

# Extract unique intron FASTA (from -3 at 5'ss to +1 at 3'ss)
# intron_key = ${species}_${chr}_${start}_${stop}_${strand}

awk -F"\t" -v OFS="\t" -v name=${species} '{if ($6=="+") print $1,($2-4),($3+1),name"_"$1"_"$2"_"$3"_F",($3-$2+5),$6; else print $1,($2-2),($3+3),name"_"$1"_"$2"_"$3"_R",($3-$2+5),$6}' ${DIR}/output/${species}/introns.bed | sort | uniq > ${DIR}/output/${species}/introns_getfasta.bed.tmp

# Use bedtools getfasta to extract the entire intronic sequence

bedtools getfasta -fi ${AssemblyFile} -bed ${DIR}/output/${species}/introns_getfasta.bed.tmp -s -name | fasta_formatter -t > ${DIR}/output/${species}/introns_raw.fa.tmp

# Extract 5'ss (-3 to +9), 3'ss (-13 to +1), potential U12 BPS (-40 to -1) window, potential U2 BPS (-44 to -18) window
# Filter sequence error introns (uncalled bases marked by 'N') and length errors

awk -F"\t" -v OFS="\t" '{len=length($2); print $1,toupper(substr($2,1,12)),toupper(substr($2,len-13,14)),toupper(substr($2,len-40,40)),toupper(substr($2,len-44,27)),toupper($2)}' ${DIR}/output/${species}/introns_raw.fa.tmp | awk -F"\t" -v OFS="\t" '{if (index($2,"N")==0 && index($3,"N")==0 && index($4,"N")==0 && length($2)==12 && length($3)==14 && length($4)==40) print $0}' | sed 's/(+)//' | sed 's/(-)//' | sort -k1,1 > ${DIR}/IntronFASTAs/${species}.${assembly}.${source}-Introns.fa

# Clean-up

rm -f ${DIR}/output/${species}/*.tmp