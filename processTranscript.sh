#!/bin/bash

# Set base path

DIR=$1

# Set species info

species=$2
assembly=$3
version=$4

# Set transcript ID

ID=$5

# Generate tmp working path

tmp=${DIR}/output/${species}/${ID}

# Extract exon coordinates

head -n 1 ${DIR}/output/${species}/exon.info > ${tmp}::exons.tmp
grep -F -w ${ID} ${DIR}/output/${species}/exon.info >> ${tmp}::exons.tmp

# Generate exon bed, ensure ID matches exactly

awk -F"\t" -v OFS="\t" -v id=${ID} 'NR==1 {for (i=1; i<=NF; i++) if ($i=="transcript_id") col=i} NR>1 && col {if ($col==id) print $1,$2,$3,$col,$4,$5}' ${tmp}::exons.tmp > ${tmp}::exons.bed.tmp

# Set chr, strand

chr=`cut -f1 ${tmp}::exons.bed.tmp | head -n 1`
strand=`cut -f6 ${tmp}::exons.bed.tmp | head -n 1`

# Determine transcript start/stop, len

start=`cut -f2 ${tmp}::exons.bed.tmp | sort -n | head -n 1`
stop=`cut -f3 ${tmp}::exons.bed.tmp | sort -n | tail -n 1`
len=`expr ${stop} - ${start} + 1`

# Generate transcript bed

echo ${chr} ${start} ${stop} ${ID} ${len} ${strand} | sed 's/ /\t/g' > ${tmp}::transcript.bed.tmp

# Reduce BED file coordinates to ensure subtractBed works & runs efficiently 

awk -F"\t" -v OFS="\t" -v n=${start} '{print $1,$2-n+1,$3-n+1,$4,$5,$6}' ${tmp}::exons.bed.tmp > ${tmp}::exons_reduced.bed.tmp
awk -F"\t" -v OFS="\t" -v n=${start} '{print $1,$2-n+1,$3-n+1,$4,$5,$6}' ${tmp}::transcript.bed.tmp > ${tmp}::transcript_reduced.bed.tmp

# Subtract exons from transcript to yield introns, fix coords

subtractBed -a ${tmp}::transcript_reduced.bed.tmp -b ${tmp}::exons_reduced.bed.tmp | awk -F"\t" -v OFS="\t" '{if (($3-1)>=($2+1)) print $1,($2+1),($3-1),$4,($3-$2-1),$6}' > ${tmp}::introns_reduced.bed.tmp

# Add start back to coords

awk -F"\t" -v OFS="\t" -v n=${start} '{print $1,$2+n-1,$3+n-1,$4,$5,$6}' ${tmp}::introns_reduced.bed.tmp > ${tmp}::introns.bed.tmp

# Check for 1nt exons (start = stop)

awk -F"\t" '{if ($2==$3) print $0}' ${tmp}::exons.bed.tmp > ${tmp}::exons_1nt.bed.tmp

if [ -s ${tmp}::exons_1nt.bed.tmp ]; then

	# Iterate over 1nt exons

	numExons=`awk 'END {print NR}' ${tmp}::exons_1nt.bed.tmp`

	for j in `seq 1 ${numExons}`; do

		# Extract exon bed

		head -n ${j} ${tmp}::exons_1nt.bed.tmp | tail -n 1 > ${tmp}::exon.bed.tmp

		# Set exon coord

		coord=`cut -f2 ${tmp}::exon.bed.tmp`

		# Fix affected introns

		awk -F"\t" -v OFS="\t" -v c=${coord} '{if (($2-2)==c) print $1,($2-1),$3,$4,($5+1),$6; else if (($3+2)==c) print $1,$2,($3+1),$4,($5+1),$6; else print $1,$2,$3,$4,$5,$6}' ${tmp}::introns.bed.tmp > ${tmp}::introns_fixed.bed.tmp
		mv ${tmp}::introns_fixed.bed.tmp ${tmp}::introns.bed.tmp

	done

fi

# Determine exon/intron rank

if [ ${strand} == "+" ]; then

	sort -k2,2n ${tmp}::exons.bed.tmp | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"EXON-"NR"::"$4,$5,$6}' > ${tmp}::exons.bed
	sort -k2,2n ${tmp}::introns.bed.tmp | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"INTRON-"NR"::"$4,$5,$6}' > ${tmp}::introns.bed

else

	sort -k2,2nr ${tmp}::exons.bed.tmp | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"EXON-"NR"::"$4,$5,$6}' > ${tmp}::exons.bed
	sort -k2,2nr ${tmp}::introns.bed.tmp | awk -F"\t" -v OFS="\t" '{print $1,$2,$3,"INTRON-"NR"::"$4,$5,$6}' > ${tmp}::introns.bed

fi

# Capture transcript entry

cat ${tmp}::transcript.bed.tmp > ${tmp}::transcript.bed

# Clean-up tmp path

rm -f ${tmp}::*.tmp