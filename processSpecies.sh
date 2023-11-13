#!/bin/bash

# Set inputs: working directory

DIR=$1			# Working directory

# Set species name and number of titrations

species=$2
assembly=$3
source=$4

t=$5

# Score introns

${DIR}/bin/scoreRefinedPWMs.sh ${DIR} ${species} ${assembly} ${source} ${t}

# Bin introns

${DIR}/bin/binIntronsAfterScoring.sh ${DIR} ${species} ${t}

# Get titration stats

total=`awk 'END {print NR}' ${DIR}/output/${species}/${t}/introns_binned.tsv`

awk -F"\t" '{if ($NF=="minor") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_minor.tmp
awk -F"\t" '{if ($NF=="minor_hybrid") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_minor_hybrid.tmp
awk -F"\t" '{if ($NF=="minor-like") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_minor_like.tmp

awk -F"\t" '{if ($NF=="major") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_major.tmp
awk -F"\t" '{if ($NF=="major_hybrid") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_major_hybrid.tmp
awk -F"\t" '{if ($NF=="major-like") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_major_like.tmp

awk -F"\t" '{if ($NF=="non-canonical") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv > ${DIR}/output/${species}/${t}/introns_non_canonical.tmp

# Echo stats to file

awk '{if (NR<=33) print $0}' ${DIR}/output/${species}/IntronClassificationStats-${t}.txt > ${DIR}/output/${species}/IntronClassificationStats-${t}.tmp
mv ${DIR}/output/${species}/IntronClassificationStats-${t}.tmp ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

echo "INTRON_CLASSIFICATION_${t}" "COUNT" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "total" ${total} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

for class in minor minor_hybrid minor_like major major_hybrid major_like non_canonical; do

	count=`cat ${DIR}/output/${species}/${t}/introns_${class}.tmp | wc -l`
	count_1=`grep -F -w "GT-AG" ${DIR}/output/${species}/${t}/introns_${class}.tmp | wc -l`
	count_2=`grep -F -w "GC-AG" ${DIR}/output/${species}/${t}/introns_${class}.tmp | wc -l`
	count_3=`grep -F -w "AT-AC" ${DIR}/output/${species}/${t}/introns_${class}.tmp | wc -l`
	count_4=`grep -F -w -v "GT-AG" ${DIR}/output/${species}/${t}/introns_${class}.tmp | grep -F -w -v "GC-AG" | grep -F -w -v "AT-AC" | wc -l`

	echo "${class}:" ${count} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
	echo "...GT-AG:" ${count_1} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
	echo "...GC-AG:" ${count_2} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
	echo "...AT-AC:" ${count_3} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
	echo "...other:" ${count_4} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

done

echo "" >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

rm -f ${DIR}/output/${species}/${t}/*.tmp

