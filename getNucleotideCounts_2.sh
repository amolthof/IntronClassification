#!/bin/bash

# Set inputs: working directory, source directory

DIR=$1		# Working directory

# Set species name, assembly, source, X, Y, titration

species=$2
assembly=$3
source=$4

t=$5
X=$6
Y=$7

# Set input file & intron FASTA
# FIELDS: intron_id, U2 5'ss, U2 5'ss score, U12 5'ss, U12 5'ss score, U2 3'ss, U2 3'ss score, U2 BPS, U2 BPS score, U12 BPS, U12 BPS score

INPUT=${DIR}/output/${species}/Initial/introns_scored.tsv

# Screen intron for putative major and minor introns based on the following criteria
	# U2_introns: U2_DS score is greater than 50 and greater than U12_DS, U2_AcceptorSite is greater than 50
	# U12_introns: U12_DS score is greater than 50 and greater than U2_DS by X, U2_DS is less than 50, U12_BPS is greater than Y

awk -F"\t" -v OFS="\t" -v x=${X} -v y=${Y} '{d=substr($2,3,2); a=substr($6,12,2); U2=$3; U12=$5; PPT=$7; BPS=$11} {if (U2>50 && U2>U12 && PPT>50) print $0,d"-"a,"U2"; else if (U12>50 && U12>(U2+x) && BPS>y) print $0,d"-"a,"U12"; else print $0,d"-"a,"discard"}' ${INPUT} | sort -k 1,1 > ${DIR}/output/${species}/introns_screened_1-${t}.tmp

# Extract 5'ss sequence

cut -f1-2 ${DIR}/IntronFASTAs/${species}.${assembly}.${source}-Introns.fa > ${DIR}/output/${species}/introns_5ss.tmp
join ${DIR}/output/${species}/introns_5ss.tmp ${DIR}/output/${species}/introns_screened_1-${t}.tmp | sed 's/ /\t/g' > ${DIR}/output/${species}/introns_screened_2-${t}.tmp

# Group by preliminary class

grep -F -w "U2" ${DIR}/output/${species}/introns_screened_2-${t}.tmp > ${DIR}/output/${species}/U2_introns-${t}.tmp
grep -F -w "U12" ${DIR}/output/${species}/introns_screened_2-${t}.tmp > ${DIR}/output/${species}/U12_introns-${t}.tmp
grep -F -w "discard" ${DIR}/output/${species}/introns_screened_2-${t}.tmp > ${DIR}/output/${species}/discard_introns-${t}.tmp

# Group by terminal dinucleotide

for type in U2 U12 discard; do

	grep -F -w "GT-AG" ${DIR}/output/${species}/${type}_introns-${t}.tmp > ${DIR}/output/${species}/${type}_GT-AG_introns-${t}.tmp
	grep -F -w "GC-AG" ${DIR}/output/${species}/${type}_introns-${t}.tmp > ${DIR}/output/${species}/${type}_GC-AG_introns-${t}.tmp
	grep -F -w "AT-AC" ${DIR}/output/${species}/${type}_introns-${t}.tmp > ${DIR}/output/${species}/${type}_AT-AC_introns-${t}.tmp
	grep -F -w -v "GT-AG" ${DIR}/output/${species}/${type}_introns-${t}.tmp | grep -F -w -v "GC-AG" | grep -F -w -v "AT-AC" > ${DIR}/output/${species}/${type}_other_introns-${t}.tmp

done

# PWM 1: U2_DonorSite PWM, -2 to +6, 8nt sequence from U2 introns

awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1))}' ${DIR}/output/${species}/U2_introns-${t}.tmp > ${DIR}/output/${species}/U2_DonorSite-${t}.tmp

for Pos in {1..8}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/U2_DonorSite-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-U2_DonorSite.tsv

done

# PWM 2: U12_DonorSite PWM, +4 to +9, 6nt from U12 introns

if [ -s ${DIR}/output/${species}/U12_introns-${t}.tmp ]; then

	awk -F"\t" -v OFS="\t" '{print (substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1))}' ${DIR}/output/${species}/U12_introns-${t}.tmp > ${DIR}/output/${species}/U12_DonorSite-${t}.tmp

	for Pos in {1..6}; do

		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/U12_DonorSite-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-U12_DonorSite.tsv

	done

fi

# PWMs 3-6: sub-type specific, full donor site, 11nt, -2 to +9

for type in U2_GT-AG U2_GC-AG U12_GT-AG U12_AT-AC; do

	if [ ! -s ${DIR}/output/${species}/${type}_introns-${t}.tmp ]; then

		continue

	fi

	# Create position matrix of nucleotides at 5'ss from -2 to +9

	awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1))}' ${DIR}/output/${species}/${type}_introns-${t}.tmp > ${DIR}/output/${species}/${type}_DonorSite-${t}.tmp

	# Count nucleotides at each position

	for Pos in {1..11}; do

		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/${type}_DonorSite-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-${type}_DonorSite.tsv

	done

done 

# PWM 7: major acceptor site, 13nt, -13 to -1

awk -F"\t" -v OFS="\t" '{print (substr($7,1,1)),(substr($7,2,1)),(substr($7,3,1)),(substr($7,4,1)),(substr($7,5,1)),(substr($7,6,1)),(substr($7,7,1)),(substr($7,8,1)),(substr($7,9,1)),(substr($7,10,1)),(substr($7,11,1)),(substr($7,12,1)),(substr($7,13,1))}' ${DIR}/output/${species}/U2_introns-${t}.tmp > ${DIR}/output/${species}/U2_AcceptorSite-${t}.tmp

for Pos in {1..13}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/U2_AcceptorSite-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-U2_AcceptorSite.tsv

done

# PWM 8: major BPS, 7nt, force A6

awk -F"\t" -v OFS="\t" '{print (substr($9,1,1)),(substr($9,2,1)),(substr($9,3,1)),(substr($9,4,1)),(substr($9,5,1)),(substr($9,6,1)),(substr($9,7,1))}' ${DIR}/output/${species}/U2_introns-${t}.tmp | awk -F"\t" '{if ($6=="A") print $0}' > ${DIR}/output/${species}/U2_BranchPoint_A6-${t}.tmp

for Pos in {1..7}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/U2_BranchPoint_A6-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-U2_BranchPoint.tsv

done

# PWM 9-10: minor BPS, 12nt, force A9 & A10

awk -F"\t" -v OFS="\t" '{print (substr($11,1,1)),(substr($11,2,1)),(substr($11,3,1)),(substr($11,4,1)),(substr($11,5,1)),(substr($11,6,1)),(substr($11,7,1)),(substr($11,8,1)),(substr($11,9,1)),(substr($11,10,1)),(substr($11,11,1)),(substr($11,12,1))}' ${DIR}/output/${species}/U12_introns-${t}.tmp | awk -F"\t" '{if ($9=="A") print $0}' > ${DIR}/output/${species}/U12_BranchPoint_A9-${t}.tmp
awk -F"\t" -v OFS="\t" '{print (substr($11,1,1)),(substr($11,2,1)),(substr($11,3,1)),(substr($11,4,1)),(substr($11,5,1)),(substr($11,6,1)),(substr($11,7,1)),(substr($11,8,1)),(substr($11,9,1)),(substr($11,10,1)),(substr($11,11,1)),(substr($11,12,1))}' ${DIR}/output/${species}/U12_introns-${t}.tmp | awk -F"\t" '{if ($10=="A") print $0}' > ${DIR}/output/${species}/U12_BranchPoint_A10-${t}.tmp

for nt in A9 A10; do

	for Pos in {1..12}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/U12_BranchPoint_${nt}-${t}.tmp >> ${DIR}/PWMs/${t}/NucleotideCounts/${species}-U12_BranchPoint_${nt}.tsv

	done

done

# Echo stats

total=`awk 'END {print NR}' ${INPUT}`
major=`awk 'END {print NR}' ${DIR}/output/${species}/U2_introns-${t}.tmp`
major_1=`awk 'END {print NR}' ${DIR}/output/${species}/U2_GT-AG_introns-${t}.tmp`
major_2=`awk 'END {print NR}' ${DIR}/output/${species}/U2_GC-AG_introns-${t}.tmp`
major_3=`awk 'END {print NR}' ${DIR}/output/${species}/U2_AT-AC_introns-${t}.tmp`
major_4=`awk 'END {print NR}' ${DIR}/output/${species}/U2_other_introns-${t}.tmp`

if [ -s ${DIR}/output/${species}/U12_introns-${t}.tmp ]; then

	minor=`awk 'END {print NR}' ${DIR}/output/${species}/U12_introns-${t}.tmp`
	minor_1=`awk 'END {print NR}' ${DIR}/output/${species}/U12_GT-AG_introns-${t}.tmp`
	minor_2=`awk 'END {print NR}' ${DIR}/output/${species}/U12_GC-AG_introns-${t}.tmp`
	minor_3=`awk 'END {print NR}' ${DIR}/output/${species}/U12_AT-AC_introns-${t}.tmp`
	minor_4=`awk 'END {print NR}' ${DIR}/output/${species}/U12_other_introns-${t}.tmp`

else

	minor=0
	minor_1=0
	minor_2=0
	minor_3=0
	minor_4=0

fi

discard=`awk 'END {print NR}' ${DIR}/output/${species}/discard_introns-${t}.tmp`
discard_1=`awk 'END {print NR}' ${DIR}/output/${species}/discard_GT-AG_introns-${t}.tmp`
discard_2=`awk 'END {print NR}' ${DIR}/output/${species}/discard_GC-AG_introns-${t}.tmp`
discard_3=`awk 'END {print NR}' ${DIR}/output/${species}/discard_AT-AC_introns-${t}.tmp`
discard_4=`awk 'END {print NR}' ${DIR}/output/${species}/discard_other_introns-${t}.tmp`

echo "species:" ${species} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "assembly:" ${assembly} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "source:" ${source} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "" >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "${t}_PWM_GROUP" "COUNT" "CONDITION" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "total_introns:" ${total} "all" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "putative_major:" ${major} "U2>50,U2>U12,PPT>50" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GT-AG:" ${major_1}  "major,GT-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GC-AG:" ${major_2} "major,GC-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...AT-AC:" ${major_3} "major,AT-AC" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...other:" ${major_4} "major,other" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "putative_minor:" ${minor} "U12>50,U12>U2+${X},BPS>${Y}" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GT-AG:" ${minor_1} "minor,GT-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GC-AG:" ${minor_2} "minor,GC-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...AT-AC:" ${minor_3} "minor,AT-AC"| sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...other:" ${minor_4} "minor,other" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "discard:" ${discard} "other" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GT-AG:" ${discard_1} "discard,GT-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...GC-AG:" ${discard_2} "discard,GC-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...AT-AC:" ${discard_3} "discard,AT-AC" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "...other:" ${discard_4} "discard,other" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt
echo "" >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

# Get consensus sequence

for type in U2_DonorSite U2_GT-AG_DonorSite U2_GC-AG_DonorSite U12_DonorSite U12_GT-AG_DonorSite U12_AT-AC_DonorSite U2_AcceptorSite U2_BranchPoint U12_BranchPoint_A9 U12_BranchPoint_A10; do

	PWM=${DIR}/PWMs/${t}/NucleotideCounts/${species}-${type}.tsv

	if [ -s ${PWM} ]; then

		echo "${type}_${t}_PWM_count_consensus:" `${DIR}/bin/getConsensus.sh ${PWM}` | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

	else

		echo "${type}_${t}_PWM_count_consensus:" "NULL" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

	fi

done

echo "" >> ${DIR}/output/${species}/IntronClassificationStats-${t}.txt

rm -f ${DIR}/output/${species}/*-${t}.tmp