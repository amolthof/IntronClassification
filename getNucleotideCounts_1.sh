#!/bin/bash

# Set base path

DIR=$1

# Set species info

species=$2
assembly=$3
source=$4

IntronFASTA=${DIR}/IntronFASTAs/${species}.${assembly}.${source}-Introns.fa

# Subtype introns by terminal dinucleotides: group GT.GC-AG introns and AT-AC introns

awk -F"\t" -v OFS="\t" '{d=substr($2,4,2); a=substr($3,12,2); if ((d=="GT" || d=="GC") && a=="AG") print $0}' ${IntronFASTA} > ${DIR}/output/${species}/GT.GC-AG_introns.tmp
awk -F"\t" -v OFS="\t" '{d=substr($2,4,2); a=substr($3,12,2); if (d=="AT" && a=="AC") print $0}' ${IntronFASTA} > ${DIR}/output/${species}/AT-AC_introns.tmp

# PWM 1: Get putative major 5'ss, -2 to +6, 8nt sequence from GT.GC-AG introns

awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1))}' ${DIR}/output/${species}/GT.GC-AG_introns.tmp > ${DIR}/output/${species}/GT.GC-AG_DonorSite.tmp

for Pos in {1..8}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/GT.GC-AG_DonorSite.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U2_DonorSite.tsv

done

# PWM 2: Generate minor 5'ss, +4 to +9, 6nt from AT-AC introns

if [ -s ${DIR}/output/${species}/AT-AC_introns.tmp ]; then

	awk -F"\t" -v OFS="\t" '{print (substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1))}' ${DIR}/output/${species}/AT-AC_introns.tmp > ${DIR}/output/${species}/AT-AC_DonorSite.tmp

	for Pos in {1..6}; do

		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/AT-AC_DonorSite.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U12_DonorSite.tsv

	done

fi

# PWM 3: Generate major 3'ss, -13 to -1, 13nt sequence from GT.GC-AG introns

awk -F"\t" -v OFS="\t" '{print (substr($3,1,1)),(substr($3,2,1)),(substr($3,3,1)),(substr($3,4,1)),(substr($3,5,1)),(substr($3,6,1)),(substr($3,7,1)),(substr($3,8,1)),(substr($3,9,1)),(substr($3,10,1)),(substr($3,11,1)),(substr($3,12,1)),(substr($3,13,1))}' ${DIR}/output/${species}/GT.GC-AG_introns.tmp > ${DIR}/output/${species}/GT.GC-AG_AcceptorSite.tmp

for Pos in {1..13}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/GT.GC-AG_AcceptorSite.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U2_AcceptorSite.tsv

done

# PWM 4: Extract highest scoring 7nt major BPS between -44 to -18 from GT.GC-AG introns

# Extract GT.GC-AG introns

join ${DIR}/output/${species}/GT.GC-AG_introns.tmp ${IntronFASTA} > ${DIR}/output/${species}/GT.GC-AG_introns_3ss.tmp

awk -v OFS="\t" '{print (substr($5,1,1)),(substr($5,2,1)),(substr($5,3,1)),(substr($5,4,1)),(substr($5,5,1)),(substr($5,6,1)),(substr($5,7,1)),(substr($5,8,1)),(substr($5,9,1)),(substr($5,10,1)),(substr($5,11,1)),(substr($5,12,1)),(substr($5,13,1)),(substr($5,14,1)),(substr($5,15,1)),(substr($5,16,1)),(substr($5,17,1)),(substr($5,18,1)),(substr($5,19,1)),(substr($5,20,1)),(substr($5,21,1)),(substr($5,22,1)),(substr($5,23,1)),(substr($5,24,1)),(substr($5,25,1)),(substr($5,26,1)),(substr($5,27,1)),$1}' ${DIR}/output/${species}/GT.GC-AG_introns_3ss.tmp > ${DIR}/output/${species}/GT.GC-AG_BranchPoint.tmp

# Generate all potential 7 nt BPS (27 nt, 21 total 7nt sequences) with A at +6

awk -F"\t" -v OFS="\t" '{n=NF-6} {for (i=1;i<n;i++) {print $i,$(i+1),$(i+2),$(i+3),$(i+4),$(i+5),$(i+6),$NF}}' ${DIR}/output/${species}/GT.GC-AG_BranchPoint.tmp | awk -F"\t" '{if ($6=="A") print $0}' > ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6.tmp

# Extract GT.GC-AG introns



for Pos in {1..7}; do

	awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U2_PotentialBPS.tsv

done

# PWM 5: Extract highest scoring 12nt minor BPS, A9 & A10, between -40 to -1 from AT-AC introns

if [ -s ${DIR}/output/${species}/AT-AC_introns.tmp ]; then

	awk -F"\t" -v OFS="\t" '{print (substr($4,1,1)),(substr($4,2,1)),(substr($4,3,1)),(substr($4,4,1)),(substr($4,5,1)),(substr($4,6,1)),(substr($4,7,1)),(substr($4,8,1)),(substr($4,9,1)),(substr($4,10,1)),(substr($4,11,1)),(substr($4,12,1)),(substr($4,13,1)),(substr($4,14,1)),(substr($4,15,1)),(substr($4,16,1)),(substr($4,17,1)),(substr($4,18,1)),(substr($4,19,1)),(substr($4,20,1)),(substr($4,21,1)),(substr($4,22,1)),(substr($4,23,1)),(substr($4,24,1)),(substr($4,25,1)),(substr($4,26,1)),(substr($4,27,1)),(substr($4,28,1)),(substr($4,29,1)),(substr($4,30,1)),(substr($4,31,1)),(substr($4,32,1)),(substr($4,33,1)),(substr($4,34,1)),(substr($4,35,1)),(substr($4,36,1)),(substr($4,37,1)),(substr($4,38,1)),(substr($4,39,1)),(substr($4,40,1)),$1}' ${DIR}/output/${species}/AT-AC_introns.tmp > ${DIR}/output/${species}/AT-AC_BranchPoint.tmp

	# Generate all potential 12nt BPS (40 nt, 29 total 12nt sequences) with A at +9 and/or +10

	awk -F"\t" -v OFS="\t" '{n=NF-11} {for (i=1;i<n;i++) {print $i,$(i+1),$(i+2),$(i+3),$(i+4),$(i+5),$(i+6),$(i+7),$(i+8),$(i+9),$(i+10),$(i+11),$NF}}' ${DIR}/output/${species}/AT-AC_BranchPoint.tmp | awk -F"\t" '{if ($9=="A" || $10=="A") print $0}' > ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10.tmp

	for Pos in {1..12}; do

		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U12_PotentialBPS.tsv

	done

fi

# Echo species stats

total=`awk 'END {print NR}' ${IntronFASTA}`
major=`awk 'END {print NR}' ${DIR}/output/${species}/GT.GC-AG_introns.tmp`
major_BPS=`awk 'END {print NR}' ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6.tmp`
other=`expr ${total} - ${major}`

if [ -s ${DIR}/output/${species}/AT-AC_introns.tmp ]; then

	minor=`awk 'END {print NR}' ${DIR}/output/${species}/AT-AC_introns.tmp`
	minor_BPS=`awk 'END {print NR}' ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10.tmp`
	other=`expr ${other} - ${minor}`

else

	minor=0
	minor_BPS=0

fi

echo "species:" ${species} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "assembly:" ${assembly} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "source:" ${source} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "INITIAL_PWM_GROUP" "COUNT" "CONDITION" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "total_introns:" ${total} "all" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "putative_major:" ${major} "GT.GC-AG" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "...PotentialBPS:" ${major_BPS} "A6" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "putative_minor:" ${minor} "AT-AC" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "...PotentialBPS:" ${minor_BPS} "A9xA10" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "discard:" ${other} "otherwise" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

# Get consensus sequence

for type in U2_DonorSite U12_DonorSite U2_AcceptorSite U2_PotentialBPS U12_PotentialBPS; do

	PWM=${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv

	if [ -s ${PWM} ]; then

		echo "${type}_initial_PWM_count_consensus:" `${DIR}/bin/getConsensus.sh ${PWM}` | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

	else

		echo "${type}_initial_PWM_count_consensus:" "NULL" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

	fi

done

echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt