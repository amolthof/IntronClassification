#!/bin/bash

# Set base path

DIR=$1

# Set species list

SpeciesList=$2

# Merge initial PWMs across all species, extract AT-AC introns again

numSpecies=`awk 'END {print NR}' ${SpeciesList}`

for type in U2_DonorSite U12_DonorSite U2_AcceptorSite U2_PotentialBPS U12_PotentialBPS; do

	for i in `seq 1 ${numSpecies}`; do

		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`

		if [ ! -s ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv ]; then

			if [ -s ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ]; then

				echo "${type} initial PWM total count + ${species} consensus (NULL):"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

			else

				echo "${type} initial PWM total count (NULL) + ${species} consensus (NULL):" "NULL" "[${i} of ${numSpecies}]"

			fi

		elif [ -s ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ]; then

			paste ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv | awk -F"\t" -v OFS="\t" '{print $1+$5,$2+$6,$3+$7,$4+$8}' > ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tmp
			mv ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tmp ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv

			echo "${type} initial PWM total count + ${species} consensus:"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		else

			cp ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv
	
			echo "${type} initial PWM total count (NULL) + ${species} consensus:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		fi

	done

	# Create a new weighted matrix using the highest scoring BPS

	awk -F"\t" -v OFS="\t" '{for (i=1;i<=NF;i++) {accum+=$i; $i=$i+0}; print $0, accum; accum=0}' ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv | awk -F"\t" -v OFS="\t" '{A=($1/$NF+0.0001); C=($2/$NF+0.0001); G=($3/$NF+0.0001); T=($4/$NF+0.0001); print A, C, G, T}' > ${DIR}/PWMs/Initial/263_Species-${type}_freq.tsv

	# Calculate LOD score

	awk -F"\t" -v OFS="\t" -v base=0.25 '{A=(log($1/base)/log(2)); C=(log($2/base)/log(2)); G=(log($3/base)/log(2)); T=(log($4/base)/log(2)); print A, C, G, T}' ${DIR}/PWMs/Initial/263_Species-${type}_freq.tsv > ${DIR}/PWMs/Initial/263_Species-${type}_LOD.tsv

	echo ""
	echo "${type} initial total PWM LOD:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/263_Species-${type}_LOD.tsv` "[${i} of ${numSpecies}]"
	echo ""

done

# Extract highest scoring U2 BPS from GT.GC-AG introns in each species

PWM=${DIR}/PWMs/Initial/263_Species-U2_PotentialBPS_LOD.tsv

for Pos in {1..7}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

for i in `seq 1 ${numSpecies}`; do

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`

	echo "Species ${i} of ${numSpecies}: ${species} - extracting U2 BPS from GT.GC-AG introns |" `date`

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}

	{print $NF,X1,X2,X3,X4,X5,X6,X7,Y1+Y2+Y3+Y4+Y5+Y6+Y7}' ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6.tmp > ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6_scored.tmp
	
	# Extract the highest scoring BPS from each transcript with positive LOD

	${DIR}/bin/key-merge  ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6_scored.tmp > ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6_merged.tmp
	awk -v OFS="\t" '{max=$9; n=(NF-1)/8; for (i=1;i<=n;i++) {j=i*8+1; if (max<=$j) {max=$j; max_i=j}}; if (max>0) print $(max_i-7),$(max_i-6),$(max_i-5),$(max_i-4),$(max_i-3),$(max_i-2),$(max_i-1),$1}' ${DIR}/output/${species}/GT.GC-AG_BranchPoint_A6_merged.tmp > ${DIR}/output/${species}/HighestScoringBPS_A6.tmp

	# Extract nucleotide counts from each position

	for Pos in {1..7}; do

		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
										{if (count["C"]>0) C=count["C"]; else C=0}
										{if (count["G"]>0) G=count["G"]; else G=0}
										{if (count["T"]>0) T=count["T"]; else T=0}
										{print A, C, G, T}}' ${DIR}/output/${species}/HighestScoringBPS_A6.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U2_BranchPoint.tsv

	done

	# Echo PWM stats to file

	major_BPS=`awk 'END {print NR}' ${DIR}/output/${species}/HighestScoringBPS_A6.tmp`

	echo "INITIAL_PWM_GROUP" "COUNT" "CONDITION" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
	echo "U2_BranchPoint:" ${major_BPS} "Potential,LOD>0" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

done

# Extract highest scoring U2 BPS from AT-AC introns in each species

PWM=${DIR}/PWMs/Initial/263_Species-U12_PotentialBPS_LOD.tsv

for Pos in {1..12}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

for i in `seq 1 ${numSpecies}`; do

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ ! -s ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10.tmp ]; then

		PWM=${DIR}/PWMs/Initial/NucleotideCounts/${species}-U2_BranchPoint.tsv

		echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
		echo "U2_BranchPoint_initial_PWM_count_consensus:" `${DIR}/bin/getConsensus.sh ${PWM}` | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt
		echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

		continue		

	fi

	echo "Species ${i} of ${numSpecies}: ${species} - extracting U12 BPS from AT-AC introns |" `date`

	# Score each position, sum LOD score

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11; X12=$12}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}
	{if (X8=="A") Y8='${A8}'; else if (X8=="C") Y8='${C8}'; else if (X8=="G") Y8='${G8}'; else if (X8=="T") Y8='${T8}'}
	{if (X9=="A") Y9='${A9}'; else if (X9=="C") Y9='${C9}'; else if (X9=="G") Y9='${G9}'; else if (X9=="T") Y9='${T9}'}
	{if (X10=="A") Y10='${A10}'; else if (X10=="C") Y10='${C10}'; else if (X10=="G") Y10='${G10}'; else if (X10=="T") Y10='${T10}'}
	{if (X11=="A") Y11='${A11}'; else if (X11=="C") Y11='${C11}'; else if (X11=="G") Y11='${G11}'; else if (X11=="T") Y11='${T11}'}
	{if (X12=="A") Y12='${A12}'; else if (X12=="C") Y12='${C12}'; else if (X12=="G") Y12='${G12}'; else if (X12=="T") Y12='${T12}'}

	{print $NF,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12}' ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10.tmp > ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10_scored.tmp

	# Extract the highest scoring BPS & score from each transcript with positive LOD

	${DIR}/bin/key-merge ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10_scored.tmp > ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10_merged.tmp
	awk -v OFS="\t" '{max=$14; n=(NF-1)/13; for (i=1;i<=n;i++) {j=i*13+1; if (max<=$j) {max=$j; max_i=j}}; if (max>0) print $(max_i-12),$(max_i-11),$(max_i-10),$(max_i-9),$(max_i-8),$(max_i-7),$(max_i-6),$(max_i-5),$(max_i-4),$(max_i-3),$(max_i-2),$(max_i-1),$1,max}' ${DIR}/output/${species}/AT-AC_BranchPoint_A9xA10_merged.tmp > ${DIR}/output/${species}/HighestScoringBPS_A9xA10.tmp

	# Separate by A9 vs. A10

	awk -F"\t" -v OFS="\t" '{if ($9=="A") print $0}' ${DIR}/output/${species}/HighestScoringBPS_A9xA10.tmp > ${DIR}/output/${species}/HighestScoringBPS_A9.tmp
	awk -F"\t" -v OFS="\t" '{if ($10=="A") print $0}' ${DIR}/output/${species}/HighestScoringBPS_A9xA10.tmp > ${DIR}/output/${species}/HighestScoringBPS_A10.tmp

	# Extract nucleotide counts from each position for each group

	for nt in A9 A10; do

		if [ ! -s ${DIR}/output/${species}/HighestScoringBPS_${nt}.tmp ]; then

			# Echo PWM stats to file

			minor_BPS=0
			echo "putative_minor_BPS_${nt}:" "LOD>0" ${minor_BPS} | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

		else

			for Pos in {1..12}; do

				awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
												{if (count["C"]>0) C=count["C"]; else C=0}
												{if (count["G"]>0) G=count["G"]; else G=0}
												{if (count["T"]>0) T=count["T"]; else T=0}
												{print A, C, G, T}}' ${DIR}/output/${species}/HighestScoringBPS_${nt}.tmp >> ${DIR}/PWMs/Initial/NucleotideCounts/${species}-U12_BranchPoint_${nt}.tsv

			done

			# Echo PWM stats to file

			minor_BPS=`awk 'END {print NR}' ${DIR}/output/${species}/HighestScoringBPS_${nt}.tmp`
			echo "U2_BPS_${nt}:" ${minor_BPS} "PotentialBPS,LOD>0" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

		fi

	done

	echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

	# Get consensus sequence

	for type in U2_BranchPoint U12_BranchPoint_A9 U12_BranchPoint_A10; do

		PWM=${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv

		if [ -s ${PWM} ]; then

			echo "${type}_initial_PWM_count_consensus:" `${DIR}/bin/getConsensus.sh ${PWM}` | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

		else

			echo "${type}_initial_PWM_count_consensus:" "NULL" | sed 's/ /\t/g' >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

		fi

	done

	echo "" >> ${DIR}/output/${species}/IntronClassificationStats-Initial.txt

done

for type in U2_BranchPoint U12_BranchPoint_A9 U12_BranchPoint_A10; do

	for i in `seq 1 ${numSpecies}`; do

		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`

		if [ ! -s ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv ]; then

			if [ -s ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ]; then

				echo "${type} initial PWM total count + ${species} consensus (NULL):"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

			else

				echo "${type} initial PWM total count (NULL) + ${species} consensus (NULL):" "NULL" "[${i} of ${numSpecies}]"

			fi

		elif [ -s ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ]; then

			paste ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv | awk -F"\t" -v OFS="\t" '{print $1+$5,$2+$6,$3+$7,$4+$8}' > ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tmp
			mv ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tmp ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv

			echo "${type} initial PWM total count + ${species} consensus:"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		else

			cp ${DIR}/PWMs/Initial/NucleotideCounts/${species}-${type}.tsv ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv
	
			echo "${type} initial PWM total count (NULL) + ${species} consensus:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		fi

	done

	# Create a new weighted matrix using the highest scoring BPS

	awk -F"\t" -v OFS="\t" '{for (i=1;i<=NF;i++) {accum+=$i; $i=$i+0}; print $0, accum; accum=0}' ${DIR}/PWMs/Initial/NucleotideCounts/263_Species-${type}.tsv | awk -F"\t" -v OFS="\t" '{A=($1/$NF+0.0001); C=($2/$NF+0.0001); G=($3/$NF+0.0001); T=($4/$NF+0.0001); print A, C, G, T}' > ${DIR}/PWMs/Initial/263_Species-${type}_freq.tsv

	# Calculate LOD score

	awk -F"\t" -v OFS="\t" -v base=0.25 '{A=(log($1/base)/log(2)); C=(log($2/base)/log(2)); G=(log($3/base)/log(2)); T=(log($4/base)/log(2)); print A, C, G, T}' ${DIR}/PWMs/Initial/263_Species-${type}_freq.tsv > ${DIR}/PWMs/Initial/263_Species-${type}_LOD.tsv

	echo ""
	echo "${type} initial total PWM LOD:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/Initial/263_Species-${type}_LOD.tsv` "[${i} of ${numSpecies}]"
	echo ""

done