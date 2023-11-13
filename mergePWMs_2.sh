#!/bin/bash

# Set inputs: working directory, source directory, list of species in source

DIR=$1		# Working directory
SpeciesList=$2	# List of species

t=$3

# Merge PWMs across all species for the given titration

numSpecies=`awk 'END {print NR}' ${SpeciesList}`

for type in U2_DonorSite U2_GT-AG_DonorSite U2_GC-AG_DonorSite  U12_DonorSite U12_GT-AG_DonorSite U12_AT-AC_DonorSite U2_AcceptorSite U2_BranchPoint U12_BranchPoint_A9 U12_BranchPoint_A10; do

	for i in `seq 1 ${numSpecies}`; do

		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`

		if [ ! -s ${DIR}/PWMs/${t}/NucleotideCounts/${species}-${type}.tsv ]; then

			if [ -s ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv ]; then

				echo "${type} ${t} PWM total count + ${species} consensus (NULL):"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv` "[${i} of ${numSpecies}]"

			else

				echo "${type} ${t} PWM total count (NULL) + ${species} consensus (NULL):" "NULL" "[${i} of ${numSpecies}]"

			fi

		elif [ -s ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv ]; then

			paste ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv ${DIR}/PWMs/${t}/NucleotideCounts/${species}-${type}.tsv | awk -F"\t" -v OFS="\t" '{print $1+$5,$2+$6,$3+$7,$4+$8}' > ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tmp
			mv ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tmp ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv

			echo "${type} ${t} PWM total count + ${species} consensus:"  `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		else

			cp ${DIR}/PWMs/${t}/NucleotideCounts/${species}-${type}.tsv ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv

			echo "${type} ${t} PWM total count (NULL) + ${species} consensus:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv` "[${i} of ${numSpecies}]"

		fi

	done

	# Create a new weighted matrix using the highest scoring BPS

	awk -F"\t" -v OFS="\t" '{for (i=1;i<=NF;i++) {accum+=$i; $i=$i+0}; print $0, accum; accum=0}' ${DIR}/PWMs/${t}/NucleotideCounts/${numSpecies}_Species-${type}.tsv | awk -F"\t" -v OFS="\t" '{A=($1/$NF+0.0001); C=($2/$NF+0.0001); G=($3/$NF+0.0001); T=($4/$NF+0.0001); print A, C, G, T}' > ${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_freq.tsv

	# Calculate LOD score

	awk -F"\t" -v OFS="\t" -v base=0.25 '{A=(log($1/base)/log(2)); C=(log($2/base)/log(2)); G=(log($3/base)/log(2)); T=(log($4/base)/log(2)); print A, C, G, T}' ${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_freq.tsv > ${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_LOD.tsv

	echo ""
	echo "${type} ${t} PWM total LOD:" `${DIR}/bin/getConsensus.sh ${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_LOD.tsv` "[${i} of ${numSpecies}]"
	echo ""

	PWM=${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_LOD.tsv

	# Determine number of positions in PWM

	n=`awk 'END {print NR}' ${PWM}`

	# Iterate over each position

	for i in `seq 1 ${n}`; do

		# Determine max and min LOD scores at the position
 
		max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if ($i>max) max=$i}; print max}' ${PWM} | head -n ${i} | tail -n 1`
		min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if ($i<min) min=$i}; print min}' ${PWM} | head -n ${i} | tail -n 1`

		# Scale each position wrt to max and min

		A=`awk -F"\t" -v max=${max} -v min=${min} '{if ($1==max) print 100; else if ($1==min) print 0; else if ($1<0) print 50*($1-min)/(0-min); else if ($1>0) print 50*$1/max+50; else print 50}' ${PWM} | head -n ${i} | tail -n 1`
		C=`awk -F"\t" -v max=${max} -v min=${min} '{if ($2==max) print 100; else if ($2==min) print 0; else if ($2<0) print 50*($2-min)/(0-min); else if ($2>0) print 50*$2/max+50; else print 50}' ${PWM} | head -n ${i} | tail -n 1`
		G=`awk -F"\t" -v max=${max} -v min=${min} '{if ($3==max) print 100; else if ($3==min) print 0; else if ($3<0) print 50*($3-min)/(0-min); else if ($3>0) print 50*$3/max+50; else print 50}' ${PWM} | head -n ${i} | tail -n 1`
		T=`awk -F"\t" -v max=${max} -v min=${min} '{if ($4==max) print 100; else if ($4==min) print 0; else if ($4<0) print 50*($4-min)/(0-min); else if ($4>0) print 50*$4/max+50; else print 50}' ${PWM} | head -n ${i} | tail -n 1`

		# Create scaled PWM

		echo ${A} ${C} ${G} ${T} | sed 's/ /\t/g' >> ${DIR}/PWMs/${t}/${numSpecies}_Species-${type}_scaled.tsv

	done

done