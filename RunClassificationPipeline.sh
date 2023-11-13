#!/bin/bash

# Set working directory, species list, calibration list

DIR=/data1/charles/Anouk/IntronClassification_git/IntronClassification		# Set working directory
SpeciesList=${DIR}/MIDB_v2.0-SpeciesList.tsv					# SpeciesList

numSpecies=`awk 'END {print NR-1}' ${SpeciesList}`

# Generate output directories

for directory in output PWMs IntronFASTAs; do

	if [ ! -d ${DIR}/${directory} ]; then

		mkdir ${DIR}/${directory}

	fi

done

echo "STARTING INTRON CLASSIFICATION PIPELINE"
echo "Last updated: 04-Nov-2023"

# PART 1: data extraction

echo ""
echo "PART 1: Extracting intron FASTA for ${numSpecies} species |" `date`
echo ""

for i in `seq 1 ${numSpecies}`; do

	i=`expr ${i} + 1`

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	source=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

	AssemblyLink=`cut -f 5 ${SpeciesList} | head -n ${i} | tail -n 1`
	AnnotationLink=`cut -f 6 ${SpeciesList} | head -n ${i} | tail -n 1`
	FileType=`cut -f 7 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ ! -d ${DIR}/output/${species} ]; then

		mkdir ${DIR}/output/${species}

	fi

	IntronFASTA=${DIR}/IntronFASTAs/${species}.${assembly}.${source}-Introns.fa

	if [ ! -s ${IntronFASTA} ]; then

		i=`expr ${i} - 1`

		echo "Species ${i} of ${numSpecies}: ${species} - extracting intron FASTA |" `date`

		# Download & unzip source files

		wget --quiet ${AssemblyLink} -O ${DIR}/output/${species}/assembly.gz
		gunzip ${DIR}/output/${species}/assembly.gz &> /dev/null
		mv ${DIR}/output/${species}/assembly* ${DIR}/output/${species}/${species}.${assembly}.${source}-genome.fa

		wget --quiet ${AnnotationLink} -O ${DIR}/output/${species}/annotation.gz
		gunzip ${DIR}/output/${species}/annotation.gz &> /dev/null
		mv ${DIR}/output/${species}/annotation* ${DIR}/output/${species}/${species}.${assembly}.${source}-genome.${FileType}

		# Process source files

		${DIR}/bin/GenerateIntronFASTA.sh ${DIR} ${species} ${assembly} ${source} ${FileType}

		gzip ${DIR}/output/${species}/${species}.${assembly}.${source}-genome.fa &> /dev/null &
		gzip ${DIR}/output/${species}/${species}.${assembly}.${source}-genome.${FileType} &> /dev/null &

	else

		echo "Species ${i} of ${numSpecies}: ${species} - introns already extracted |" `date`

	fi

done

# PART 2: initial GT.GC-AG/AT-AC grouping & nt counting across species

if [ ! -d ${DIR}/PWMs/Initial ]; then

	echo ""
	echo "PART 2: Getting nucleotide counts at GT.GC-AG/AT-AC splice-sites across all ${numSpecies} species |" `date`
	echo ""

	mkdir ${DIR}/PWMs/Initial
	mkdir ${DIR}/PWMs/Initial/NucleotideCounts

	for i in `seq 1 ${numSpecies}`; do

		i=`expr ${i} + 1`

		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
		assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
		source=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`
		group=`cut -f 8 ${SpeciesList} | head -n ${i} | tail -n 1`

		if [ ${group} == "263" ]; then
		
			${DIR}/bin/getNucleotideCounts_1.sh ${DIR} ${species} ${assembly} ${source} &
			head -n ${i} ${SpeciesList} | tail -n 1 >> ${DIR}/MIDB_v2.0-SpeciesList_PWMs.tsv

		fi
	
	done

	# PAUSE: wait for all species to be processed

	wait

	# Merge PWMS

	${DIR}/bin/mergePWMs_1.sh ${DIR} ${DIR}/MIDB_v2.0-SpeciesList_PWMs.tsv

else

	echo ""
	echo "PART 2: GT.GC-AG/AT-AC splice-site nucleotides already counted |" `date`
	echo ""

fi

# PART 3: initial scoring & BPS extraction

echo ""
echo "PART 3: Scoring introns against initial PWMs |" `date`
echo ""

for i in `seq 1 ${numSpecies}`; do

	i=`expr ${i} + 1`

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	source=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`
	group=`cut -f 8 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ ! -d ${DIR}/output/${species}/Initial ]; then

		mkdir ${DIR}/output/${species}/Initial

		${DIR}/bin/scoreInitialPWMs.sh ${DIR} ${species} ${assembly} ${source} &

	fi

done

# PAUSE: wait for all species to be processed

wait

# PART 4: generate titration matrix, 169 total

if [ ! -s ${DIR}/TitrationMatrix.tsv ]; then

	echo ""
	echo "PART 4: Generating titration matrix |" `date`
	echo ""

	# U12 5'SS > U2 5'SS + X; for X = { 0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, }

	X=0

	for x in A B C D E F G H I J K L M; do

		# U12 BPS > Y; for Y = { 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0 }

		Y=50

		for y in 1 2 3 4 5 6 7 8 9 10 11 12 13; do

			t=${x}${y} 

			echo ${t} ${X} ${Y} | sed 's/ /\t/g' >> ${DIR}/TitrationMatrix.tsv

			Y=`awk -v y=${Y} 'BEGIN {print y+2.5}'`

		done

		X=`awk -v x=${X} 'BEGIN {print x+2.5}'`

	done

else

	echo ""
	echo "PART 4: Titration matrix generated |" `date`
	echo ""

fi

# PART 5: get nucleotide counts for each titration

numTitrations=`awk 'END {print NR}' ${DIR}/TitrationMatrix.tsv`

echo ""
echo "PART 5: Get nucleotide counts for ${numTitrations} titrations |" `date`
echo ""

for i in `seq 1 ${numTitrations}`; do

	# Set titration parameters

	t=`cut -f 1 ${DIR}/TitrationMatrix.tsv | head -n ${i} | tail -n 1`
	X=`cut -f 2 ${DIR}/TitrationMatrix.tsv | head -n ${i} | tail -n 1`
	Y=`cut -f 3 ${DIR}/TitrationMatrix.tsv | head -n ${i} | tail -n 1`

	if [ ! -d ${DIR}/PWMs/${t} ]; then

		echo "Titration ${t}: 5'ss + ${X}, BPS > ${Y} | PROCESSING"

		mkdir ${DIR}/PWMs/${t}
		mkdir ${DIR}/PWMs/${t}/NucleotideCounts

		for j in `seq 1 ${numSpecies}`; do

			species=`cut -f 1 ${SpeciesList} | head -n ${j} | tail -n 1`
			assembly=`cut -f 2 ${SpeciesList} | head -n ${j} | tail -n 1`
			source=`cut -f 3 ${SpeciesList} | head -n ${j} | tail -n 1`
			group=`cut -f 8 ${SpeciesList} | head -n ${i} | tail -n 1`

			if [ ${group} == "263" ]; then
	
				${DIR}/bin/getNucleotideCounts_2.sh ${DIR} ${species} ${assembly} ${source} ${t} ${X} ${Y} &

			fi
	
		done

		# PAUSE: wait for all species to be processed

		wait

		# Merge PWMS

		${DIR}/bin/mergePWMs_2.sh ${DIR} ${DIR}/MIDB_v2.0-SpeciesList_PWMs.tsv ${t}

		# Run titration

		if [ -d ${DIR}/PipelineCalibration ]; then

			${DIR}/bin/runTitration.sh ${DIR} ${SpeciesList} ${t} &> ${DIR}/PipelineCalibration/RunTitration-${t}.log &

		fi

	else

		echo "Titration ${t}: 5'ss + ${X}, BPS > ${Y} | COMPLETE"

	fi

done

# PAUSE: wait for all titrations to be processed

wait

# Determine titration by max F1 score

if [ ! -d ${DIR}/PipelineCalibration ]; then

	echo "PipelineCalibration not found - cannot determine best titration. Exiting..."

	exit

fi

t=`grep -F -w "RESPONSIVE" ${DIR}/PipelineCalibration/PipelineCalibration-*.tsv | cut -f3,9 | sort -k1,1 | key-merge | awk '{sum=0; for (i=2;i<=NF;i++) sum+=$i; print $1,sum/(NF-1)}' | sort -k2,2n | tail -n1 | cut -f1 -d' '`

echo ""
echo "PART 9: Deploying ${numSpecies} for scoring & binning against ${t} PWMs |" `date`
echo ""

for i in `seq 1 ${numSpecies}`; do

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	source=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ ! -d ${DIR}/output/${species}/${t} ]; then

		echo "Deploying species ${i} of ${numSpecies}: ${species}, ${t}"

		mkdir ${DIR}/output/${species}/${t}

		${DIR}/bin/processSpecies.sh ${DIR} ${species} ${assembly} ${source} ${t} &

	fi

done

echo ""
echo "INTRON CLASSIFICATION PIPELINE FULLY DEPLOYED |" `date`
echo ""
