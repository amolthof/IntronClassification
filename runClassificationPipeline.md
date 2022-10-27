## Intron Classification Pipeline

The purpose of this bash script is to run the intron classification pipeline and bin all annotated introns based on PWM scores.  

Dependencies:

awk: https://www.gnu.org/software/gawk/gawk.html <br>
___

**Set inputs: working directory, source directory, source path name, list of species in source**

    workdir="/path/to/my/working_directory"				
    fastadir="/path/to/my/fasta_directory"				
    SpeciesList="/path/to/my/speciesList.tsv			# List of species in fasta directory = ${species}, ${assembly}, ${version}

**Make output directories**

	if [ ! -d ${workdir}/output ]; then
		mkdir ${workdir}/output
	fi

	if [ ! -d ${workdir}/PWMs ]; then
		mkdir ${workdir}/PWMs
	fi

	numSpecies=`awk 'END {print NR}' ${SpeciesList}`

	echo "STARTING INTRON CLASSIFICATION PIPELINE"

**PART 1: data extraction**

	echo ""
	echo "PART 1: Extracting intron FASTA for ${numSpecies} species |" `date`
	echo ""

	for i in `seq 1 ${numSpecies}`; do

		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
		assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
		version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`
		IntronFASTA=${fastadir}/IntronFASTAs/${species}.${assembly}.${version}.introns.fa

		if [ ! -d ${workdir}/output/${species} ]; then

			echo "Species ${i} of ${numSpecies}: ${species} - extracting introns |" `date`
			mkdir ${workdir}/output/${species}

			awk -F"\t" -v name=${species} '{if ($1==name) print $0}' ${SpeciesList} >> ${workdir}/speciesList_classify.tsv
			sort -k1,1 ${IntronFASTA} > ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns.fa

			IntronFASTA=${workdir}/output/${species}/${species}.${assembly}.${version}-Introns.fa

			# Get 5'ss sequence FASTA file
			cut -f 1-2 ${IntronFASTA} | sort -k 1,1 > ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa

			# Get 3'ss sequence FASTA file
			${workdir}/bin/extractIntronSequence_U2-BPS.sh ${workdir} ${fastadir} ${species} ${assembly} ${version}
		else
			echo "Species ${i} of ${numSpecies}: ${species} - introns already extracted |" `date`
		fi

	done

	sort -k1,1 ${workdir}/speciesList_classify.tsv | uniq > ${workdir}/speciesList_classify.tmp
	mv ${workdir}/speciesList_classify.tmp ${workdir}/speciesList_classify.tsv

**PART 2: initial GT.GC-AG/AT-AC grouping & nt counting across species**

	if [ ! -d ${workdir}/PWMs/Initial ]; then

		echo ""
		echo "PART 2: Getting nucleotide counts at GT.GC-AG/AT-AC splice-sites across all ${numSpecies} species |" `date`
		echo ""

		mkdir ${workdir}/PWMs/Initial
		mkdir ${workdir}/PWMs/Initial/NucleotideCounts

		for i in `seq 1 ${numSpecies}`; do
			species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
			assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
			version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`
			${workdir}/bin/getNucleotideCounts_1.sh ${workdir} ${fastadir} ${species} ${assembly} ${version} &
		done

		echo ""
		echo "EXITING: Re-run RunClassificationPipeline.sh when getNucleotideCounts_1.sh is complete across all species |" `date`
		echo ""
		exit
	else
		echo ""
		echo "PART 2: GT.GC-AG/AT-AC splice-site nucleotides already counted |" `date`
		echo ""
	fi

**PART 3: merge initial PWMs across species**

	if [ ! -s ${workdir}/PWMs/Initial/${numSpecies}_Species-U2_DonorSite_LOD.tsv ]; then
		echo ""
		echo "PART 3: Generating initial PWMs across ${numSpecies} species |" `date`
		echo ""

		${workdir}/bin/mergePWMs_1.sh ${workdir} ${SpeciesList} 
	else
		echo ""
		echo "PART 3: Initial PWMs already generated |" `date`
		echo ""
	fi

**PART 4: initial scoring & BPS extraction**

	for i in `seq 1 ${numSpecies}`; do
		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
		assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
		version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

		if [ ! -d ${workdir}/output/${species}/Initial ]; then
			if [ ${i} -eq 1 ]; then
				echo ""
				echo "PART 4: Scoring introns against initial PWMs |" `date`
				echo ""
			fi

			mkdir ${workdir}/output/${species}/Initial
			${workdir}/bin/scoreInitialPWMs.sh ${workdir} ${SpeciesList} ${species} ${assembly} ${version} &

			if [ ${i} -eq ${numSpecies} ]; then
				echo ""
				echo "EXITING: Re-run RunClassificationPipeline.sh when scoreInitialPWMs.sh is complete across all species |" `date`
				echo ""
				exit
			fi
		else
			echo ""
			echo "PART 4: Initial scoring complete |" `date`
			echo ""
			break
		fi
	done

**PART 6: merge refined PWMs across species**

	if [ ! -s ${workdir}/PWMs/${numSpecies}_Species-U2_GT-AG_DonorSite_LOD.tsv ]; then
		echo ""
		echo "PART 6: Merge PWMs across species |" `date`
		echo ""
		
		nice -n 19 ${workdir}/bin/mergePWMs_2.sh ${workdir} ${SpeciesList} &> ${workdir}/PipelineCalibration/mergePWMs_2.log &
	else
		echo ""
		echo "PART 6: Refined PWMs already generated |" `date`
		echo ""
	fi

**PART 7: Score and bin all introns against refined PWMs**

	echo ""
	echo "PART 7: Deploy ${numSpecies} for scoring & binning against PWMs |" `date`
	echo ""

	for i in `seq 1 ${numSpecies}`; do
		species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
		assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
		version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`
		
		if [ ! -s ${workdir}/output/${species}/introns_binned.tsv ]; then
			echo "Deploying species ${i} of ${numSpecies}: ${species}"
			${workdir}/bin/processSpecies.sh ${workdir} ${SpeciesList} ${species} ${assembly} ${version} &
		else
			echo "Species ${i} of ${numSpecies}: ${species} - already processed"
		fi
	done

	echo ""
	echo "INTRON CLASSIFICATION PIPELINE FULLY DEPLOYED |" `date`
	echo ""
