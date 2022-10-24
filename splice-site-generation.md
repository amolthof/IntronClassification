# Set inputs: working directory, source directory, source path name, list of species in source

DIR=/data1/charles/IntronClassification					# Working directory
dir=/data1/charles/MIDB_extended					# Source directory
SpeciesList=${DIR}/MIDB_extended_speciesList.tsv			# List of species in source directory = ${species}, ${assembly}, ${version}
CalibrationList=${DIR}/MIDB_extended_speciesList_calibration.tsv	# List of species to score across titrations

# Make output directories

if [ ! -d ${DIR}/output ]; then

	mkdir ${DIR}/output

fi

if [ ! -d ${DIR}/PWMs ]; then

	mkdir ${DIR}/PWMs

fi

numSpecies=`awk 'END {print NR}' ${SpeciesList}`
numSpeciesEval=`awk 'END {print NR}' ${CalibrationList}`

echo "STARTING INTRON CLASSIFICATION PIPELINE"
echo "Last updated: 30-Jul-2022, Charles Schwoerer"

# PART 1: data extraction

echo ""
echo "PART 1: Extracting intron FASTA for ${numSpecies} species |" `date`
echo ""

for i in `seq 1 ${numSpecies}`; do

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

	IntronFASTA=${dir}/IntronFASTAs/${species}.${assembly}.${version}.introns.fa

	if [ ${species} != "Ambystoma_mexicanum" ] && [ ${species} != "Physarum_polycephalum" ] && [ -s ${IntronFASTA} ]; then

		if [ ! -d ${DIR}/output/${species} ]; then

			echo "Species ${i} of ${numSpecies}: ${species} - extracting introns |" `date`

			mkdir ${DIR}/output/${species}

			awk -F"\t" -v name=${species} '{if ($1==name) print $0}' ${SpeciesList} >> ${DIR}/MIDB_extended_speciesList_classify.tsv
			sort -k1,1 ${IntronFASTA} > ${DIR}/output/${species}/${species}.${assembly}.${version}-Introns.fa

			IntronFASTA=${DIR}/output/${species}/${species}.${assembly}.${version}-Introns.fa

			# Get 5'ss sequence FASTA file

			cut -f 1-2 ${IntronFASTA} | sort -k 1,1 > ${DIR}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa

			# Get 3'ss sequence FASTA file

			${DIR}/bin/extractIntronSequence_U2-BPS.sh ${DIR} ${dir} ${species} ${assembly} ${version}

		else

			echo "Species ${i} of ${numSpecies}: ${species} - introns already extracted |" `date`

		fi

	fi

done

sort -k1,1 ${DIR}/MIDB_extended_speciesList_classify.tsv | uniq > ${DIR}/MIDB_extended_speciesList_classify.tmp
mv ${DIR}/MIDB_extended_speciesList_classify.tmp ${DIR}/MIDB_extended_speciesList_classify.tsv
