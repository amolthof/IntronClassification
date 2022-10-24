## Extracting 5'SS, BPS and 3'SS sequences from intron FASTA files.

The purpose of this bash script is to extract splice sites from all annotated introns.  

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
