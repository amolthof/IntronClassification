#!/bin/bash

# Set inputs: working directory, list of species

DIR=$1			# Working directory
SpeciesList=$2		# List of species

# Set titration, num species

t=$3

numSpecies=`awk 'END {print NR}' ${SpeciesList}`

# Iterate over calibration list

for i in `seq 1 ${numSpecies}`; do

	i=`expr ${i} + 1`

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	source=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ -s ${DIR}/PipelineCalibration/${species}-ResponsiveIntrons.tsv ]; then

		${DIR}/bin/processSpecies.sh ${DIR} ${species} ${assembly} ${source} ${t}

	fi

done

# Calculate precision/recall

${DIR}/bin/calculatePrecisionRecall.sh ${DIR} ${SpeciesList} ${t}
