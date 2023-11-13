#!/bin/bash

# Set inputs: working directory, list of species in source

DIR=$1			# Working directory
SpeciesList=$2		# List of species

# Set titration and num species

t=$3

numSpecies=`awk 'END {print NR}' ${SpeciesList}`

# Generate stat files

echo "species" "dataset" "titration" "minor" "minor_hybrid" "minor_like" "major" "major_hybrid" "major_like" "non-canonical" | sed 's/ /\t/g' > ${DIR}/PipelineCalibration/ResponsiveIntronCounts-${t}.tsv
echo "species" "dataset" "titration" "TP" "FP" "FN" "precision" "recall" "F1" "F2" | sed 's/ /\t/g' > ${DIR}/PipelineCalibration/PipelineCalibration-${t}.tsv

# Iterate over calibration list

for i in `seq 1 ${numSpecies}`; do

	i=`expr ${i} + 1`

	species=`cut -f 1 ${SpeciesList} | head -n ${i} | tail -n 1`
	assembly=`cut -f 2 ${SpeciesList} | head -n ${i} | tail -n 1`
	version=`cut -f 3 ${SpeciesList} | head -n ${i} | tail -n 1`

	if [ ! -s ${DIR}/PipelineCalibration/${species}-ResponsiveIntrons.tsv ]; then

		continue

	fi

	echo "GETTING RESPONSIVE INTRON STATS | ${species}, ${assembly}, ${version} | ${t} TITRATION |" `date`
	echo ""

	all_minor=`awk -F"\t" '{if ($NF=="minor") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	minor_hybrid=`awk -F"\t" '{if ($NF=="minor_hybrid") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	minor_like=`awk -F"\t" '{if ($NF=="minor-like") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	major=`awk -F"\t" '{if ($NF=="major") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	major_hybrid=`awk -F"\t" '{if ($NF=="major_hybrid") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	major_like=`awk -F"\t" '{if ($NF=="major-like") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`
	non_canonical=`awk -F"\t" '{if ($NF=="non-canonical") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tsv | wc -l`

	echo ${species} "ALL_INTRONS" ${t} ${all_minor} ${minor_hybrid} ${minor_like} ${major} ${major_hybrid} ${major_like} ${non_canonical} | sed 's/ /\t/g' >> ${DIR}/PipelineCalibration/ResponsiveIntronCounts-${t}.tsv

	# Set responsive introns input

	INPUT=${DIR}/PipelineCalibration/input/${species}-ResponsiveIntrons.tsv

	# Join to titration output

	join ${INPUT} ${DIR}/output/${species}/${t}/introns_binned.tsv | sed 's/ /\t/g' > ${DIR}/PipelineCalibration/${species}-introns_binned_${t}.tmp

	echo "INPUT:" `awk 'END {print NR}' ${INPUT}`
	echo "binned:" `awk 'END {print NR}' ${DIR}/PipelineCalibration/${species}-introns_binned_${t}.tmp`

	cut -f 2 ${INPUT} | sort | uniq > ${DIR}/PipelineCalibration/${species}-datasets_${t}.tmp

	numDatasets=`awk 'END {print NR}' ${DIR}/PipelineCalibration/${species}-datasets_${t}.tmp`

	for j in `seq 1 ${numDatasets}`; do

		dataset=`head -n ${j} ${DIR}/PipelineCalibration/${species}-datasets_${t}.tmp| tail -n 1`

		awk -F"\t" -v OFS="\t" -v name=${dataset} '{if ($2==name) print $1,$NF}' ${DIR}/PipelineCalibration/${species}-introns_binned_${t}.tmp | sort -k1,1 | uniq > ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp

		echo "${dataset}:" `awk 'END {print NR}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp`

		minor=`awk -F"\t" '{if ($NF=="minor") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		minor_hybrid=`awk -F"\t" '{if ($NF=="minor_hybrid") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		minor_like=`awk -F"\t" '{if ($NF=="minor-like") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		major=`awk -F"\t" '{if ($NF=="major") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		major_hybrid=`awk -F"\t" '{if ($NF=="major_hybrid") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		major_like=`awk -F"\t" '{if ($NF=="major-like") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
		non_canonical=`awk -F"\t" '{if ($NF=="non-canonical") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`

		echo ${species} ${dataset} ${t} ${minor} ${minor_hybrid} ${minor_like} ${major} ${major_hybrid} ${major_like} ${non_canonical} | sed 's/ /\t/g' >> ${DIR}/PipelineCalibration/ResponsiveIntronCounts-${t}.tsv

		TP=${minor}
		FP=`expr ${all_minor} - ${minor}`
		FN=`awk -v a=${minor_hybrid} -v b=${minor_like} -v c=${major} -v d=${major_hybrid} -v e=${major_like} -v f=${non_canonical} 'BEGIN {print a+b+c+d+e+f}'`

		precision=`awk -v tp=${TP} -v fp=${FP} 'BEGIN {if ((tp+fp)!=0) print tp/(tp+fp); else print 0}'`
		recall=`awk -v tp=${TP} -v fn=${FN} 'BEGIN {if ((tp+fn)!=0) print tp/(tp+fn); else print 0}'`

		F1=`awk -v p=${precision} -v r=${recall} 'BEGIN {if ((p+r)!=0) print (2*p*r)/(p+r); else print 0}'`
		F2=`awk -v p=${precision} -v r=${recall} 'BEGIN {if (((1+2**2)*p+r)!=0) print (1+2**2)*p*r/((1+2**2)*p+r); else print 0}'`

		echo ${species} ${dataset} ${t} ${TP} ${FP} ${FN} ${precision} ${recall} ${F1} ${F2}

		echo ${species} ${dataset} ${t} ${TP} ${FP} ${FN} ${precision} ${recall} ${F1} ${F2} | sed 's/ /\t/g' >> ${DIR}/PipelineCalibration/PipelineCalibration-${t}.tsv

	done

	dataset="RESPONSIVE"

	awk -F"\t" -v OFS="\t" '{print $1,$NF}' ${DIR}/PipelineCalibration/${species}-introns_binned_${t}.tmp | sort -k1,1 | uniq > ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp

	minor=`awk -F"\t" '{if ($NF=="minor") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	minor_hybrid=`awk -F"\t" '{if ($NF=="minor_hybrid") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	minor_like=`awk -F"\t" '{if ($NF=="minor-like") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	major=`awk -F"\t" '{if ($NF=="major") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	major_hybrid=`awk -F"\t" '{if ($NF=="major_hybrid") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	major_like=`awk -F"\t" '{if ($NF=="major-like") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`
	non_canonical=`awk -F"\t" '{if ($NF=="non-canonical") print $0}' ${DIR}/PipelineCalibration/${species}-${dataset}_binned_${t}.tmp | wc -l`

	echo ${species} ${dataset} ${t} ${minor} ${minor_hybrid} ${minor_like} ${major} ${major_hybrid} ${major_like} ${non_canonical} | sed 's/ /\t/g' >> ${DIR}/PipelineCalibration/ResponsiveIntronCounts-${t}.tsv

	TP=${minor}
	FP=`expr ${all_minor} - ${minor}`
	FN=`awk -v a=${minor_hybrid} -v b=${minor_like} -v c=${major} -v d=${major_hybrid} -v e=${major_like} -v f=${non_canonical} 'BEGIN {print a+b+c+d+e+f}'`

	precision=`awk -v tp=${TP} -v fp=${FP} 'BEGIN {if ((tp+fp)!=0) print tp/(tp+fp); else print 0}'`
	recall=`awk -v tp=${TP} -v fn=${FN} 'BEGIN {if ((tp+fn)!=0) print tp/(tp+fn); else print 0}'`

	F1=`awk -v p=${precision} -v r=${recall} 'BEGIN {if ((p+r)!=0) print (2*p*r)/(p+r); else print 0}'`
	F2=`awk -v p=${precision} -v r=${recall} 'BEGIN {if (((1+2**2)*p+r)!=0) print (1+2**2)*p*r/((1+2**2)*p+r); else print 0}'`

	echo ${species} ${dataset} ${t} ${TP} ${FP} ${FN} ${precision} ${recall} ${F1} ${F2}

	echo ${species} ${dataset} ${t} ${TP} ${FP} ${FN} ${precision} ${recall} ${F1} ${F2} | sed 's/ /\t/g' >> ${DIR}/PipelineCalibration/PipelineCalibration-${t}.tsv

	echo ""

done

rm -f ${DIR}/PipelineCalibration/*_${t}.tmp