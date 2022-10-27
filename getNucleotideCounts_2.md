***Set inputs: working directory, species***

	workdir=$1
	species=$2
	assembly=$3
	version=$4
	X=22.5
	Y=75

	# Set input file & intron FASTA
	# FIELDS: intron_id, U2 5'ss, U2 5'ss score, U12 5'ss, U12 5'ss score, U2 3'ss, U2 3'ss score, U2 BPS, U2 BPS score, U12 BPS, U12 BPS score

	INPUT=${workdir}/output/${species}/Initial/introns_scored.tsv
	IntronFASTA=${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa

***Screen for putative minor and major introns***
	awk -F"\t" -v OFS="\t" -v x=${X} -v y=${Y} '{d=substr($2,3,2); a=substr($6,12,2); U2=$3; U12=$5; PPT=$7; BPS=$11} {if (U2>50 && U2>U12 && PPT>50) print $0,d"-"a,"U2"; else if (U12>50 && U12>(U2+x) && BPS>y) print $0,d"-"a,"U12"; else print $0,d"-"a,"discard"}' ${INPUT} | sort -k 1,1 > ${workdir}/output/${species}/introns_screened_1.tmp

	# Extract 5'ss sequence
	join ${IntronFASTA} ${workdir}/output/${species}/introns_screened_1.tmp | sed 's/ /\t/g' > ${workdir}/output/${species}/introns_screened_2.tmp

	# Group by preliminary class
	grep -F -w "U2" ${workdir}/output/${species}/introns_screened_2.tmp > ${workdir}/output/${species}/U2_introns.tmp
	grep -F -w "U12" ${workdir}/output/${species}/introns_screened_2.tmp > ${workdir}/output/${species}/U12_introns.tmp
	grep -F -w "discard" ${workdir}/output/${species}/introns_screened_2.tmp > ${workdir}/output/${species}/discard_introns.tmp

	# Group by terminal dinucleotide

	for type in U2 U12 discard; do
		grep -F -w "GT-AG" ${workdir}/output/${species}/${type}_introns.tmp > ${workdir}/output/${species}/${type}_GT-AG_introns.tmp
		grep -F -w "GC-AG" ${workdir}/output/${species}/${type}_introns.tmp > ${workdir}/output/${species}/${type}_GC-AG_introns.tmp
		grep -F -w "AT-AC" ${workdir}/output/${species}/${type}_introns.tmp > ${workdir}/output/${species}/${type}_AT-AC_introns.tmp
		grep -F -w -v "GT-AG" ${workdir}/output/${species}/${type}_introns.tmp | grep -F -w -v "GC-AG" | grep -F -w -v "AT-AC" > ${workdir}/output/${species}/${type}_other_introns.tmp
	done

***Obtain nucleotide counts for different splice sites***

	# PWM 1: U2_DonorSite PWM, -2 to +6, 8nt sequence from U2 introns
	awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1))}' ${workdir}/output/${species}/U2_introns.tmp > ${workdir}/output/${species}/U2_DonorSite.tmp

	for Pos in {1..8}; do
		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${workdir}/output/${species}/U2_DonorSite.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-U2_DonorSite.tsv
	done

	# PWM 2: U12_DonorSite PWM, +4 to +9, 6nt from U12 introns

	if [ -s ${workdir}/output/${species}/U12_introns.tmp ]; then
		awk -F"\t" -v OFS="\t" '{print (substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1))}' ${workdir}/output/${species}/U12_introns.tmp > ${workdir}/output/${species}/U12_DonorSite.tmp

		for Pos in {1..6}; do
			awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${workdir}/output/${species}/U12_DonorSite.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-U12_DonorSite.tsv
		done
	fi

	# PWMs 3-6: sub-type specific, full donor site, 11nt, -2 to +9
	
	for type in U2_GT-AG U2_GC-AG U12_GT-AG U12_AT-AC; do
		if [ ! -s ${workdir}/output/${species}/${type}_introns.tmp ]; then
			continue
		fi

		# Create position matrix of nucleotides at 5'ss from -2 to +9
		awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1))}' ${workdir}/output/${species}/${type}_introns.tmp > ${workdir}/output/${species}/${type}_DonorSite.tmp

		# Count nucleotides at each position
		for Pos in {1..11}; do
			awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
										{if (count["C"]>0) C=count["C"]; else C=0}
										{if (count["G"]>0) G=count["G"]; else G=0}
										{if (count["T"]>0) T=count["T"]; else T=0}
										{print A, C, G, T}}' ${workdir}/output/${species}/${type}_DonorSite.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-${type}_DonorSite.tsv
		done
	done 

	# PWM 7: major acceptor site, 13nt, -13 to -1
	awk -F"\t" -v OFS="\t" '{print (substr($7,1,1)),(substr($7,2,1)),(substr($7,3,1)),(substr($7,4,1)),(substr($7,5,1)),(substr($7,6,1)),(substr($7,7,1)),(substr($7,8,1)),(substr($7,9,1)),(substr($7,10,1)),(substr($7,11,1)),(substr($7,12,1)),(substr($7,13,1))}' ${workdir}/output/${species}/U2_introns.tmp > ${workdir}/output/${species}/U2_AcceptorSite.tmp

	for Pos in {1..13}; do
		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${workdir}/output/${species}/U2_AcceptorSite.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-U2_AcceptorSite.tsv
	done

	# PWM 8: major BPS, 7nt, force A6
	awk -F"\t" -v OFS="\t" '{print (substr($9,1,1)),(substr($9,2,1)),(substr($9,3,1)),(substr($9,4,1)),(substr($9,5,1)),(substr($9,6,1)),(substr($9,7,1))}' ${workdir}/output/${species}/U2_introns.tmp | awk -F"\t" '{if ($6=="A") print $0}' > ${workdir}/output/${species}/U2_BranchPoint_A6.tmp

	for Pos in {1..7}; do
		awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${workdir}/output/${species}/U2_BranchPoint_A6.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-U2_BranchPoint.tsv
	done

	# PWM 9-10: minor BPS, 12nt, force A9 & A10
	awk -F"\t" -v OFS="\t" '{print (substr($11,1,1)),(substr($11,2,1)),(substr($11,3,1)),(substr($11,4,1)),(substr($11,5,1)),(substr($11,6,1)),(substr($11,7,1)),(substr($11,8,1)),(substr($11,9,1)),(substr($11,10,1)),(substr($11,11,1)),(substr($11,12,1))}' ${workdir}/output/${species}/U12_introns.tmp | awk -F"\t" '{if ($9=="A") print $0}' > ${workdir}/output/${species}/U12_BranchPoint_A9.tmp
	awk -F"\t" -v OFS="\t" '{print (substr($11,1,1)),(substr($11,2,1)),(substr($11,3,1)),(substr($11,4,1)),(substr($11,5,1)),(substr($11,6,1)),(substr($11,7,1)),(substr($11,8,1)),(substr($11,9,1)),(substr($11,10,1)),(substr($11,11,1)),(substr($11,12,1))}' ${workdir}/output/${species}/U12_introns.tmp | awk -F"\t" '{if ($10=="A") print $0}' > ${workdir}/output/${species}/U12_BranchPoint_A10.tmp

	for nt in A9 A10; do
		for Pos in {1..12}; do
			awk -F"\t" -v OFS="\t" -v pos=${Pos} '{count[$(pos)]++} END {{if (count["A"]>0) A=count["A"]; else A=0}
									{if (count["C"]>0) C=count["C"]; else C=0}
									{if (count["G"]>0) G=count["G"]; else G=0}
									{if (count["T"]>0) T=count["T"]; else T=0}
									{print A, C, G, T}}' ${workdir}/output/${species}/U12_BranchPoint_${nt}.tmp >> ${workdir}/PWMs/NucleotideCounts/${species}-U12_BranchPoint_${nt}.tsv
		done
	done

***Output stats***

	total=`awk 'END {print NR}' ${INPUT}`
	major=`awk 'END {print NR}' ${workdir}/output/${species}/U2_introns.tmp`
	major_1=`awk 'END {print NR}' ${workdir}/output/${species}/U2_GT-AG_introns.tmp`
	major_2=`awk 'END {print NR}' ${workdir}/output/${species}/U2_GC-AG_introns.tmp`
	major_3=`awk 'END {print NR}' ${workdir}/output/${species}/U2_AT-AC_introns.tmp`
	major_4=`awk 'END {print NR}' ${workdir}/output/${species}/U2_other_introns.tmp`

	if [ -s ${workdir}/output/${species}/U12_introns.tmp ]; then
		minor=`awk 'END {print NR}' ${workdir}/output/${species}/U12_introns.tmp`
		minor_1=`awk 'END {print NR}' ${workdir}/output/${species}/U12_GT-AG_introns.tmp`
		minor_2=`awk 'END {print NR}' ${workdir}/output/${species}/U12_GC-AG_introns.tmp`
		minor_3=`awk 'END {print NR}' ${workdir}/output/${species}/U12_AT-AC_introns.tmp`
		minor_4=`awk 'END {print NR}' ${workdir}/output/${species}/U12_other_introns.tmp`
	else
		minor=0
		minor_1=0
		minor_2=0
		minor_3=0
		minor_4=0
	fi

	discard=`awk 'END {print NR}' ${workdir}/output/${species}/discard_introns.tmp`
	discard_1=`awk 'END {print NR}' ${workdir}/output/${species}/discard_GT-AG_introns.tmp`
	discard_2=`awk 'END {print NR}' ${workdir}/output/${species}/discard_GC-AG_introns.tmp`
	discard_3=`awk 'END {print NR}' ${workdir}/output/${species}/discard_AT-AC_introns.tmp`
	discard_4=`awk 'END {print NR}' ${workdir}/output/${species}/discard_other_introns.tmp`

	echo "species:" ${species} | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "assembly:" ${assembly} | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "version:" ${version} | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "" >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "${t}_PWM_GROUP" "COUNT" "CONDITION" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "total_introns:" ${total} "all" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "putative_major:" ${major} "U2>50,U2>U12,PPT>50" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GT-AG:" ${major_1}  "major,GT-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GC-AG:" ${major_2} "major,GC-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...AT-AC:" ${major_3} "major,AT-AC" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...other:" ${major_4} "major,other" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "putative_minor:" ${minor} "U12>50,U12>U2+${X},BPS>${Y}" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GT-AG:" ${minor_1} "minor,GT-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GC-AG:" ${minor_2} "minor,GC-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...AT-AC:" ${minor_3} "minor,AT-AC"| sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...other:" ${minor_4} "minor,other" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "discard:" ${discard} "other" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GT-AG:" ${discard_1} "discard,GT-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...GC-AG:" ${discard_2} "discard,GC-AG" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...AT-AC:" ${discard_3} "discard,AT-AC" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "...other:" ${discard_4} "discard,other" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
	echo "" >> ${workdir}/output/${species}/IntronClassificationStats.txt

	# Get consensus sequence
	
	for type in U2_DonorSite U2_GT-AG_DonorSite U2_GC-AG_DonorSite U12_DonorSite U12_GT-AG_DonorSite U12_AT-AC_DonorSite U2_AcceptorSite U2_BranchPoint U12_BranchPoint_A9 U12_BranchPoint_A10; do
			PWM=${workdir}/PWMs/NucleotideCounts/${species}-${type}.tsv
			
			if [ -s ${PWM} ]; then
				echo "${type}_PWM_count_consensus:" `${workdir}/bin/getConsensus.sh ${PWM}` | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
			else
				echo "${type}_PWM_count_consensus:" "NULL" | sed 's/ /\t/g' >> ${workdir}/output/${species}/IntronClassificationStats.txt
			fi
	done

	echo "" >> ${workdir}/output/${species}/IntronClassificationStats.txt

	rm -f ${workdir}/output/${species}/*.tmp
