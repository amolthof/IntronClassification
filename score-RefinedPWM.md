***Set inputs: working workdirectory, list of species in source***

	workdir=$1
	SpeciesList=$2
	species=$3
	assembly=$4
	version=$5

	INPUT=${workdir}/output/${species}/Initial/introns_scored.tsv
	IntronFASTA=${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa

	# FIELDS: intron_id, 5'ss, U2 5'ss, U2 5'ss score, U12 5'ss, U12 5'ss score, U2 3'ss, U2 3'ss score, U2 BPS, U2 BPS score, U12 BPS, U12 BPS score

	sort -k 1,1 ${INPUT} > ${workdir}/output/${species}/introns_splice_sites.tmp
	join ${IntronFASTA} ${workdir}/output/${species}/introns_splice_sites.tmp | sed 's/ /\t/g' > ${workdir}/output/${species}/introns.tmp

	INPUT=${workdir}/output/${species}/introns.tmp

	numSpecies=`awk 'END {print NR}' ${SpeciesList}`

***Extract nucleotide positions for scoring***

	# Extract major 5'ss nucleotide positions from -2 to +6 for scoring
	PWM=${workdir}/PWMs/Initial/${numSpecies}_Species-U2_DonorSite_LOD.tsv
	awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),$1}' ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa > ${workdir}/output/${species}/U2_DonorSite.tmp

	# Extract minor 5'ss nucleotide positions from +4 to +9 for scoring
	awk -F"\t" -v OFS="\t" '{print (substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1)),$1}' ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa > ${workdir}/output/${species}/U12_DonorSite.tmp

	# Extract 5'ss nucleotide positions from -2 to +9 for scoring
	awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1)),$1}' ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_5SS.fa > ${workdir}/output/${species}/DonorSite.tmp

	# Extract acceptor site nucleotide positions from -13 to -1 at 3'ss for scoring
	awk -F"\t" -v OFS="\t" '{print (substr($2,1,1)),(substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1)),(substr($2,13,1)),$1}' ${workdir}/output/${species}/Initial/U2_AcceptorSite_scaled.tsv > ${workdir}/output/${species}/U2_AcceptorSite.tmp

	# Extract 7nt U2 BPS for scoring
	awk -F"\t" -v OFS="\t" '{print (substr($2,1,1)),(substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),$1}' ${workdir}/output/${species}/Initial/U2_BranchPoint_scaled.tsv > ${workdir}/output/${species}/U2_BranchPoint.tmp

	# Extract 12nt U12 BPS for scoring
	awk -F"\t" -v OFS="\t" '{print (substr($2,1,1)),(substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1)),$1}' ${workdir}/output/${species}/Initial/U12_BranchPoint_scaled.tsv > ${workdir}/output/${species}/U12_BranchPoint.tmp

***Score U2 donor site***

	# Score U2_DonorSite from -2 to +6 for scoring (8 nt)

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_DonorSite_LOD.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..8}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}
	{if (X8=="A") Y8='${A8}'; else if (X8=="C") Y8='${C8}'; else if (X8=="G") Y8='${G8}'; else if (X8=="T") Y8='${T8}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)/8; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2)/8} {print $NF,$1$2$3$4$5$6$7$8,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8,mean,var}' ${workdir}/output/${species}/U2_DonorSite.tmp > ${workdir}/output/${species}/U2_DonorSite_scored.tmp 

	# Scale scores

	max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
	min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

	# Scale, calculate mean, calculate variance

	awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/U2_DonorSite_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_DonorSite_scaled.tsv

	# Score using scaled matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_DonorSite_scaled.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..8}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	# Score each position, sum LOD score, get mean, calculate variance

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}
	{if (X8=="A") Y8='${A8}'; else if (X8=="C") Y8='${C8}'; else if (X8=="G") Y8='${G8}'; else if (X8=="T") Y8='${T8}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)/8; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2)/8} {print $NF,$1$2$3$4$5$6$7$8,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8,mean,var}' ${workdir}/output/${species}/U2_DonorSite.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_DonorSite_match.tsv

***Score U12 donor site***

	# Score U12_DonorSite from +4 to +9 for scoring (6 nt)

	PWM=${workdir}/PWMs/${numSpecies}_Species-U12_DonorSite_LOD.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..6}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6)/6; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2)/6} {print $NF,$1$2$3$4$5$6,Y1+Y2+Y3+Y4+Y5+Y6,mean,var}' ${workdir}/output/${species}/U12_DonorSite.tmp > ${workdir}/output/${species}/U12_DonorSite_scored.tmp 

	# Scale scores

	max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
	min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

	# Scale, calculate mean, calculate variance

	awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/U12_DonorSite_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/U12_DonorSite_scaled.tsv

	# Score using scaled matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U12_DonorSite_scaled.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..6}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	# Score each position, sum LOD score, get mean, calculate variance

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6)/6; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2)/6} {print $NF,$1$2$3$4$5$6,Y1+Y2+Y3+Y4+Y5+Y6,mean,var}' ${workdir}/output/${species}/U12_DonorSite.tmp | sort -k 1,1 > ${workdir}/output/${species}/U12_DonorSite_match.tsv

***Score full donor site***
	# Score full donor site

	for type in U2_GT-AG_DonorSite U2_GC-AG_DonorSite U12_GT-AG_DonorSite U12_AT-AC_DonorSite; do

		# Score donor site using LOD matrix

		PWM=${workdir}/PWMs/${numSpecies}_Species-${type}_LOD.tsv

		# Assign LOD for each nucleotide at each position

		for Pos in {1..11}; do
		for N in A C G T; do

			eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

		done 
		done

		awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11}

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

		{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11)/11; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2)/11} {print $NF,$1$2$3$4$5$6$7$8$9$10$11,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11,mean,var}' ${workdir}/output/${species}/DonorSite.tmp > ${workdir}/output/${species}/${type}_scored.tmp 

		# Scale scores

		max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
		min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

		# Scale, calculate mean, calculate variance

		awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/${type}_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/${type}_scaled.tsv

		# Score using scaled matrix

		PWM=${workdir}/PWMs/${numSpecies}_Species-${type}_scaled.tsv

		# Assign LOD for each nucleotide at each position

		for Pos in {1..11}; do
		for N in A C G T; do

			eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

		done 
		done

		# Score each position, sum LOD score, get mean, calculate variance

		awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11}

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

		{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11)/11; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2)/11} {print $NF,$1$2$3$4$5$6$7$8$9$10$11,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11,mean,var}' ${workdir}/output/${species}/DonorSite.tmp | sort -k 1,1 > ${workdir}/output/${species}/${type}_match.tsv

	done

***Select highest scaled donor site***

	# Pick highest scaled donor site score for major 

	mv ${workdir}/output/${species}/U2_DonorSite_scaled.tsv ${workdir}/output/${species}/U2_DonorSite_8nt_scaled.tsv
	mv ${workdir}/output/${species}/U2_DonorSite_match.tsv ${workdir}/output/${species}/U2_DonorSite_8nt_match.tsv

	join ${workdir}/output/${species}/U2_GT-AG_DonorSite_scaled.tsv ${workdir}/output/${species}/U2_GC-AG_DonorSite_scaled.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' | sort -k1,1 > ${workdir}/output/${species}/U2_DonorSite_11nt_scaled.tsv
	join ${workdir}/output/${species}/U2_GT-AG_DonorSite_match.tsv ${workdir}/output/${species}/U2_GC-AG_DonorSite_match.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' | sort -k1,1 > ${workdir}/output/${species}/U2_DonorSite_11nt_match.tsv

	join ${workdir}/output/${species}/U2_DonorSite_11nt_scaled.tsv ${workdir}/output/${species}/U2_DonorSite_8nt_scaled.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' | sort -k1,1 > ${workdir}/output/${species}/U2_DonorSite_scaled.tsv
	join ${workdir}/output/${species}/U2_DonorSite_11nt_match.tsv ${workdir}/output/${species}/U2_DonorSite_8nt_match.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' | sort -k1,1 > ${workdir}/output/${species}/U2_DonorSite_match.tsv

	# Pick highest scaled donor site score for minor

	mv ${workdir}/output/${species}/U12_DonorSite_scaled.tsv ${workdir}/output/${species}/U12_DonorSite_6nt_scaled.tsv
	mv ${workdir}/output/${species}/U12_DonorSite_match.tsv ${workdir}/output/${species}/U12_DonorSite_6nt_match.tsv

	join ${workdir}/output/${species}/U12_GT-AG_DonorSite_scaled.tsv ${workdir}/output/${species}/U12_AT-AC_DonorSite_scaled.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' > ${workdir}/output/${species}/U12_DonorSite_11nt_scaled.tsv
	join ${workdir}/output/${species}/U12_GT-AG_DonorSite_match.tsv ${workdir}/output/${species}/U12_AT-AC_DonorSite_match.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' > ${workdir}/output/${species}/U12_DonorSite_11nt_match.tsv

	join ${workdir}/output/${species}/U12_DonorSite_11nt_scaled.tsv ${workdir}/output/${species}/U12_DonorSite_6nt_scaled.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' > ${workdir}/output/${species}/U12_DonorSite_scaled.tsv
	join ${workdir}/output/${species}/U12_DonorSite_11nt_match.tsv ${workdir}/output/${species}/U12_DonorSite_6nt_match.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9,$10}' > ${workdir}/output/${species}/U12_DonorSite_match.tsv

***Score U2 acceptor site***

	# Score acceptor site using LOD matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_AcceptorSite_LOD.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..13}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	# Score each position, sum LOD score, get mean, calculate variance

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11; X12=$12; X13=$13}

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
	{if (X13=="A") Y13='${A13}'; else if (X13=="C") Y13='${C13}'; else if (X13=="G") Y13='${G13}'; else if (X13=="T") Y13='${T13}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13)/13; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2+(Y13-mean)**2)/13} {print $NF,$1$2$3$4$5$6$7$8$9$10$11$12$13,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13,mean,var}' ${workdir}/output/${species}/U2_AcceptorSite.tmp > ${workdir}/output/${species}/U2_AcceptorSite_scored.tmp

	# Scale scores

	max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
	min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

	# Scale, calculate mean, calculate variance

	awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/U2_AcceptorSite_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_AcceptorSite_scaled.tsv

	# Score using scaled matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_AcceptorSite_scaled.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..13}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	# Score each position, sum LOD score, get mean, calculate variance

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11; X12=$12; X13=$13}

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
	{if (X13=="A") Y13='${A13}'; else if (X13=="C") Y13='${C13}'; else if (X13=="G") Y13='${G13}'; else if (X13=="T") Y13='${T13}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13)/13; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2+(Y13-mean)**2)/13} {print $NF,$1$2$3$4$5$6$7$8$9$10$11$12$13,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13,mean,var}' ${workdir}/output/${species}/U2_AcceptorSite.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_AcceptorSite_match.tsv

***Score U2 BPS***

	# Score U2 BPS using LOD matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_BranchPoint_LOD.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..7}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done
	done

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7)/7; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2)/7} {print $NF,$1$2$3$4$5$6$7,Y1+Y2+Y3+Y4+Y5+Y6+Y7,mean,var}' ${workdir}/output/${species}/U2_BranchPoint.tmp > ${workdir}/output/${species}/U2_BranchPoint_scored.tmp

	# Scale scores

	max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
	min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

	# Scale, calculate mean, calculate variance

	awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/U2_BranchPoint_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_BranchPoint_scaled.tsv

	# Score U2 BPS using scaled matrix

	PWM=${workdir}/PWMs/${numSpecies}_Species-U2_BranchPoint_scaled.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..7}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done
	done

	awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7}

	{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
	{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
	{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
	{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
	{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
	{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
	{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7)/7; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2)/7} {print $NF,$1$2$3$4$5$6$7,Y1+Y2+Y3+Y4+Y5+Y6+Y7,mean,var}' ${workdir}/output/${species}/U2_BranchPoint.tmp | sort -k 1,1 > ${workdir}/output/${species}/U2_BranchPoint_match.tsv

***Score U12 BPS***

	for nt in A9 A10; do

		# Score BPS using LOD matrix

		PWM=${workdir}/PWMs/${numSpecies}_Species-U12_BranchPoint_${nt}_LOD.tsv

		# Assign LOD for each nucleotide at each position

		for Pos in {1..12}; do
		for N in A C G T; do

			eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

		done 
		done

		# Score each position, sum LOD score, get mean, calculate variance

		awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11; X12=$12}

		{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
		{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G1}'; else if (X2=="T") Y2='${T2}'}
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

		{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y12)/12} {print $NF,$1$2$3$4$5$6$7$8$9$10$11$12,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12,mean,((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2)/12}' ${workdir}/output/${species}/U12_BranchPoint.tmp > ${workdir}/output/${species}/U12_BranchPoint_${nt}_scored.tmp

		# Scale scores

		max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
		min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

		# Scale, calculate mean, calculate variance

		awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${workdir}/output/${species}/U12_BranchPoint_${nt}_scored.tmp | sort -k 1,1 > ${workdir}/output/${species}/U12_BranchPoint_${nt}_scaled.tsv

		# Score using scaled matrix

		PWM=${workdir}/PWMs/${numSpecies}_Species-U12_BranchPoint_${nt}_scaled.tsv

		# Assign LOD for each nucleotide at each position

		for Pos in {1..12}; do
		for N in A C G T; do

			eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

		done 
		done

		# Score each position, sum LOD score, get mean, calculate variance

		awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8; X9=$9; X10=$10; X11=$11; X12=$12}

		{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
		{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G1}'; else if (X2=="T") Y2='${T2}'}
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

		{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y12)/12} {print $NF,$1$2$3$4$5$6$7$8$9$10$11$12,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12,mean,((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2)/12}' ${workdir}/output/${species}/U12_BranchPoint.tmp | sort -k1,1 > ${workdir}/output/${species}/U12_BranchPoint_${nt}_match.tsv

	done

	for file in scaled match; do

		join ${workdir}/output/${species}/U12_BranchPoint_A9_${file}.tsv ${workdir}/output/${species}/U12_BranchPoint_A10_${file}.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,"A9"; else print $1,"A10"}' > ${workdir}/output/${species}/U12_BranchPoint_pick.tmp

		for type in A9 A10; do

			grep -F -w ${type} ${workdir}/output/${species}/U12_BranchPoint_pick.tmp | cut -f 1 > ${workdir}/output/${species}/U12_BranchPoint_${type}.tmp
			join ${workdir}/output/${species}/U12_BranchPoint_${type}.tmp ${workdir}/output/${species}/U12_BranchPoint_${type}_${file}.tsv | sed 's/ /\t/g' >> ${workdir}/output/${species}/U12_BranchPoint_${file}.tmp

		done

		sort -k1,1 ${workdir}/output/${species}/U12_BranchPoint_${file}.tmp > ${workdir}/output/${species}/U12_BranchPoint_${file}.tsv

	done

***Combine splice site scores to one file***
	# Merge scores

	for type in U2_DonorSite_11nt U12_DonorSite_11nt U2_AcceptorSite U2_BranchPoint U12_BranchPoint; do

		if [ -s ${workdir}/output/${species}/introns_scored.tsv ]; then

			awk -F"\t" -v OFS="\t" '{print $2,$3}' ${workdir}/output/${species}/${type}_scaled.tsv > ${workdir}/output/${species}/${type}_scaled.tmp
			paste ${workdir}/output/${species}/introns_scored.tsv ${workdir}/output/${species}/${type}_scaled.tmp > ${workdir}/output/${species}/introns_scored.tmp
			mv ${workdir}/output/${species}/introns_scored.tmp ${workdir}/output/${species}/introns_scored.tsv

		else

			awk -F"\t" -v OFS="\t" '{print $1,$2,$3}' ${workdir}/output/${species}/${type}_scaled.tsv > ${workdir}/output/${species}/introns_scored.tsv

		fi

	done

	sort -k 1,1 ${workdir}/output/${species}/introns_scored.tsv > ${workdir}/output/${species}/introns_scored.tmp
	mv ${workdir}/output/${species}/introns_scored.tmp ${workdir}/output/${species}/introns_scored.tsv

	rm -f ${workdir}/output/${species}/*.tmp
