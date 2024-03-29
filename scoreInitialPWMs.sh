#!/bin/bash

# Set inputs: working directory, list of species in source

DIR=$1		# Working directory

# Set species name, assembly, source

species=$2
assembly=$3
source=$4

# Set input file

IntronFASTA=${DIR}/IntronFASTAs/${species}.${assembly}.${source}-Introns.fa

# Extract major 5'ss nucleotide positions from -2 to +6 for scoring

PWM=${DIR}/PWMs/Initial/263_Species-U2_DonorSite_LOD.tsv

awk -F"\t" -v OFS="\t" '{print (substr($2,2,1)),(substr($2,3,1)),(substr($2,4,1)),(substr($2,5,1)),(substr($2,6,1)),(substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),$1}' ${IntronFASTA} > ${DIR}/output/${species}/U2_DonorSite.tmp

# Assign LOD for each nucleotide at each position

for Pos in {1..8}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

# Score each position, sum LOD score

awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7; X8=$8}

{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}
{if (X8=="A") Y8='${A8}'; else if (X8=="C") Y8='${C8}'; else if (X8=="G") Y8='${G8}'; else if (X8=="T") Y8='${T8}'}
		
{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8)/8; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2)/8} {print $NF,$1$2$3$4$5$6$7$8,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8,mean,var}' ${DIR}/output/${species}/U2_DonorSite.tmp > ${DIR}/output/${species}/U2_DonorSite_scored.tmp

# Scale scores

max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${DIR}/output/${species}/U2_DonorSite_scored.tmp | sort -k 1,1 > ${DIR}/output/${species}/Initial/U2_DonorSite_scaled.tsv

# Extract minor 5'ss nucleotide positions from +4 to +9 for scoring

awk -F"\t" -v OFS="\t" '{print (substr($2,7,1)),(substr($2,8,1)),(substr($2,9,1)),(substr($2,10,1)),(substr($2,11,1)),(substr($2,12,1)),$1}' ${IntronFASTA} > ${DIR}/output/${species}/U12_DonorSite.tmp

PWM=${DIR}/PWMs/Initial/263_Species-U12_DonorSite_LOD.tsv

# Assign LOD for each nucleotide at each position

for Pos in {1..6}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

# Score each position, sum LOD score

awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6}

{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G2}'; else if (X2=="T") Y2='${T2}'}
{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
		
{mean=(Y1+Y2+Y3+Y4+Y5+Y6)/6; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2)/6} {print $NF,$1$2$3$4$5$6,Y1+Y2+Y3+Y4+Y5+Y6,mean,var}' ${DIR}/output/${species}/U12_DonorSite.tmp > ${DIR}/output/${species}/U12_DonorSite_scored.tmp

# Scale scores

max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${DIR}/output/${species}/U12_DonorSite_scored.tmp | sort -k 1,1 > ${DIR}/output/${species}/Initial/U12_DonorSite_scaled.tsv

# Score acceptor site using LOD matrix

PWM=${DIR}/PWMs/Initial/263_Species-U2_AcceptorSite_LOD.tsv

# Assign LOD for each nucleotide at each position

for Pos in {1..13}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

# Extract acceptor site nucleotide positions from -13 to -1 at 3'ss for scoring

awk -v OFS="\t" '{print (substr($3,1,1)),(substr($3,2,1)),(substr($3,3,1)),(substr($3,4,1)),(substr($3,5,1)),(substr($3,6,1)),(substr($3,7,1)),(substr($3,8,1)),(substr($3,9,1)),(substr($3,10,1)),(substr($3,11,1)),(substr($3,12,1)),(substr($3,13,1)),$1}' ${IntronFASTA} > ${DIR}/output/${species}/U2_AcceptorSite.tmp

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

{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13)/13; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2+(Y13-mean)**2)/13} {print $NF,$1$2$3$4$5$6$7$8$9$10$11$12$13,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13,mean,var}' ${DIR}/output/${species}/U2_AcceptorSite.tmp > ${DIR}/output/${species}/U2_AcceptorSite_scored.tmp

# Scale scores

max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

# Scale, calculate mean, calculate variance

awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${DIR}/output/${species}/U2_AcceptorSite_scored.tmp | sort -k 1,1 > ${DIR}/output/${species}/Initial/U2_AcceptorSite_scaled.tsv

# Score U2 BPS using LOD matrix

PWM=${DIR}/PWMs/Initial/263_Species-U2_BranchPoint_LOD.tsv

# Assign LOD for each nucleotide at each position

for Pos in {1..7}; do
for N in A C G T; do

	eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

done 
done

# Get sliding window from -44 to -4 at 3'ss

awk -F"\t" -v OFS="\t" '{print (substr($5,1,1)),(substr($5,2,1)),(substr($5,3,1)),(substr($5,4,1)),(substr($5,5,1)),(substr($5,6,1)),(substr($5,7,1)),(substr($5,8,1)),(substr($5,9,1)),(substr($5,10,1)),(substr($5,11,1)),(substr($5,12,1)),(substr($5,13,1)),(substr($5,14,1)),(substr($5,15,1)),(substr($5,16,1)),(substr($5,17,1)),(substr($5,18,1)),(substr($5,19,1)),(substr($5,20,1)),(substr($5,21,1)),(substr($5,22,1)),(substr($5,23,1)),(substr($5,24,1)),(substr($5,25,1)),(substr($5,26,1)),(substr($5,27,1)),$1}' ${IntronFASTA} > ${DIR}/output/${species}/U2_BranchPoint.tmp

# Generate all potential 7 nt BPS from stretch

awk -F"\t" -v OFS="\t" '{n=NF-6} {for (i=1;i<n;i++) {print $i,$(i+1),$(i+2),$(i+3),$(i+4),$(i+5),$(i+6),$NF}}' ${DIR}/output/${species}/U2_BranchPoint.tmp > ${DIR}/output/${species}/U2_PotentialBPS.tmp

# Accumulate each intron score from LOD score

awk -F"\t" -v OFS="\t" '{X1=$1; X2=$2; X3=$3; X4=$4; X5=$5; X6=$6; X7=$7}

{if (X1=="A") Y1='${A1}'; else if (X1=="C") Y1='${C1}'; else if (X1=="G") Y1='${G1}'; else if (X1=="T") Y1='${T1}'}
{if (X2=="A") Y2='${A2}'; else if (X2=="C") Y2='${C2}'; else if (X2=="G") Y2='${G1}'; else if (X2=="T") Y2='${T2}'}
{if (X3=="A") Y3='${A3}'; else if (X3=="C") Y3='${C3}'; else if (X3=="G") Y3='${G3}'; else if (X3=="T") Y3='${T3}'}
{if (X4=="A") Y4='${A4}'; else if (X4=="C") Y4='${C4}'; else if (X4=="G") Y4='${G4}'; else if (X4=="T") Y4='${T4}'}
{if (X5=="A") Y5='${A5}'; else if (X5=="C") Y5='${C5}'; else if (X5=="G") Y5='${G5}'; else if (X5=="T") Y5='${T5}'}
{if (X6=="A") Y6='${A6}'; else if (X6=="C") Y6='${C6}'; else if (X6=="G") Y6='${G6}'; else if (X6=="T") Y6='${T6}'}
{if (X7=="A") Y7='${A7}'; else if (X7=="C") Y7='${C7}'; else if (X7=="G") Y7='${G7}'; else if (X7=="T") Y7='${T7}'}

{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7)/7; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2)/7} {print $NF,X1,X2,X3,X4,X5,X6,X7"::"mean"::"var,Y1+Y2+Y3+Y4+Y5+Y6+Y7}' ${DIR}/output/${species}/U2_PotentialBPS.tmp > ${DIR}/output/${species}/U2_PotentialBPS_scored.tmp

# Extract the highest scoring BPS from each transcript

${DIR}/bin/key-merge ${DIR}/output/${species}/U2_PotentialBPS_scored.tmp > ${DIR}/output/${species}/U2_PotentialBPS_merged.tmp
awk -v OFS="\t" '{max=$9; n=(NF-1)/8; for (i=1;i<=n;i++) {j=i*8+1; if (max<=$j) {max=$j; max_i=j}}; print $(max_i-7),$(max_i-6),$(max_i-5),$(max_i-4),$(max_i-3),$(max_i-2),$(max_i-1),$1,max}' ${DIR}/output/${species}/U2_PotentialBPS_merged.tmp | awk -F"\t" -v OFS="\t" '{print $8,$1$2$3$4$5$6$7,$9}' | awk -F"\t|::" -v OFS="\t" '{print $1,$2,$5,$3,$4}' > ${DIR}/output/${species}/U2_HighestScoringBPS_scored.tmp
	
# Determine max and min scores

max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

# Scale scores

awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${DIR}/output/${species}/U2_HighestScoringBPS_scored.tmp | sort -k 1,1 > ${DIR}/output/${species}/Initial/U2_BranchPoint_scaled.tsv

# Extract potential U12 BPS nucleotides from -40 to -1

awk -F"\t" -v OFS="\t" '{print (substr($4,1,1)),(substr($4,2,1)),(substr($4,3,1)),(substr($4,4,1)),(substr($4,5,1)),(substr($4,6,1)),(substr($4,7,1)),(substr($4,8,1)),(substr($4,9,1)),(substr($4,10,1)),(substr($4,11,1)),(substr($4,12,1)),(substr($4,13,1)),(substr($4,14,1)),(substr($4,15,1)),(substr($4,16,1)),(substr($4,17,1)),(substr($4,18,1)),(substr($4,19,1)),(substr($4,20,1)),(substr($4,21,1)),(substr($4,22,1)),(substr($4,23,1)),(substr($4,24,1)),(substr($4,25,1)),(substr($4,26,1)),(substr($4,27,1)),(substr($4,28,1)),(substr($4,29,1)),(substr($4,30,1)),(substr($4,31,1)),(substr($4,32,1)),(substr($4,33,1)),(substr($4,34,1)),(substr($4,35,1)),(substr($4,36,1)),(substr($4,37,1)),(substr($4,38,1)),(substr($4,39,1)),(substr($4,40,1)),$1}' ${IntronFASTA} > ${DIR}/output/${species}/U12_BranchPoint.tmp

# Generate potential 12nt BPS using a sliding window

awk -F"\t" -v OFS="\t" '{n=NF-11} {for (i=1;i<n;i++) {print $i,$(i+1),$(i+2),$(i+3),$(i+4),$(i+5),$(i+6),$(i+7),$(i+8),$(i+9),$(i+10),$(i+11),$NF}}' ${DIR}/output/${species}/U12_BranchPoint.tmp > ${DIR}/output/${species}/U12_PotentialBPS.tmp

for nt in A9 A10; do

	PWM=${DIR}/PWMs/Initial/263_Species-U12_BranchPoint_${nt}_LOD.tsv

	# Assign LOD for each nucleotide at each position

	for Pos in {1..12}; do
	for N in A C G T; do

		eval "${N}${Pos}"="`cat ${PWM} | head -n ${Pos} | tail -n 1 | awk -F"\t" -v n=${N} '{if (n=="A") print $1; else if (n=="C") print $2; else if (n=="G") print $3; else if (n=="T") print $4}'`"

	done 
	done

	# Score each position, sum LOD score

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

	{mean=(Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12)/12; var=((Y1-mean)**2+(Y2-mean)**2+(Y3-mean)**2+(Y4-mean)**2+(Y5-mean)**2+(Y6-mean)**2+(Y7-mean)**2+(Y8-mean)**2+(Y9-mean)**2+(Y10-mean)**2+(Y11-mean)**2+(Y12-mean)**2)/12} {print $NF,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12"::"mean"::"var,Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12}' ${DIR}/output/${species}/U12_PotentialBPS.tmp > ${DIR}/output/${species}/U12_PotentialBPS_${nt}_scored.tmp

	# Extract the highest scoring BPS & score from each transcript

	${DIR}/bin/key-merge ${DIR}/output/${species}/U12_PotentialBPS_${nt}_scored.tmp > ${DIR}/output/${species}/U12_PotentialBPS_${nt}_merged.tmp
	awk -v OFS="\t" '{max=$14; n=(NF-1)/13; for (i=1;i<=n;i++) {j=i*13+1; if (max<=$j) {max=$j; max_i=j}}; print $(max_i-12),$(max_i-11),$(max_i-10),$(max_i-9),$(max_i-8),$(max_i-7),$(max_i-6),$(max_i-5),$(max_i-4),$(max_i-3),$(max_i-2),$(max_i-1),$1,max}' ${DIR}/output/${species}/U12_PotentialBPS_${nt}_merged.tmp | awk -F"\t" -v OFS="\t" '{print $13,$1$2$3$4$5$6$7$8$9$10$11$12,$14}' | awk -F"\t|::" -v OFS="\t" '{print $1,$2,$5,$3,$4}' > ${DIR}/output/${species}/U12_HighestScoringBPS_${nt}_scored.tmp

	# Scale scores

	max=`awk -F"\t" '{max=0; for (i=1;i<=NF;i++) {if (max<$i) max=$i}; print max}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`
	min=`awk -F"\t" '{min=0; for (i=1;i<=NF;i++) {if (min>$i) min=$i}; print min}' ${PWM} | awk 'BEGIN {sum=0} {sum+=$1} END {print sum}'`

	awk -F"\t" -v OFS="\t" -v min=${min} -v max=${max} '{sum=$(NF-2); mean=$(NF-1); var=$NF} {if (sum<0) score=(sum-min)/(0-min)*50; else if (sum>0) score=(sum/max)*50+50; else score=50} {print $1,$2,score,mean,var}' ${DIR}/output/${species}/U12_HighestScoringBPS_${nt}_scored.tmp | sort -k 1,1 > ${DIR}/output/${species}/Initial/U12_BranchPoint_${nt}_scaled.tsv

done

join ${DIR}/output/${species}/Initial/U12_BranchPoint_A9_scaled.tsv ${DIR}/output/${species}/Initial/U12_BranchPoint_A10_scaled.tsv | awk -v OFS="\t" '{if ($3>=$7) print $1,$2,$3,$4,$5; else print $1,$6,$7,$8,$9}' | sort -k 1,1 > ${DIR}/output/${species}/Initial/U12_BranchPoint_scaled.tsv

for type in U2_DonorSite U12_DonorSite U2_AcceptorSite U2_BranchPoint U12_BranchPoint; do

	if [ -s ${DIR}/output/${species}/Initial/introns_scored.tsv ]; then

		awk -F"\t" -v OFS="\t" '{print $2,$3}' ${DIR}/output/${species}/Initial/${type}_scaled.tsv > ${DIR}/output/${species}/Initial/${type}_scaled.tmp
		paste ${DIR}/output/${species}/Initial/introns_scored.tsv ${DIR}/output/${species}/Initial/${type}_scaled.tmp > ${DIR}/output/${species}/Initial/introns_scored.tmp
		mv ${DIR}/output/${species}/Initial/introns_scored.tmp ${DIR}/output/${species}/Initial/introns_scored.tsv

	else

		awk -F"\t" -v OFS="\t" '{print $1,$2,$3}' ${DIR}/output/${species}/Initial/${type}_scaled.tsv > ${DIR}/output/${species}/Initial/introns_scored.tsv

	fi

done

rm -f ${DIR}/output/${species}/*.tmp
rm -f ${DIR}/output/${species}/Initial/*.tmp