#!/bin/bash

# Set inputs: working directory

DIR=$1		# Working directory

# Set species name, titration

species=$2

t=$3

# Set input

INPUT=${DIR}/output/${species}/${t}/introns_scored.tsv

# Separate major and minor introns from scored introns
# FIELDS: intron_id, 5'ss, U2 5'ss score, 5'ss, U12 5'ss score, 3'ss, U2 3'ss score, U2 BPS, U2 BPS score, U12 BPS, U12 BPS score

awk -F"\t" -v OFS="\t" '{d=substr($2,3,2); a=substr($6,12,2); U2=$3; U12=$5; PPT=$7; BPS=$11; if (U2>=80 && U2>=(U12+20)) print $0,d"-"a,"strong_major"; else if (U12>=80 && U12>=(U2+20)) print $0,d"-"a,"strong_minor"; else if (U2>=60 && U2>=(U12+20)) print $0,d"-"a,"major"; else if (U12>=60 && U12>=(U2+20)) print $0,d"-"a,"minor"; else if (U2>=50 && U2>=U12) print $0,d"-"a,"amb_major"; else if (U12>=50 && U12>=U2) print $0,d"-"a,"amb_minor"; else print $0,d"-"a,"non-canonical"}' ${DIR}/output/${species}/${t}/introns_scored.tsv > ${DIR}/output/${species}/${t}/introns_binned.tmp

awk '{if ($NF=="strong_major") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/strong_major_introns.tmp
awk '{if ($NF=="strong_minor") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/strong_minor_introns.tmp
awk '{if ($NF=="major") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/major_introns.tmp
awk '{if ($NF=="minor") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/minor_introns.tmp
awk '{if ($NF=="amb_major") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/amb_major_introns.tmp
awk '{if ($NF=="amb_minor") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/amb_minor_introns.tmp
awk '{if ($NF=="non-canonical") print $0}' ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/non-canonical_introns.tmp

cat ${DIR}/output/${species}/${t}/strong_major_introns.tmp ${DIR}/output/${species}/${t}/major_introns.tmp ${DIR}/output/${species}/${t}/amb_major_introns.tmp > ${DIR}/output/${species}/${t}/U2_introns.tmp
cat ${DIR}/output/${species}/${t}/strong_minor_introns.tmp ${DIR}/output/${species}/${t}/minor_introns.tmp ${DIR}/output/${species}/${t}/amb_minor_introns.tmp > ${DIR}/output/${species}/${t}/U12_introns.tmp

# Bin major introns & non-canonical

awk -F"\t" -v OFS="\t" '{U2=$3; U12=$5; PPT=$7; BPS=$11; if (PPT>=80) print $0,"strong_PPT"; else if (PPT>=60) print $0,"weak_PPT"; else if (PPT<=50 && BPS>=70) print $0,"major_hybrid"; else print $0,"degenerate_major"}' ${DIR}/output/${species}/${t}/U2_introns.tmp > ${DIR}/output/${species}/${t}/U2_introns_binned.tmp
awk -F"\t" -v OFS="\t" '{U2=$3; U12=$5; PPT=$7; BPS=$11; if (PPT>=80) print $0,"strong_PPT"; else if (PPT>=60) print $0,"weak_PPT"; else if (PPT<=50 && BPS>=70) print $0,"major_hybrid"; else print $0,"degenerate_major"}' ${DIR}/output/${species}/${t}/non-canonical_introns.tmp > ${DIR}/output/${species}/${t}/non-canonical_introns_binned_1.tmp

# Bin minor introns & non-canonical

awk -F"\t" -v OFS="\t" '{U2=$3; U12=$5; PPT=$7; BPS=$11; if (BPS>=80) print $0,"strong_BPS"; else if (BPS>=60) print $0,"weak_BPS"; else if (BPS<=50 && PPT>=70) print $0,"minor_hybrid"; else print $0,"degenerate_minor"}' ${DIR}/output/${species}/${t}/U12_introns.tmp > ${DIR}/output/${species}/${t}/U12_introns_binned.tmp
awk -F"\t" -v OFS="\t" '{U2=$3; U12=$5; PPT=$7; BPS=$11; if (BPS>=80) print $0":strong_BPS"; else if (BPS>=60) print $0":weak_BPS"; else if (BPS<=50 && PPT>=70) print $0":minor_hybrid"; else print $0":degenerate_minor"}' ${DIR}/output/${species}/${t}/non-canonical_introns_binned_1.tmp > ${DIR}/output/${species}/${t}/non-canonical_introns_binned_2.tmp

# Make final calls
# major: ($13=="strong_major" || $13=="major") && ($14=="strong_PPT" || $14=="weak_PPT")
# major hybrid: (($13=="strong_major" || $13=="major") && $14=="major_hybrid")
# major-like: U2, otherwise
# minor: ($13=="strong_minor" || $13=="minor") && ($14=="strong_BPS" || $14=="weak_BPS")
# minor_hybrid: ($13=="strong_minor" || $13=="minor") && $14=="minor_hybrid"
# minor-like: U12, otherwise
# non-canonical: non-canonical

# FIELDS: intron_id, 5'ss, U2 5'ss score, 5'ss, U12 5'ss score, 3'ss, U2 3'ss score, U2 BPS, U2 BPS score, U12 BPS, U12 BPS score, term dint, 5'ss class, 3'ss class + final call

awk -F"\t" -v OFS="\t" '{if (($13=="strong_major" || $13=="major") && ($14=="strong_PPT" || $14=="weak_PPT")) print $0,"major"; else if (($13=="strong_major" || $13=="major") && $14=="major_hybrid") print $0,"major_hybrid"; else print $0,"major-like"}' ${DIR}/output/${species}/${t}/U2_introns_binned.tmp > ${DIR}/output/${species}/${t}/introns_binned.tmp
awk -F"\t" -v OFS="\t" '{if (($13=="strong_minor" || $13=="minor") && ($14=="strong_BPS" || $14=="weak_BPS")) print $0,"minor"; else if (($13=="strong_minor" || $13=="minor") && $14=="minor_hybrid") print $0,"minor_hybrid"; else  print $0,"minor-like"}' ${DIR}/output/${species}/${t}/U12_introns_binned.tmp >> ${DIR}/output/${species}/${t}/introns_binned.tmp
awk -F"\t" -v OFS="\t" '{print $0,"non-canonical"}' ${DIR}/output/${species}/${t}/non-canonical_introns_binned_2.tmp >> ${DIR}/output/${species}/${t}/introns_binned.tmp

sort -k 1,1 ${DIR}/output/${species}/${t}/introns_binned.tmp > ${DIR}/output/${species}/${t}/introns_binned.tsv

rm -f ${DIR}/output/${species}/${t}/*.tmp