***Set inputs: working directory, source directory***

	workdir=$1
	fastadir=$2	
	species=$3
	assembly=$4
	version=$5

	BioMartFile=${fastadir}/Species/${species}/${species}.${assembly}.${version}.BioMart.txt
	TranscriptInfo=${fastadir}/Species/${species}/${species}.${assembly}.${version}.BioMart_transcriptInfo.txt
	GenomeFile=${fastadir}/Species/${species}/${species}.${assembly}.dna.*.fa

	BedFiles=${fastadir}/Species/${species}/${species}.${assembly}.${version}.BedFiles

	# Set BedFile paths

	BedFileTranscripts=${BedFiles}/${species}.${assembly}.${version}.transcripts.bed
	BedFileExons=${BedFiles}/${species}.${assembly}.${version}.exons.bed
	BedFileIntrons=${BedFiles}/${species}.${assembly}.${version}.introns.bed

***Extract putative BP sequence from intron FASTA file***
	
	# Generate BedFile to extract intron FASTA (need to sequence from -3 at 5'ss to +1 at 3'ss), only include unique introns

	awk -F"\t" -v OFS="\t" -v species=${species} '{if ($6=="+") print $1,($2-4),($3+1),species"_"$1"_"$2"_"$3"_F",($3-$2+5),$6; else print $1,($2-2),($3+3),species"_"$1"_"$2"_"$3"_R",($3-$2+5),$6}' ${BedFileIntrons} | sort | uniq > ${BedFiles}/${species}_introns_getfasta.bed.tmp

	# Use bedtools getfasta to extract the entire intronic sequence

	bedtools getfasta -fi ${GenomeFile} -bed ${BedFiles}/${species}_introns_getfasta.bed.tmp -s -name | fasta_formatter -t > ${BedFiles}/${species}_introns_raw.fa.tmp

	# Extract potential BPS (-44 to -1), remove sequence error introns and introns less than 40bp

	awk -F"\t" -v OFS="\t" '{print substr($1,1,(index($1,"(")-1)),substr($2,1,12),substr($2,(length($2)-13),14),substr($2,(length($2)-44),45)}' ${BedFiles}/${species}_introns_raw.fa.tmp > ${BedFiles}/${species}_introns_ss.fa.tmp
	awk -F"\t" -v OFS="\t" '{if (index($2,"N")==0 && index($3,"N")==0 && index($4,"N")==0 && length($2)==12 && length($3)==14 && length($4)>=41) print $1,$4}' ${BedFiles}/${species}_introns_ss.fa.tmp | sort -k1,1 | uniq > ${workdir}/output/${species}/${species}.${assembly}.${version}-Introns_3SS.fa

	rm -f ${BedFiles}/*.tmp
