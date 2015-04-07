###############################################################
# This script is design to obtain the BAM files from SRA files.
###############################################################

set -e
set -u

# <<ASSUMPTION>>: the relevant files have already been obtained from GEO:
#	Dekker et al. = GSE18199
#	Sofueva et al. = GSE49017
#	Rickman et al. = GSE37752

dekker=(SRR027957.sra SRR027958.sra)
sofueva=(SRR941267.sra) # SRR941268.sra SRR941269.sra SRR941270.sra SRR941271.sra SRR941272.sra SRR941273.sra SRR941274.sra SRR941275.sra SRR941276.sra SRR941277.sra SRR941278.sra SRR941279.sra SRR941280.sra SRR941281.sra SRR941282.sra) 
rickman=(SRR493818.sra) # SRR493819.sra SRR493820.sra SRR493821.sra SRR493822.sra SRR493823.sra SRR493824.sra  SRR493825.sra SRR493826.sra SRR493827.sra)

# <<ASSUMPTION>>: Bowtie2 and cutadapt are installed.

bowcmd=bowtie2
cutcmd=cutadapt

# <<ASSUMPTION>>: hg19 and mm10 indices have been built.
#	This can be done by making a FASTA file from a BSGenome object:
#
#   > bs <- BSGenome.Mmusculus.UCSC.mm10
#	> outfile <- "mm10.fa"
#	> for (chr in seqnames(bs)) {
#	+     y <- getSeq(bs, names = chr, start = 1, end = length(bs[[chr]]))
#	+     y <- DNAStringSet(y)
#	+     names(y) <- chr
#	+     writeXStringSet(filepath = outfile, y, append = TRUE)
#	+ }
#
#   ... and then running "bowtie2-build mm10.fa mm10_index/mm10" in the shell.
#	Alternatively, you can use your own FASTA file.

hg19=hg19_index/hg19
mm10=mm10_index/mm10

# <<ASSUMPTION>>: FixMateInformation and MarkDuplicates are available from the Picard suite.

fixcmd=FixMateInformation
markcmd=MarkDuplicates

# <<ASSUMPTION>>: samtools has been installed.

samcmd=samtools

# <<ASSUMPTION>>: diffHic has been installed on the relevant R installation.

rcmd=R

# <<ASSUMPTION>>: fastq-dump (from NCBI's SRA toolkit) has been installed.

fqdcmd=fastq-dump

###############################################################
# We pull out the Hi-C mapping script from diffHic.

mapper=hicmap.py
echo "require(diffHic); file.copy(system.file('python', 'presplit_map.py', package='diffHic'), '$mapper');" | ${rcmd} --no-save
if [ ! -e $mapper ]; then
	echo "Extraction of Hi-C mapping script failed."
	exit 1
fi

vtmp=`mktemp -d --tmpdir=.`
tmpfile=blah.txt

###############################################################
# Running through each set of files.

ligation=AAGCTAGCTT # all of them use the same ligation signature.

for i in {1..3}
do 
	if [[ $i -eq 1 ]]; then
		allfiles=(${dekker[@]})
		genome=$hg19
		phred=33
	elif [[ $i -eq 2 ]]; then
		allfiles=(${sofueva[@]})
		genome=$mm10
		phred=64
	else 
		allfiles=(${rickman[@]})
		genome=$hg19
		phred=33
	fi

	for sra in "${allfiles[@]}"
	do 
		${fqdcmd} --split-files $sra
		prefix=`echo $sra | sed "s/\.sra$//g"`
		read1=${prefix}_1.fastq
		read2=${prefix}_2.fastq

		rawbam=temp.bam
		python $mapper -o $rawbam -G $genome -1 $read1 -2 $read2 --sig $ligation -P $phred --cmd "${bowcmd} -p 8" --cut $cutcmd

		fixbam=fixed.bam
		${fixcmd} I=$rawbam O=$fixbam TMP_DIR=$vtmp VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate

		markbam=marked.bam
		${markcmd} I=$fixbam O=$markbam M=$tmpfile TMP_DIR=$vtmp AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
		${samcmd} sort -n $markbam $prefix

		rm $rawbam $fixbam $markbam
		rm $read1 $read2
	done
done

###############################################################
# Mopping up.

rm $tmpfile
rm -r $vtmp

###############################################################
