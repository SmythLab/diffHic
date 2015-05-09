##################################################################################################
##################################################################################################
##################################################################################################
# iter_map.py: written by Aaron Lun
# written 5 June, 2014
# last modified 7 August, 2014

import argparse
parser = argparse.ArgumentParser(description="""Perform alignment of Hi-C read pairs with iterative mapping using Bowtie2, a la Imakaev et al.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-G", type=str, dest="genome", help="path to bowtie2 genome index, prefix only", required=True)
parser.add_argument("-1", type=str, dest="fq1", help="path to FASTQ file containing the first read of the pair", required=True)
parser.add_argument("-2", type=str, dest="fq2", help="path to FASTQ file containing the second read of the pair", required=True)
parser.add_argument("-P", type=int, dest="phred", help="phred offset for the quality scores in the FASTQ files", choices=[33, 64], default=33)
parser.add_argument("-o", type=str, dest="output", help="path to output BAM file", required=True)
parser.add_argument("--cmd", type=str, dest="cmd", help="path to/command for bowtie2 executable, with or without flags", default="bowtie2")
args = parser.parse_args()

import os
from subprocess import call
import pysam
import tempfile
tmpdir=tempfile.mkdtemp(dir=".")

outdir=args.output
bottom=25
step=5

# Processing the bowtie2 flag.
import re
bwt2cmd = re.split("\s+", args.cmd.strip()) + ["-x", args.genome, "--phred"+str(args.phred)]

####################################################################################################

allnames=[]
dumpf=open(os.devnull, "w")
for x, curf in enumerate([args.fq1, args.fq2]):
	# Running through and checking for read length uniformity (also get rid of the newline).
	counter=0
	toplen=0
	for line in open(curf, 'rU'):
		counter+=1
		if counter==2:
			toplen = len(line)
		elif counter%4==2 and toplen!=len(line):
			raise ValueError, "varying read lengths are not supported"
	toplen-=1 

	# Resetting flags for first and second.
	if x==0:
		add=0x41
	else:
		add=0x81

	# Defining the assorted files.
	fname, fext = os.path.splitext(os.path.basename(curf))
	currentfq=os.path.join(tmpdir, "temp.fq")
	os.link(curf, currentfq)
	newfq=os.path.join(tmpdir, "next.fq")
	newsam=os.path.join(tmpdir, "temp.sam")
	allbam=os.path.join(tmpdir, fname+".bam")

	starting=True
	for newlen in xrange(bottom, toplen+step, step):
		thrown_away=max(toplen-newlen, 0)
		if call(bwt2cmd+["--reorder", "--very-sensitive", "-U", currentfq, "-S", newsam, "--trim3", str(thrown_away)], stdout=dumpf, stderr=dumpf):
			raise SystemError, "bowtie2 failed for iterative alignment"

		fqout=open(newfq, mode="w")
		fqin=open(currentfq, mode="rU")
		samin=pysam.Samfile(newsam, 'r')
		if starting:
			samout=pysam.Samfile(allbam, "wb", template=samin)
			starting=False

		# Gradually extending the read from the 3' end. We compare the old
		# FASTQ with the old sam, see which ones aligned and which ones didn't.
		for curalign in samin:
			try:
				curname=fqin.readline()
				curseq=fqin.readline()	
				curother=fqin.readline()
				curqual=fqin.readline()
			except StopIteration:
				raise ValueError, "ran out of reads in the FASTQ file"
			if re.split("\s+", curname.strip())[0][1:]!=curalign.qname:
				raise ValueError, "mismatch between read names in FASTQ and BAM file"
			
			# Keeping only mapped alignments that have positive MAPQ and no XS tag.	Full length
			# reads also require uniqueness, otherwise they are added as empties.
			if curalign.flag&0x4==0 and 'XS' not in dict(curalign.tags) and curalign.mapq > 0:
				if thrown_away:
					if curalign.flag&0x10:
						curalign.cigar = [(5, thrown_away)] + curalign.cigar
					else:
						curalign.cigar += [(5, thrown_away)]
				curalign.flag += add
				samout.write(curalign)
			elif thrown_away:
				fqout.write(curname)
				fqout.write(curseq)
				fqout.write(curother)
				fqout.write(curqual)
			else:
				curalign.flag = 0x4 + add
				curalign.mapq = 0
				samout.write(curalign)
		
		# Shuffling some files around.	
		fqout.close()
		fqin.close()
		samin.close()
		os.remove(currentfq)
		os.rename(newfq, currentfq)

	# Cleaning up.
	samout.close()
	os.remove(currentfq)
	os.remove(newsam)
	bsorted=os.path.join(tmpdir, "sorted")
	pysam.sort("-n", allbam, bsorted)
	os.rename(bsorted+".bam", allbam)
	allnames.append(allbam)

####################################################################################################
# Running through and merging the same BAM files together.

cmd=["-n", "-f", args.output]+allnames
pysam.merge(*cmd)
for x in allnames:
	os.remove(x)

# Mopping up.
import shutil
shutil.rmtree(tmpdir)
dumpf.close()

####################################################################################################
