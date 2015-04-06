##################################################################################################
##################################################################################################
##################################################################################################
# presplit_map.py: written by Aaron Lun
# written 7 February, 2014
# last modified 7 August, 2014

import argparse
class myformat(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
	pass

parser = argparse.ArgumentParser(description="""Align chimeric reads from a Hi-C experiment using a presplitting approach: 
 a) Reads are split into segments at a known 'ligation signature' using cutadapt
 b) Segments are aligned with Bowtie2 in end-to-end mode with very high sensitivity 
 c) All results for each pair of libraries are stitched into a single BAM file
""", 
epilog="""Notes:
 - No processing is done to remove non-uniquely mapped reads. Filtering on MAPQ 
   or other associated tags (e.g., XS) is recommended before downstream analyses.
 - The segment corresponding to the 5' end of each read will always be considered the 
   primary alignment (i.e., all other segments are marked as secondary alignments).
""", formatter_class=myformat)
parser.add_argument("-G", type=str, dest="genome", help="path to bowtie2 genome index, prefix only", required=True)
parser.add_argument("-1", type=str, dest="fq1", help="path to FASTQ file containing the first read of the pair", required=True)
parser.add_argument("-2", type=str, dest="fq2", help="path to FASTQ file containing the second read of the pair", required=True)
parser.add_argument("-P", type=int, dest="phred", help="phred offset for the quality scores in the FASTQ files", choices=[33, 64], default=33)
parser.add_argument("--sig", type=str, dest="ligation", help="sequence of ligation signature (i.e. after blunt-end ligation of filled-in overhangs)", required=True)
parser.add_argument("-o", type=str, dest="output", help="path to the output BAM file", required=True)
parser.add_argument("--cmd", type=str, dest="cmd", help="path to/command for bowtie2 executable, with or without flags", default="bowtie2")
parser.add_argument("--cut", type=str, dest="cut", help="path to/command for cutadapt executable, with or without flags", default="cutadapt")
args = parser.parse_args()

## Setting up odds and ends which might be necessary.
import pysam
from subprocess import call
import os
import tempfile
tmpdir=tempfile.mkdtemp(dir=".")

## Additional processing of arguments before proceeding.
ligseq=args.ligation
liglen=len(ligseq)
if liglen%2!=0:
	raise ValueError, "ligation signature must be even in length"

## Assembling the bowtie2 command. 
import re
bwtcmd=re.split("\s+", args.cmd)+["--phred"+str(args.phred), "-x", args.genome]
cutcmd=re.split("\s+", args.cut)
	
####################################################################################################
## Writing a function which decides where to add the hard clipping (passes by reference).

def add_hard_clip(xread, is5, altlen):
	if not xread.is_unmapped:
		xcigar=list(xread.cigar)
		if xread.is_reverse:
			addfront=is5
		else:
			addfront=not is5
		if addfront:
			xcigar.insert(0, (5, altlen))
		else:
			xcigar.append((5, altlen))	
		xread.cigar=tuple(xcigar)
	return

####################################################################################################

allnames=[]
dumpf=open(os.devnull, "w")
for x, curf in enumerate([args.fq1, args.fq2]):

	# Setting up custom names for the output file.
	fname, fext = os.path.splitext(os.path.basename(curf))
	is_first=(x==0)
	outbam=os.path.join(tmpdir, fname+".bam")
	allnames.append(outbam)
			
	# Trimming the FASTQ file prior to input.
	init_split=os.path.join(tmpdir, "splitted_temp.fq")
	if call(cutcmd+["-o", init_split, "-a", ligseq, curf], stdout=dumpf, stderr=dumpf):
		raise SystemError, "cutadapt failed to trim reads"

	# Runs through the temporary and the original, and adds back half the ligation signature length.
	# Also records the 3' sequence for realignment. If it wasn't split, then we record that 
	# in a separate file as well. We also keep a log telling us what order everything was in before
	# we started cutting things up.
	original=open(curf, "rU")
	tempUs=os.path.join(tmpdir, "unsplit.fq")
	temp5=os.path.join(tmpdir, "5.fq")
	temp3=os.path.join(tmpdir, "3.fq")
	outputUs=open(tempUs, "w")
	output5=open(temp5, "w")
	output3=open(temp3, "w")
	counter=1

	totalsplit=0
	for line in open(init_split, "rU"):
		line=line.strip()
		oline=original.next().strip()
		if counter%2==1:
			extra=line
		else:
			junction=len(line)+liglen/2
			if junction < len(oline):
				output5.write(extra+"\n"+oline[:junction]+"\n")
   			 	output3.write(extra+"\n"+oline[junction:]+"\n")
				totalsplit+=1
			else:
			 	outputUs.write(extra+"\n"+oline+"\n")
			if counter==4:
				counter=0
		counter+=1
	totalsplit/=2

	try:
		original.next()
		raise IOError, "mismatching number of lines in cutadapt output and original file"	   			   
	except StopIteration:
		pass
	 	
	original.close()
	output5.close()
	output3.close()
	outputUs.close()
	os.remove(init_split)

	# Mapping both ends, separately. Need reordering to stitch split reads back together.
	mapped5=os.path.join(tmpdir, "5.sam")
	mapped3=os.path.join(tmpdir, "3.sam")
	for dex in xrange(2):
		if not dex: 
			curfile=temp5
			curout=mapped5
		else:
			curfile=temp3
			curout=mapped3
		if call(bwtcmd+["--reorder", "--very-sensitive", "-U", curfile, "-S", curout], stdout=dumpf, stderr=dumpf):
			raise SystemError, "bowtie2 failed for presplit alignment"
		os.remove(curfile)

	# Mapping any unsplit reads with local alignment.		 	
	mappedUs=os.path.join(tmpdir, "us.sam")
	if call(bwtcmd+["--local", "--very-sensitive-local", "-U", tempUs, "-S", mappedUs], stdout=dumpf, stderr=dumpf):
		raise SystemError, "bowtie2 failed for unsplit alignment"

	# Running through and re-integerating everything into a single file.
	# We need to set appropriate flags for first/second read of a pair, 
	# and for primary/secondary alignments. 5' reads are marked as primaries.
	# We also have to re-integrate it with respect to the unsplit reads, 
	# so it's effectively a merger of three files into one.
	sin3=pysam.Samfile(mapped3, "r")
	sin5=pysam.Samfile(mapped5, "r")
	sout=pysam.Samfile(outbam, "wb", template=sin3)

	for idx in xrange(totalsplit):
		try:
			read3=sin3.next()
			read5=sin5.next()	
		except StopIteration:
			raise IOError, "ran out of split reads to add"
		if read5.qname!=read3.qname:
			raise ValueError, "names of split segments do not match"
		len5=len(read5.seq)
		len3=len(read3.seq)

		read3.is_paired=read5.is_paired=True
		read3.is_read1=read5.is_read1=is_first
		read5.is_read2=read5.is_read2=not is_first

		add_hard_clip(read5, True, len3)
		sout.write(read5)
		read3.is_secondary=True
		add_hard_clip(read3, False, len5)
		sout.write(read3)
	sin5.close()
	sin3.close()

	for okread in pysam.Samfile(mappedUs, "r"):
		okread.is_paired=True	
		okread.is_read1=is_first
		okread.is_read2=not is_first	
		sout.write(okread)
	sout.close()

	os.remove(mapped5)
	os.remove(mapped3)
	os.remove(mappedUs)

	# Sorting the resulting file. We do it this way because the original
	# file wasn't guaranteed to be sorted. To try to merge the three files
	# (split, 5' and 3'), we'd need to assume some sorting order and I'm
	# not willing to do that. Sorting afterward enforces the 'samtools'
	# name ordering over anything that might have been there originally.
	bsorted=os.path.join(tmpdir, "sorted")
	pysam.sort("-n", outbam, bsorted)
	os.rename(bsorted+".bam", outbam)
	
####################################################################################################
# Stitching two files together to reform a single BAM file.

cmd=["-n", "-f", args.output]+allnames
pysam.merge(*cmd)
for x in allnames:
	os.remove(x)

# Mopping up.
import shutil
shutil.rmtree(tmpdir)
dumpf.close()

####################################################################################################


