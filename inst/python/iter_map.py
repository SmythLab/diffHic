##################################################################################################
##################################################################################################
##################################################################################################
# iter_map.py: written by Aaron Lun
# written 5 June, 2014

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

# Setting up imports
import os
from subprocess import Popen, PIPE

from Bio import SeqIO
import pysam
import tempfile
tmpdir=tempfile.mkdtemp(dir=".")

import gzip

# Processing the bowtie2 flag.
import re
bwt2cmd = re.split("\s+", args.cmd.strip()) + ["-x", args.genome, "--phred"+str(args.phred)]

####################################################################################################

is_gz=re.compile("\\.gz$")
def open_handle(fname, write=False):
    """A function that returns a Gzip file handle or a normal file handle, depending on the file name."""
    if (is_gz.search(fname)):
        if write:
            mode='wb'
        else:
            mode='rb'

        return gzip.open(fname, mode)
    else:
        if write:
            mode='w'
        else:
            mode='rU'

        return open(fname, mode)

####################################################################################################

allnames=[]
dumpf=open(os.devnull, "w")
bottom=25
step=5

for x, curf in enumerate([args.fq1, args.fq2]):
    # Running through and checking for read length uniformity (also get rid of the newline).
    counter=0
    toplen=0

    with open_handle(curf) as in_handle: 
        for record in SeqIO.parse(in_handle, "fastq"):
            if not counter:
                toplen=len(record)
            elif toplen!=len(record):
                raise ValueError("varying read lengths are not supported")

    # Resetting flags for first and second.
    if x==0:
        add=0x41
    else:
        add=0x81

    # Defining the assorted files.
    fname, fext = os.path.splitext(os.path.basename(curf))
    currentfq=os.path.join(tmpdir, "temp.fq")
    newfq=os.path.join(tmpdir, "next.fq")
    newsam=os.path.join(tmpdir, "temp.sam")
    allbam=os.path.join(tmpdir, fname+".bam")

    starting=True
    for newlen in xrange(bottom, toplen+step+1, step):
        thrown_away=max(toplen-newlen, 0)

        # Aligning reads.
        if starting:
            read_src=curf
        else:
            read_src=currentfq

        map_proc=Popen(bwt2cmd+["--reorder", "--very-sensitive", "-U", read_src, "-S", newsam, "--trim3", str(thrown_away)], stdout=dumpf, stderr=PIPE)
        map_out, map_err=map_proc.communicate()
        if map_proc.returncode:
            raise SystemError("bowtie2 failed for iterative alignment\n"+map_err)

        # Running through the reads and figuring out which ones aligned (or did not).
        with open_handle(read_src) as fqin, \
                open(newfq, mode="w") as fqout, \
                pysam.Samfile(newsam, 'r') as samin:

            inparse=SeqIO.parse(fqin, "fastq")
            if starting:
                samout=pysam.Samfile(allbam, "wb", template=samin)
                starting=False

            # Gradually extending the read from the 3' end. We compare the old
            # FASTQ with the old sam, see which ones aligned and which ones didn't.
            for curalign in samin:
                record=next(inparse)
                if record.id!=curalign.qname:
                    raise ValueError("mismatch between read names in FASTQ and BAM file")
                
                # Keeping only mapped alignments that have positive MAPQ and no XS tag.    Full length
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
                    SeqIO.write(record, fqout, "fastq")
                else:
                    curalign.flag = 0x4 + add
                    curalign.mapq = 0
                    samout.write(curalign)
            
        # Shuffling some files around.    
        os.rename(newfq, currentfq)

    # Cleaning up.
    samout.close()
    os.remove(currentfq)
    os.remove(newsam)
    bsorted=os.path.join(tmpdir, "sorted.bam")
    pysam.sort("-n", allbam, "-o", bsorted)
    os.rename(bsorted, allbam)
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
