from Bio import SeqIO
import gzip

with open("tmpxcmM8P/splitted_temp.fastq", "rU") as splitted:
    original=SeqIO.parse(gzip.open("SRR027957_2.fastq.gz", "rb"), "fastq")
    for record in SeqIO.parse(splitted, "fastq"):
        orecord=original.next()

