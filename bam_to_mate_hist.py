#!/usr/bin/python

# takes a bam file and makes a histogram of distances between mate alignments to the 
# reference assembly
# takes the first 1M read pairs by default
### USAGE:
# python bam_to_mate_hist.py -b <BAM_FILE> -n <NUM_READS_TO_USE>
# creates a file Read_mate_dist.pdf in the working directory with relevant plots.

import sys
import pysam
import numpy as np
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args(desc):
	'''parse command-line args
	
	Args: 
		desc(str): program description, e.g. __file__
	
	Returns:
		args (dict): dict of the form {arg_name: arg_value}
	'''
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument("-n", "--num_reads", default=1000000,
                        help="Number of reads from bam file to use. default: 1000000")
	parser.add_argument("-b", "--bam_file", required=True, 
	help = "BAM file to evaluate for QC")
	args = parser.parse_args()
	return vars(args)

def parse_bam_file(bamfile_handle, num_reads):
	'''Parse a bam file, collect distances between mates in each read pair
	Args:
		bamfile_handle (str): path to bamfile to read
		num_reads (int): maximum number of reads to use from file (if there are so many 
		reads in the file). 
	Returns:
		dists (numpy array of ints): distances between mates.
	'''
	num = 0
	dists = np.empty([num_reads,1], dtype=int)
	
	last_read = ""
	bamfile = pysam.AlignmentFile(bamfile_handle, 'rb')
	for read in bamfile:
		if num >= num_reads:
			break
		if read.qname == last_read:
			continue
		last_read = read.qname
		ref1 = read.reference_id
		ref2 = read.next_reference_id
		if ref1 != ref2:
			continue
		read1_pos = read.reference_start
		read2_pos = read.next_reference_start
		dist = str(abs(read1_pos - read2_pos))
		dists[num] = dist
		num += 1

	dists = dists[0:num+1]
	return dists

def make_histograms(dists, bamfile_handle):
	'''make the read distance histograms
	Args:
		dists (numpy array of ints): Distances to plot in histogram.
		bamfile_handle (str): path to bamfile of dists
		'''
	num_reads = len(dists)
	with PdfPages("Read_mate_dist.pdf") as pdf:
		fig1 = plt.figure()
		plt.hist(dists, bins=40)
		plt.yscale("log", nonposy="clip")
		plt.title("Mate distance distribution for sample\n" + bamfile_handle+"\nfor first " + str(num_reads)+ " reads")
		plt.xlabel("Distance between read pair mates in Hi-C mapping")
		plt.ylabel("Number of reads")
		pdf.savefig(fig1)
		plt.close()

		fig2 = plt.figure()
		plt.hist(dists, bins=xrange(0,20000, 500))
		ax = fig2.add_subplot(111)
		ax.set_xlim(0, 20000)
		plt.yscale("log", nonposy="clip")
		plt.title("Mate distance distribution for sample\n" + bamfile_handle+"\nfor first " + str(num_reads)+ " reads")
		plt.xlabel("Distance between read pair mates in Hi-C mapping")
		plt.ylabel("Number of reads")
		pdf.savefig(fig2)
		plt.close()

if __name__ == "__main__":
	c_args = parse_args(__file__)
	num_reads = int(c_args["num_reads"])
	bamfile_handle = c_args["bam_file"]
	print "parsing the first {0} reads in bam file {1} to QC Hi-C library quality, plots"\
	" are written to ./Read_mate_dist.pdf.".format(
		num_reads, bamfile_handle
		)
	dists = parse_bam_file(num_reads=num_reads, bamfile_handle=bamfile_handle)
	print "Counts of zero distances (many is a sign of bad prep)"
	unique, counts = np.unique(dists, return_counts=True)
	print dict(zip(unique, counts))[0], "of total", len(dists)
	make_histograms(dists=dists, bamfile_handle=bamfile_handle)
