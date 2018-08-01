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
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
import pdfkit
import markdown as md


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
                        help="BAM file to evaluate for QC")
    parser.add_argument("--count_diff_refname_stub", default=False, action="store_true",
                        help="Can be used to QC differently given reference name formatting. For internal use.")
    parser.add_argument("--outfile_name", "-o", default="Read_mate_dist", type=str,
                        help="Path to which to write plots to (PNG suffix will be attached).")
    parser.add_argument("--make_report", "-r", default=False, action="store_true",
                        help="Whether to export results in a PDF report. Requires that the QC script be" \
                             "in the same directory as the QC repo's collateral directory. Default: False.")

    args = parser.parse_args()

    return vars(args)


def parse_bam_file(bamfile_handle, num_reads, count_diff_refname_stub=False):
    '''Parse a bam file, collect distances between mates in each read pair
    Args:
        bamfile_handle (str): path to bamfile to read
        num_reads (int): maximum number of reads to use from file (if there are so many
        reads in the file).
        count_from_stub (bool): whether to
    Returns:
        diff_chr (int): number of reads mapping between contigs/chromosomes
        dists (numpy array of ints): distances between mates.
    '''
    num = 0
    dists = np.empty([num_reads, 1], dtype=int)
    last_read = ""
    bamfile = pysam.AlignmentFile(bamfile_handle, 'rb')
    refs = bamfile.references
    #n50, total_len = calc_n50_from_header(bamfile.header)
    #print "assembly N50: {0} bp, assembly length: {1} bp".format(n50, total_len)
    diff_chr = 0
    diff_stub = 0  # if reference name is trimmed back to "." delim, how many among such?
    split_reads = 0
    dupe_reads = 0
    for read in bamfile:
        if num >= num_reads:
            break
        # if (num % 10000) == 0:
        #	print num
        # read 2 gives info only for split reads
        if read.qname == last_read:
            if "S" in read.cigarstring:
                split_reads += 1
            if read.flag == "1024":
                dupe_reads += 1
            continue
        last_read = read.qname
        ref1 = read.reference_id
        ref2 = read.next_reference_id
        # count split reads
        if "S" in read.cigarstring:
            split_reads += 1

        if ref1 != ref2:
            diff_chr += 1
            if count_diff_refname_stub:
                ref1_stub = refs[ref1].split(".")[0]
                ref2_stub = refs[ref2].split(".")[0]
                # print ref1_stub, ref2_stub
                if ref1_stub != ref2_stub:
                    diff_stub += 1

        if read.flag == "1024":
            dupe_reads += 1

        else:
            read1_pos = read.reference_start
            read2_pos = read.next_reference_start
            dist = str(abs(read1_pos - read2_pos))
            dists[num] = dist
        num += 1

    dists = dists[0:num + 1]
    return diff_chr, dists, diff_stub, split_reads, dupe_reads


def calc_n50_from_header(header, xx=50.0):
    '''calculate the N50 of the starting assembly from the information in a pysam header object.

        Args:
            header (pysam.AlignmentFile.header): the header from which to extract reference sequence lengths
            xx (float): the XX of NXX. we assume it is N50 so the default is 50.0

        Returns:
            contig_len (int): the NXX (probably N50) of the assembly based on the header
            total (int): the total length of the assembly

    '''
    frac = xx / 100.0
    lens = [contig["LN"] for contig in header["SQ"]]
    lens.sort()
    total = sum(lens)
    nxx_len = total * frac

    cumsum = 0
    contig_len = 0
    for length in reversed(lens):
        contig_len = length
        cumsum += length
        if cumsum >= nxx_len:
            break
    return contig_len, total


def make_histograms(dists, bamfile_handle, num_reads, outfile_name):
    '''make the read distance histograms using matplotlib and write them to disk.
	Args:
		dists (numpy array of ints): Distances to plot in histogram.
		bamfile_handle (str): path to bamfile of dists
	'''
    num_reads = len(dists)
    # with PdfPages(outfile_name) as pdf:
    fig1 = plt.figure()
    plt.hist(dists, bins=40)
    ax = fig1.add_subplot(111)
    ax.set_ylim(0.5, num_reads * 2)
    plt.yscale("log", nonposy="clip")
    plt.title("\nMate distance distribution for first " + str(num_reads) + " reads for sample\n" + bamfile_handle)
    plt.xlabel("Distance between read pair mates in Hi-C mapping")
    plt.ylabel("Number of reads")
    fig1.savefig(outfile_name + "_long.png")
    plt.close(fig1)

    fig2 = plt.figure()
    plt.hist(dists, bins=xrange(0, 20000, 500))
    ax = fig2.add_subplot(111)
    ax.set_xlim(0, 20000)
    ax.set_ylim(0.5, num_reads * 2)
    plt.yscale("log", nonposy="clip")
    plt.title("Mate distance distribution for first " + str(num_reads) + " reads for sample\n" + bamfile_handle)
    plt.xlabel("Distance between read pair mates in Hi-C mapping")
    plt.ylabel("Number of reads")
    fig2.savefig(outfile_name + "_short.png")
    plt.close(fig2)


def make_pdf_report(qc_repo_path, stat_dict, outfile_name):
    '''Make the pdf report using the template, the

    Args:
        qc_repo_path (str): path the qc tool repo, which holds markdown and style templates
        stat_dict ({str: str/float}): a dictionary of relevant data (stats, file paths) regarding the lib to be used in the report.
        outfile_name (str): path to which the report shall be written.

    '''
    options = {
        'page-size': 'A4',
        'dpi': 350,
        'margin-top': '0.75in',
        'margin-right': '0.75in',
        'margin-bottom': '0.75in',
        'margin-left': '0.75in',
        'encoding': "UTF-8",
        'custom-header': [
            ('Accept-Encoding', 'gzip')
        ],
        'no-outline': None
    }

    template_path = os.path.join(qc_repo_path, "collateral", "HiC_QC_report_template.md")
    style_path = os.path.join(qc_repo_path, "collateral", "style.css")
    if not os.path.exists(template_path):
        UserWarning("can't find markdown template at {0}! skipping making template.".format(
            qc_repo_path)
        )
        sys.exit()

    with open(template_path) as file:
        str = file.read()
        sub_str = str.format(**stat_dict)  # splat the statistics and path into the markdown, render as html
        # print sub_str
        html = md.markdown(sub_str, extensions=['tables', 'nl2br'])

        # write out just html
        html_out = open(outfile_name + "_qc_report.html", 'w')
        html_out.write(html)
        html_out.close()

        # print html
        pdfkit.from_string(html, outfile_name + "_qc_report.pdf", options=options, css=style_path)


def extract_stats(stat_list, bamfile_handle, outfile_name, count_diff_refname_stub=False):
    '''Make a dict of results and data suitable to be passed to the report template for splatting.

    Args:
        stat_list ([str]): the four statistics estimated from the bam file
                           (stat_list = [diff_chr, dists, diff_stub, split_reads])
        bamfile_handle (str): path to bam file
        outfile_handle (str): path to where output files will be written
        count_diff_refname_stub (bool): whether we are counting the contig name stub differences.

    Returns:
        stat_dict ({str:float/str}: mappings of the stats and data suitable to be consumed by report generator
    '''
    stat_dict = {}
    stat_dict["BAM_FILE_PATH"] = os.path.split(bamfile_handle)[-1]
    num_pairs = len(stat_list[1])
    stat_dict["NUM_PAIRS"] = num_pairs

    print "Histograms written to:", os.path.abspath(outfile_name + "_long.png"), os.path.abspath(outfile_name + "_short.png")
    stat_dict["PATH_TO_LONG_HIST"] = os.path.abspath(outfile_name + "_long.png")
    stat_dict["PATH_TO_SHORT_HIST"] = os.path.abspath(outfile_name + "_short.png")

    print "Counts of zero distances (many is a sign of bad prep):"
    unique, counts = np.unique(stat_list[1], return_counts=True)  # tabulates the distances, with indices as the dists
    # print unique[-100:-1]
    zero_dist = dict(zip(unique, counts))[0]  # first element is the zero-distances
    stat_dict["ZERO_DIST_PAIRS"] = "{0:.3f}".format(float(zero_dist) / num_pairs)
    print zero_dist, "of total", num_pairs, "fraction ", stat_dict["ZERO_DIST_PAIRS"]

    above_10k = len([dist for dist in stat_list[1] if dist > 10000])
    stat_dict["NUM_10KB_PAIRS"] = "{0:.3f}".format(float(above_10k) / num_pairs)
    print "Count of read pairs with distance > 10KB (many is a sign of good prep):"
    print above_10k, "of total", len(dists), ", fraction ", stat_dict["NUM_10KB_PAIRS"]

    stat_dict["NUM_DIFF_CONTIG_PAIRS"] = "{0:.3f}".format(float(stat_list[0]) / num_pairs)
    print "Count of read pairs with mates mapping to different chromosomes/contigs (sign of good prep IF same genome):"
    print stat_list[0], "of total", num_pairs, ", fraction ", stat_dict["NUM_DIFF_CONTIG_PAIRS"]

    stat_dict["NUM_SPLIT_READS"] = "{0:.3f}".format(stat_list[3] / float(num_pairs * 2))
    print "Count of split reads (more is usually good, as indicates presence of Hi-C junction in read):"
    print stat_list[3], "of total", num_pairs * 2, ", fraction ", stat_dict["NUM_SPLIT_READS"]

    stat_dict["NUM_DUPE_READS"] = "{0:.3f}".format(stat_list[4] / float(num_pairs * 2))
    print "Count of duplicate reads (bad; WILL ALWAYS BE ZERO UNLESS BAM FILE IS PREPROCESSED TO SET THE DUPLICATES FLAG):"
    print stat_list[4], "of total", num_pairs * 2, ", fraction ", stat_dict["NUM_DUPE_READS"]

    if count_diff_refname_stub:
        print "Count of read pairs with mates mapping to different reference groupings, e.g. genomes (sign of bad " \
              "prep potentially):"
        print stat_list[2], "of total", num_pairs, ", fraction", float(stat_list[2]) / num_pairs

    return stat_dict


if __name__ == "__main__":
    c_args = parse_args(__file__)
    num_reads = int(c_args["num_reads"])
    bamfile_handle = c_args["bam_file"]
    outfile_name = c_args["outfile_name"]
    make_report = c_args["make_report"]

    count_diff_refname_stub = c_args["count_diff_refname_stub"]
    print "parsing the first {0} reads in bam file {1} to QC Hi-C library quality, plots" \
          " are written to {2}*".format(
        num_reads, bamfile_handle, outfile_name
    )
    diff_chr, dists, diff_stub, split_reads, dupe_reads = parse_bam_file(num_reads=num_reads, bamfile_handle=bamfile_handle,
                                                             count_diff_refname_stub=count_diff_refname_stub)

    stat_list = [diff_chr, dists, diff_stub, split_reads, dupe_reads]
    script_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

    stat_dict = extract_stats(stat_list=stat_list, bamfile_handle=bamfile_handle, outfile_name=outfile_name,
                              count_diff_refname_stub=count_diff_refname_stub)

    make_histograms(dists=dists, bamfile_handle=bamfile_handle, num_reads=num_reads, outfile_name=outfile_name)

    if make_report:
        make_pdf_report(qc_repo_path=script_path, stat_dict=stat_dict, outfile_name=outfile_name)
