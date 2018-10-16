#!/usr/bin/python

# takes a bam file and makes a histogram of distances between mate alignments to the
# reference assembly
# takes the first 1M read pairs by default

### USAGE:
# python bam_to_mate_hist.py -b <BAM_FILE> -n <NUM_READS_TO_USE> -o <outfile_stub>
# creates files in the working directory with relevant plots, also text files of statistics.
# flip -r flag  (assuming you have dependencies) to make a PDF report with everything together.

from __future__ import print_function
from __future__ import division

import sys
import pysam
import numpy as np
import argparse
import os
import matplotlib
from collections import Counter

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pdfkit
import markdown as md
from scipy import optimize

def saturation(x, V, K):
    return V * x / (x + K)

def calc_nxx(header, xx=50):
    '''Calculate the NXX (typically N50) of an assembly given a pysam.AlignmentHeader object.

        Args:
            header (pysam.AlignmentHeader): the header from which to extract reference sequence lengths
            xx (float): the XX of NXX. we assume it is N50 so the default is 50.

        Returns:
            contig_len (int): the NXX (probably N50) of the assembly
            total (int): The total length of contigs in the assembly
    '''

    frac = xx / 100.0
    lens = [contig['LN'] for contig in header['SQ']]
    lens.sort()
    total = sum(lens)
    nxx_len = total * frac

    cumsum = 0
    contig_len = 0
    for length in reversed(lens):
        contig_len = int(length)
        cumsum += length
        if cumsum >= nxx_len:
            break

    return contig_len, total

class HiCQC(object):
    def __init__(self):
        self.per_read_metrics = set(['total_reads', 'unmapped_reads', 'split_reads', 'duplicate_reads', 'mapq0_reads'])
        self.per_pair_metrics = set(['intercontig_pairs',
                                     'different_ref_stub_pairs',
                                     'pairs_greater_10k',
                                     'zero_dist_pairs',
                                     'total_pairs_on_same_contig',
                                     'pairs_on_contigs_greater_10k',
                                     'pairs_greater_10k_on_contigs_greater_10k',
                                     'pairs_on_same_strand',
                                     'total_read_pairs'
                                     ])
        # Dictionary of key --> numerator, denominator pairs for stringify_stats
        self.to_percents =    {
                               'perc_zero_dist_pairs': ('zero_dist_pairs', 'total_read_pairs'),
                               'perc_pairs_greater_10k': ('pairs_greater_10k', 'total_read_pairs'),
                               'perc_pairs_greater_10k_on_contigs_greater_10k': ('pairs_greater_10k_on_contigs_greater_10k', 'pairs_on_contigs_greater_10k'),
                               'perc_intercontig_pairs': ('intercontig_pairs', 'total_read_pairs'),
                               'perc_split_reads': ('split_reads', 'total_reads'),
                               'perc_duplicate_reads': ('duplicate_reads', 'total_reads'),
                               'perc_mapq0_reads': ('mapq0_reads', 'total_reads'),
                               'perc_unmapped_reads': ('unmapped_reads', 'total_reads'),
                               'perc_different_ref_stub_pairs': ('different_ref_stub_pairs', 'total_read_pairs'),
                               'perc_pairs_on_same_strand': ('pairs_on_same_strand', 'total_pairs_on_same_contig')
                               }

        self.convert_to_pairs = set(['unmapped_reads', 'split_reads', 'duplicate_reads', 'mapq0_reads'])

        self.paths = {'script_dir': os.path.dirname(__file__)}

        self.N50 = None
        self.stats = Counter()
        self.dists = Counter()
        self.total_array = []
        self.non_dup_array = []

    def parse_bam(self, bamfile, max_read_pairs=-1):
        self.paths['bamfile'] = bamfile
        self.paths['bamname'] = os.path.basename(bamfile)

        a = None
        b = None
        i = 0

        with pysam.AlignmentFile(self.paths['bamfile']) as bam_fh:
            self.extract_header_info(bam_fh.header)

            for read in bam_fh:
                if read.is_secondary or read.is_supplementary:
                    continue
                if a is None:
                    a = read
                    continue
                if max_read_pairs != -1 and i / 2 > max_read_pairs:
                    break

                if read.query_name == a.query_name:
                    b = read
                    self.process_pair(a, b)

                    if i % 1000 == 0:
                        self.update_dup_stats()
                else:
                    a = read
                i += 1

        self.finalize_stats()

    def extract_header_info(self, header, xx=50):
        '''Extract reference names, calculate N50, get total assembly length, and get set of contigs > 10kbp from a pysam header.
        Also checks if input bamfile is labeled as coordinate sorted and throws a ValueError if True.

            Args:
                header (pysam.AlignmentFile.header): the header from which to extract reference sequence lengths
                xx (float): the XX of NXX. we assume it is N50 so the default is 50

            Sets:
                self.refs dict(int-->str): ref_id --> ref_name mappings for the assembly
                self.contig_len (int): the NXX (probably N50) of the assembly based on the header
                self.total (int): the total length of the assembly
                self.contigs_greater_10k (set(str)): The set of names of contigs with length > 10Kbp
            Raises:
                ValueError if header labels bamfile as coordinate sorted
        '''

        if 'HD' in header and 'SO' in header['HD'] and header['HD']['SO'].lower().strip() == 'coordinate':
            raise ValueError('Error: bam_to_mate_hist.py requires read name sorted input, but bamfile {} is coordinate sorted'.format(self.paths['bamfile']))

        self.refs = header.references
        self.stats['num_refs'] = len(self.refs)

        self.N50, self.total_length = calc_nxx(header)
        self.contigs_greater_10k = set([contig['SN'] for contig in header['SQ'] if contig['LN'] > 10000])

    def process_pair(self, a, b):
        self.stats['total_read_pairs'] += 1

        self.update_read_stats(a)
        self.update_read_stats(b)

        if not a.is_unmapped and not b.is_unmapped:
            self.update_mapped_pair_stats(a, b)

    def update_read_stats(self, read):
        self.stats['total_reads'] += 1
        if read.is_unmapped:
            self.stats['unmapped_reads'] += 1
        elif read.mapping_quality == 0:
            self.stats['mapq0_reads'] += 1
        if read.has_tag('SA'):
            self.stats['split_reads'] += 1
        if read.is_duplicate:
            self.stats['duplicate_reads'] += 1

    def update_mapped_pair_stats(self, a, b):
        if a.reference_name != b.reference_name:
            self.stats['intercontig_pairs'] += 1

            refa_stub = a.reference_name.split('.')[0]
            refb_stub = b.reference_name.split('.')[0]

            if refa_stub != refb_stub:
                self.stats['different_ref_stub_pairs'] += 1

        else:
            self.stats['total_pairs_on_same_contig'] += 1
            if (a.is_reverse and b.is_reverse) or (not a.is_reverse and not b.is_reverse):
                self.stats['pairs_on_same_strand'] += 1
            dist = abs(a.reference_start - b.reference_start)
            self.dists[dist] += 1

            if dist > 10000:
                self.stats['pairs_greater_10k'] += 1
            if dist == 0:
                self.stats['zero_dist_pairs'] += 1
            if a.reference_name in self.contigs_greater_10k:
                self.stats['pairs_on_contigs_greater_10k'] += 1
                if dist > 10000:
                    self.stats['pairs_greater_10k_on_contigs_greater_10k'] += 1

    def update_dup_stats(self):
        self.total_array.append(self.stats['total_reads'])
        self.non_dup_array.append(self.stats['total_reads'] - self.stats['duplicate_reads'])

    def finalize_stats(self):
        self.total_array = np.array(self.total_array)
        self.non_dup_array = np.array(self.non_dup_array)
        self.stats['proportion_pairs_greater_10k_on_contigs_greater_10k'] = self.stats['pairs_greater_10k_on_contigs_greater_10k'] / \
                                                                            self.stats['pairs_on_contigs_greater_10k']

    def plot_dup_saturation(self, outfile_prefix, target_x=100000000, min_sample=100000, target_y=None):
        '''Fit and plot a saturation curve from cumulative total and non-dup read counts.

            Args:
                outfile_prefix (str): prefix of output plot file.
                self.total_array (np.array): Numpy array of total read counts.
                self.non_dup_array (np.array): Numpy array of non-duplicate read counts (must be same shape as self.total_array)
                target_x (int): Total reads to extrapolate to.
                target_y (int): Non-dup reads to extrapolate to. If specified with target_x, plots a point.

            Returns:
                Observed duplication rate, extrapolated duplication rate, and target total reads
        '''

        self.stats['observed_dup_rate'] = -1
        self.stats['extrapolated_dup_rate'] = -1
        self.stats['target_read_total'] = target_x
        self.stats['dup_sat_V'] = -1
        self.stats['dup_sat_K'] = -1

        outfile = outfile_prefix + '.dup_saturation.png'
        self.paths['dup_sat_curve'] = outfile

        if self.total_array[-1] < min_sample:
            UserWarning('too few reads to estimate duplication rate (<{0})!!'.format(min_sample))
            fig, ax = plt.subplots(1)
            plt.title('Insufficient reads to estimate duplication rate!!!')
            plt.savefig(outfile)
            plt.close()
            return 0

        try:
            params, params_cov = optimize.curve_fit(saturation,
                                                    self.total_array,
                                                    self.non_dup_array,
                                                    p0=[self.total_array[-1],
                                                        self.non_dup_array[-1]/2],
                                                    maxfev=6000
                                                    )
        except RuntimeError as e:
            UserWarning('Convergence failed for duplicate curve fitting')
            fig, ax = plt.subplots(1)
            plt.title('Convergence failed for duplicate curve fitting!!!')
            plt.savefig(outfile)
            plt.close()
            return 0

        self.stats['dup_sat_V'] = params[0]
        self.stats['dup_sat_K'] = params[1]

        if target_x is not None:
            # print(args.target_x.dtype)
            coord_max = target_x
            t = np.linspace(0, target_x, 1000)
        elif target_y is not None:
            coord_max = target_y
            x_temp = self.total_array[-1]
            y_temp = saturation(x_temp, *params)
            while y_temp < args.target_y:
                x_temp += 1000000
                y_temp = saturation(x_temp, *params)
            t = np.linspace(0, x_temp, 1000)
        else:
            coord_max = self.total_array[-1]
            t = np.linspace(0, self.total_array[-1], 1000)

        fig, ax = plt.subplots(1)
        plt.plot(self.total_array, self.non_dup_array, 'k')
        plt.plot(t, saturation(t, *params), 'r-')

        self.stats['observed_dup_rate'] = 1-float(self.non_dup_array[-1])/self.total_array[-1]
        non_dup_rate = None

        if target_x is not None:
            print('At {} reads, estimated {:.0f} non-dup reads'.format(target_x, saturation(target_x, *params)))
            print('True non-dup reads: {}'.format(target_y))
            non_dup_rate = saturation(target_x, *params) / float(target_x)
            if target_y is not None:
                plt.plot(target_x, target_y, 'bo')
        self.stats['extrapolated_dup_rate'] = 1-non_dup_rate if non_dup_rate is not None else None

        patch = matplotlib.patches.Rectangle((0, 0), self.total_array[-1], self.non_dup_array[-1], fill=False, color='k')
        ax.add_patch(patch)
        plt.ylim(0, coord_max)
        plt.xlim(0, coord_max)
        if non_dup_rate is not None:
            plt.title('{}\nproportion duplicated (sampled): {:.2f}\nproportion duplicated (extrapolated): {:.2f}'.format(
                outfile, self.stats['observed_dup_rate'], 1-non_dup_rate))
        else:
            plt.title('{}\nproportion duplicated (sampled): {:.2f}'.format(
                outfile, self.stats['observed_dup_rate']))
        plt.xlabel('Total reads')
        plt.ylabel('Non-duplicate reads')
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

        print('Best V = {}, best K = {}'.format(self.stats['dup_sat_V'], self.stats['dup_sat_K']))

        return 0

    def plot_histograms(self, outfile_prefix):
        '''make the read distance histograms using matplotlib and write them to disk.
        Args:
            dists (dictionary of mate distances and counts): Distances to plot in histogram.
            num_pairs (int): number of read pairs analyzed
            bamfile (str): path to bamfile of dists
        '''

        num_dists = sum(self.dists.values())
        num_pairs = self.stats['total_read_pairs']
        title_string = '\nMate distance distribution for first {} read pairs for sample\n{}'.format(num_pairs,
                                                                                                    self.paths['bamname'])
        fig1 = plt.figure()
        plt.hist(list(self.dists.keys()), weights=list(self.dists.values()), bins=50)
        ax = fig1.add_subplot(111)
        ax.set_ylim(0.5, num_dists * 2)
        plt.yscale('log', nonposy='clip')
        plt.title(title_string)
        plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig)')
        plt.ylabel('Number of reads')
        long_hist_path = outfile_prefix + '_long.png'
        fig1.savefig(long_hist_path)
        plt.close(fig1)
        self.paths['long_hist'] = long_hist_path

        fig2 = plt.figure()
        plt.hist(list(self.dists.keys()), weights=list(self.dists.values()), bins=range(0, 20000, 500))
        ax = fig2.add_subplot(111)
        ax.set_xlim(0, 20000)
        ax.set_ylim(0.5, num_pairs * 2)
        plt.yscale('log', nonposy='clip')
        plt.title(title_string)
        plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig)')
        plt.ylabel('Number of reads')
        short_hist_path = outfile_prefix + '_short.png'
        fig2.savefig(short_hist_path)
        plt.close(fig2)
        self.paths['short_hist'] = short_hist_path

        fig3 = plt.figure()
        offset_dists = {}
        for key, value in self.dists.items():
            offset_dists[key+1] = value
        min_dist = min(offset_dists.keys())
        max_dist = max(offset_dists.keys())
        plt.hist(list(offset_dists.keys()),
                 weights=list(offset_dists.values()),
                 bins=np.logspace(np.log10(min_dist),
                                  np.log10(max_dist),
                                  50),
                 log=True)
        ax = fig3.add_subplot(111)
        ax.set_ylim(0.5, num_pairs * 2)
        plt.yscale('log', nonposy='clip')
        plt.xscale('log')
        plt.xlim(xmin=1)
        plt.title(title_string)
        plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig, log scale)')
        plt.ylabel('Number of reads (log scale)')
        plt.tight_layout()
        log_log_hist_path = outfile_prefix + '_log_log.png'
        fig3.savefig(log_log_hist_path)
        plt.close(fig3)
        self.paths['log_log_hist'] = log_log_hist_path

    def html_from_judgement(self):
        '''return a formatted HTML string based on two judgment bool values
            Args:
                good (bool): does the hi-c library show characteristics of 'goodness', e.g. many long-distance contacts etc.
                bad (bool): does the hi-c library show 'bad' characteristics, e.g. zero-distance reads or too many duplicates.
            Returns:
                str: an HTML string to be substituted into the report to subjectively grade the assembly. 4 possibilities.
            Raises:
                ValueError: if impossible logical situations occur given two bools.

        '''

        if self.judge_good and not self.judge_bad:
            self.judge_html = '<span class="pass">PASS</span>'
        elif not self.judge_good and self.judge_bad:
            self.judge_html = '<span class="fail">FAIL</span>'
        elif self.judge_good and self.judge_bad:
            self.judge_html = '<span class="mixed-results">MIXED RESULTS</span>'
        elif not self.judge_good and not self.judge_bad:
            self.judge_html = '<span class="low-signal">LOW SIGNAL</span>'
        else:
            raise ValueError('logical impossibility!')

    def pass_judgement(self):
        '''Pass judgement on the library according to certain mostly subjective ideas about what is good

        Args:
            stat_dict ({str: float/str}): mapping of lib characteristics to their values

        Returns:
            judgement (str): an HTML string to put into pass/fail box
        '''
        long_contacts = self.stats['pairs_greater_10k'] / self.stats['total_read_pairs'] > 0.05
        long_floor = self.stats['pairs_greater_10k'] / self.stats['total_read_pairs'] > 0.01
        useful_contacts = self.stats['intercontig_pairs'] / self.stats['total_read_pairs'] > 0.3
        low_contiguity = self.N50 < 100000
        many_zero_pairs = self.stats['zero_dist_pairs'] / self.stats['total_read_pairs'] > 0.1
        many_many_zero_pairs = self.stats['zero_dist_pairs'] / self.stats['total_read_pairs'] > 0.2
        high_dupe = (self.stats['duplicate_reads'] / self.stats['total_reads'] > 0.05 and self.stats['total_read_pairs'] <= 1e6) \
                    or (self.stats['duplicate_reads'] / self.stats['total_reads'] > 0.3 and self.stats['total_read_pairs'] <= 1e8)

        #print long_contacts, long_floor, useful_contacts, low_contiguity, many_zero_pairs, many_many_zero_pairs, high_dupe
        if (long_contacts or useful_contacts) and (low_contiguity or long_floor):
            good = True
        else:
            good = False

        bad = False
        if low_contiguity:
            if many_many_zero_pairs or high_dupe:
                bad = True
        else:
            if many_zero_pairs or high_dupe or not long_floor:
                bad = True
        self.judge_good = good
        self.judge_bad = bad
        self.html_from_judgement()

    def stringify_stats(self):
        '''Convert stats to output dictionary with pretty strings and percents.'''

        # Dict of key --> (value, fmt) pairs for items that aren't counts
        self.other_stats = {
                            'N50': (self.N50, '{:,}'),
                            'contigs': (len(self.refs), '{:,}'),
                            'contigs_greater_10k': (len(self.contigs_greater_10k), '{:,}'),
                            'total_length': (self.total_length, '{:,}'),
                            'total_reads': (self.stats['total_reads'], '{:,}'),
                            'target_read_total': (self.stats['target_read_total'], '{:,}'),
                            'extrapolated_dup_rate': (self.stats['extrapolated_dup_rate'] * 100 \
                                                      if self.stats['extrapolated_dup_rate'] > 0 \
                                                      else self.stats['extrapolated_dup_rate'], '{:.2f}%'),
                            'judgment': (self.judge_html, '{}')
                            }
        self.out_stats = {}
        for key, (num, denom) in self.to_percents.items():
            try:
                self.out_stats[key] = '{:.2f}%'.format((self.stats[num] / self.stats[denom]) * 100)
            except ZeroDivisionError as e:
                self.out_stats[key] = 'NaN'

        for item in self.convert_to_pairs:
            self.out_stats[item] = '{:,}'.format(self.stats[item] // 2)

        for item in self.per_pair_metrics:
            self.out_stats[item] = '{:,}'.format(self.stats[item])

        for item in self.paths:
            if item == "bamname":
                self.out_stats[item] = self.paths[item]
            else:
                self.out_stats[item] = os.path.abspath(self.paths[item])

        for key, (value, fmt) in self.other_stats.items():
            try:
                self.out_stats[key] = fmt.format(value)
            except:
                print(key, value, fmt)
                sys.exit()


    def print_stats(self, count_diff_refname_stub=False):
        '''Process a dict of results and data suitable to be passed to the report template for splatting.
        Now that dict processing is happening mostly outside of here,

        Args:
            stat_list ([str]): the four statistics estimated from the bam file
                               (stat_list = [diff_chr, dists, diff_stub, split_reads, dupe_reads, refs, zero_dists, num_reads, extrap_dup_rate, target_reads])
            bamfile (str): path to bam file
            outfile_prefix (str): path to where output files will be written
            count_diff_refname_stub (bool): whether we are counting the contig name stub differences.

        Returns:
            stat_dict ({str:float/str}: mappings of the stats and data suitable to be consumed by report generator
        '''

        print('Histograms written to:', self.paths['long_hist'], self.paths['short_hist'], self.paths['log_log_hist'])
        print('Duplicate saturation curve written to: {}'.format(self.paths['dup_sat_curve']))

        print('Number of contigs (more is harder):')
        print(self.out_stats['contigs'])

        print('Number of contigs greater than 10KB (longer contigs are better):')
        print(self.out_stats['contigs_greater_10k'])

        print('N50 of input assembly (longer contigs are better):')
        print(self.out_stats['N50'])

        print('Length of input assembly (bigger is harder):')
        print(self.out_stats['total_length'])

        print('Counts of zero distances (many is a sign of bad prep):')
        print(self.out_stats['zero_dist_pairs'],
              'of total',
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_zero_dist_pairs']
              )

        print('Count of same-contig read pairs with distance > 10KB (many is a sign of good prep):')
        print(self.out_stats['pairs_greater_10k'],
              'of total',
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_pairs_greater_10k']
              )

        print('Proportion of reads mapping to contigs > 10 Kbp with inserts > 10 Kbp:')
        print(self.out_stats['pairs_greater_10k_on_contigs_greater_10k'],
              'of total',
              self.out_stats['pairs_on_contigs_greater_10k'],
              self.out_stats['perc_pairs_greater_10k_on_contigs_greater_10k']
              )

        print('Count of read pairs with mates mapping to different chromosomes/contigs (sign of good prep IF same genome):')
        print(self.out_stats['intercontig_pairs'],
              'of total',
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_intercontig_pairs']
              )

        print('Count of split reads (more is usually good, as indicates presence of Hi-C junction in read):')
        print(self.out_stats['split_reads'],
              'of total',
              self.stats['total_reads'],
              self.out_stats['perc_split_reads']
              )

        print('Count of MAPQ zero reads (bad, ambiguously mapped):')
        print(self.out_stats['mapq0_reads'],
              'of total',
              self.out_stats['total_reads'],
              self.out_stats['perc_mapq0_reads']
              )

        print('Count of duplicate reads (-1 if insufficient to estimate; duplicates are bad; WILL ALWAYS BE ZERO UNLESS BAM FILE IS PREPROCESSED TO SET THE DUPLICATES FLAG):')
        print(self.out_stats['duplicate_reads'],
              'of total',
              self.out_stats['total_reads'],
              self.out_stats['perc_duplicate_reads']
              )

        print('Percent duplicated at {} reads: {} (-1 if insufficient to estimate)'.format(self.out_stats['target_read_total'], self.out_stats['extrapolated_dup_rate']))

        if count_diff_refname_stub:
            print('Count of read pairs with mates mapping to different reference groupings, e.g. genomes (sign of bad ' \
                  'prep potentially):')
            print(self.out_stats['different_ref_stub_pairs'],
                  'of total',
                  self.out_stats['total_read_pairs'],
                  self.out_stats['perc_different_ref_stub_pairs']
                  )

    def write_stat_table(self, outfile_prefix):
        '''Writes the stats as a plain text file.

        Args:
            stat_dict ({str:str/float}): dict mapping stat labels to their values and other info.
            outfile_name (str): a path to which to write the data.

        '''

        tsv = outfile_prefix + '.tsv'

        with open(tsv, "w") as outfile:
            for k, v in self.out_stats.items():
                if k == "refs" or k == "dists":
                    #skip long metadata fields we don't really need
                    continue
                print(k, v, sep="\t", file=outfile)

    def write_dists_file(self, outfile_prefix):
        '''Writes the dists as a plain text file.

        Args:
            stat_dict ({str:str/float}): dict mapping stat labels to their values and other info.
            outfile_name (str): a path to which to write the data.

        '''
        if not outfile_prefix.endswith(".dists"):
            outfile_name = outfile_prefix + ".dists"
        else:
            outfile_name = outfile_prefix

        with open(outfile_name, "w") as outfile:
            for k, v in self.dists.items():
                print(k, v, sep="\t", file=outfile)

    def write_pdf_report(self, outfile_prefix):
        '''Make the pdf report using the template.

        Args:
            qc_repo_path (str): path the qc tool repo, which holds markdown and style templates
            stat_dict ({str: str/float}): a dictionary of relevant data (stats, file paths) regarding the lib to be used in the report.
            outfile_prefix (str): path to which the report shall be written.

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

        template_path = os.path.join(self.paths['script_dir'], "collateral", "HiC_QC_report_template.md")
        style_path = os.path.join(self.paths['script_dir'], "collateral", "style.css")
        commit_path = os.path.join(self.paths['script_dir'], "collateral", "commit_id")
        if not os.path.exists(template_path):
            UserWarning("can't find markdown template at {0}! skipping making template.".format(
                qc_repo_path)
            )
            sys.exit(1)

        if os.path.exists(commit_path):
            with open(commit_path) as commit_file:
                commit_id = commit_file.read().strip()
        else:
            commit_id = "unversioned"


        with open(template_path) as template_fh:
            str = template_fh.read()
            sub_str = str.replace("COMMIT_VERSION", commit_id)  # versions report
            sub_str = sub_str.format(**self.out_stats)  # splat the statistics and path into the markdown, render as html
            html = md.markdown(sub_str, extensions=['tables', 'nl2br'])

            # write out just html
            with open(outfile_prefix + "_qc_report.html", 'w') as html_out:
                html_out.write(html)

            # print html
            pdfkit.from_string(html, outfile_prefix + "_qc_report.pdf", options=options, css=style_path)

def parse_args():
    '''parse command-line args

    Args:
        desc(str): program description, e.g. __file__

    Returns:
        args (dict): dict of the form {arg_name: arg_value}
    '''
    parser = argparse.ArgumentParser(description=__file__)
    parser.add_argument('-n', '--num_reads', default=1000000, type=int,
                        help='Number of reads from bam file to use. Use -1 to use all reads. Default: %(default)s')
    parser.add_argument('-b', '--bam_file', required=True, type=str,
                        help='BAM file to evaluate for QC')
    parser.add_argument('--count_diff_refname_stub', action='store_true',
                        help='Can be used to QC differently given reference name formatting. For internal use.')
    parser.add_argument('--outfile_prefix', '-o', default='Read_mate_dist', type=str,
                        help='Path to which to write plots to (PNG suffix will be attached).')
    parser.add_argument('--make_report', '-r', default=False, action='store_true',
                        help='Whether to export results in a PDF report. Requires that the QC script be' \
                             'in the same directory as the QC repo\'s collateral directory. Default: False.')
    parser.add_argument('--target_read_total', type=int, default=100000000, help='Total read count for duplicate read extrapolation (Default: %(default)s)')

    args = parser.parse_args()

    return args

def estimate_required_num_reads(diff_contig, refs, num_pairs, target=600.0):
    '''estimate the desired number of reads total based on the observed reads'''
    target_num = len(refs) * target  # desired among-contig
    diff_rate = float(diff_contig) / float(num_pairs)  # rate of among-contig per read pair
    total_num = target_num / diff_rate
    return int(total_num)

if __name__ == "__main__":
    args = parse_args()
    dirname = os.path.dirname(args.outfile_prefix)
    if dirname != '' and not os.path.exists(dirname):
        os.makedirs(dirname)

    QC = HiCQC()
    if args.num_reads != -1:
        print('parsing the first {} read pairs in bam file {} to QC Hi-C library quality'.format(args.num_reads, args.bam_file))
    else:
        print('parsing all read pairs in bam file {} to QC Hi-C library quality'.format(args.bam_file))

    QC.parse_bam(args.bam_file, max_read_pairs=args.num_reads)
    QC.plot_dup_saturation(args.outfile_prefix)
    QC.pass_judgement()
    QC.html_from_judgement()
    QC.plot_histograms(args.outfile_prefix)
    QC.stringify_stats()
    QC.print_stats()

    QC.write_stat_table(args.outfile_prefix)
    QC.write_dists_file(args.outfile_prefix)

    QC.write_pdf_report(args.outfile_prefix)
