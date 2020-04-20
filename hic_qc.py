#!/usr/bin/env python

# takes a bam file and makes a histogram of distances between mate alignments to the
# reference assembly
# takes the first 1M read pairs by default

### USAGE:
# python hic_qc.py -b <BAM_FILE> -n <NUM_READS_TO_USE> -o <outfile_stub>
# creates files in the working directory with relevant plots, also text files of statistics.
# flip -r flag  (assuming you have dependencies) to make a PDF report with everything together.

from __future__ import print_function
from __future__ import division

import sys
import pysam
import numpy as np
import argparse
import os
import logging
import matplotlib
from collections import Counter
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pdfkit
import markdown as md
from scipy import optimize
import re

from _version import get_versions
__version__ = get_versions()['version']

# default QC thresholds if there is no thresholds file
DEFAULT_MIN_SAME_STRAND_HQ_PERCENTAGE           =   0.015
DEFAULT_MIN_INFORMATIVE_READ_PAIRS_PERCENTAGE   =   0.05
DEFAULT_MAX_NONINFORMATIVE_READ_PAIR_PERCENTAGE =   0.50
DEFAULT_MIN_LONG_CONTACT_PERCENTAGE             =   0.03
DEFAULT_MIN_INTERCONTIG_CONTACT_PERCENTAGE      =   0.025
DEFAULT_MIN_USABLE_READS_PER_CONTIG             =   600
DEFAULT_MAX_DUPE_PERCENTAGE                     =   0.40
DEFAULT_MAX_ZERO_DIST_PERCENTAGE                =   0.20
DEFAULT_MAX_ZERO_MAPQ0_PERCENTAGE               =   0.20
DEFAULT_MAX_UNMAPPED_PERCENTAGE                 =   0.10

def saturation(x, V, K):
    '''Computes non-duplicate read count given x reads and parameters V and K.
    Intended for use within scipy optimize.

    Args:
        x (int): Read count
        V (float): Optimization parameter
        K (float): Optimization parameter
    Returns:
        (float) Estimated non-duplicate read count
    '''

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
    '''Class for extracting QC metrics from bam files and using them
    to create plots, print results, save tables, and create pdf reports.
    '''

    def __init__(self, outfile_prefix='Read_mate_dist', sample_type='genome', thresholds_file=None,
                 rp_stats=None, mq_stats=None, edist_stats=None, lib_enzyme=None):
        '''Initialize metrics for later extraction and conversion.
        '''
        logging.basicConfig(format="[%(name)s - %(asctime)s] %(message)s", level=logging.INFO)
        self.logger = logging.getLogger("hic_qc")

        self.sample_type = sample_type.lower()
        self.qc_purpose = 'Unknown'
        self.lib_enzyme = lib_enzyme if lib_enzyme is not None else ['undefined']

        if self.sample_type == 'metagenome':
            self.qc_purpose = 'Metagenome Deconvolution'
        elif self.sample_type == 'genome':
            self.qc_purpose = 'Genome Scaffolding'

        if thresholds_file is None:
            self.min_same_strand_hq_percentage              =   DEFAULT_MIN_SAME_STRAND_HQ_PERCENTAGE
            self.min_informative_read_pairs_percentage      =   DEFAULT_MIN_INFORMATIVE_READ_PAIRS_PERCENTAGE
            self.max_noninformative_read_pair_percentage    =   DEFAULT_MAX_NONINFORMATIVE_READ_PAIR_PERCENTAGE
            self.min_long_contact_percentage                =   DEFAULT_MIN_LONG_CONTACT_PERCENTAGE
            self.min_intercontig_contact_percentage         =   DEFAULT_MIN_INTERCONTIG_CONTACT_PERCENTAGE
            self.min_usable_reads_per_contig                =   DEFAULT_MIN_USABLE_READS_PER_CONTIG
            self.max_dupe_percentage                        =   DEFAULT_MAX_DUPE_PERCENTAGE
            self.max_zero_dist_percentage                   =   DEFAULT_MAX_ZERO_DIST_PERCENTAGE
            self.max_zero_mapq0_percentage                  =   DEFAULT_MAX_ZERO_MAPQ0_PERCENTAGE
            self.max_unmapped_percentage                    =   DEFAULT_MAX_UNMAPPED_PERCENTAGE
        else:
            import json
            with open(thresholds_file) as f:
                thresholds = json.load(f)
            self.min_same_strand_hq_percentage              =   float(thresholds[sample_type]['MIN_SAME_STRAND_HQ_PERCENTAGE'])
            self.min_informative_read_pairs_percentage      =   float(thresholds[sample_type]['MIN_INFORMATIVE_READ_PAIRS_PERCENTAGE'])
            self.max_noninformative_read_pair_percentage    =   float(thresholds[sample_type]['MAX_NONINFORMATIVE_READ_PAIR_PERCENTAGE'])
            self.min_long_contact_percentage                =   float(thresholds[sample_type]['MIN_LONG_CONTACT_PERCENTAGE'])
            self.min_intercontig_contact_percentage         =   float(thresholds[sample_type]['MIN_INTERCONTIG_CONTACT_PERCENTAGE'])
            self.min_usable_reads_per_contig                =   float(thresholds[sample_type]['MIN_USABLE_READS_PER_CONTIG'])
            self.max_dupe_percentage                        =   float(thresholds[sample_type]['MAX_DUPE_PERCENTAGE'])
            self.max_zero_dist_percentage                   =   float(thresholds[sample_type]['MAX_ZERO_DIST_PERCENTAGE'])
            self.max_zero_mapq0_percentage                  =   float(thresholds[sample_type]['MAX_ZERO_MAPQ0_PERCENTAGE'])
            self.max_unmapped_percentage                    =   float(thresholds[sample_type]['MAX_UNMAPPED_PERCENTAGE'])

        self.per_read_metrics = set(['duplicate_reads', 'mapq0_reads', 'split_reads', 'total_reads', 'unmapped_reads',])
        self.per_pair_metrics = set(['different_ref_stub_pairs',
                                     'informative_pairs',
                                     'intercontig_pairs',
                                     'intercontig_pairs_hq',
                                     'pairs_greater_10k',
                                     'pairs_greater_10k_on_contigs_greater_10k',
                                     'pairs_intracontig_hq',
                                     'pairs_intracontig_hq_gt10kbp',
                                     'pairs_on_contigs_greater_10k',
                                     'pairs_on_same_strand',
                                     'pairs_on_same_strand_hq',
                                     'proximo_usable_rp',
                                     'proximo_usable_rp_hq',
                                     'total_pairs_on_same_contig',
                                     'total_read_pairs',
                                     'total_read_pairs_hq',
                                     'zero_dist_pairs',
                                     ])
        # Dictionary of key --> numerator, denominator pairs for stringify_stats
        self.to_percents =    {
                               # Top Table
                               'perc_pairs_on_same_strand_hq': ('pairs_on_same_strand_hq', 'pairs_intracontig_hq'),
                               'perc_informative_read_pairs': ('informative_pairs', 'total_read_pairs'),
                               # Other Good Metrics
                               'perc_pairs_intra_hq_gt10kbp': ('pairs_intracontig_hq_gt10kbp', 'pairs_intracontig_hq'),
                               'perc_intercontig_pairs_hq_gt10kbp': ('pairs_intercontig_hq_gt10kbp', 'total_read_pairs_hq'),
                               # Usable reads number is in self.to_round
                               # Noninformative Read Breakdown
                               'perc_noninformative_read_pairs': ('noninformative_read_pairs', 'total_read_pairs'),
                               'perc_duplicate_reads': ('duplicate_reads', 'total_reads'),
                               'perc_zero_dist_pairs': ('zero_dist_pairs', 'total_read_pairs'),
                               'perc_unmapped_reads': ('unmapped_reads', 'total_reads'),
                               'perc_mapq0_reads': ('mapq0_reads', 'total_reads'),
                               # Extended Metrics
                               'perc_pairs_greater_10k': ('pairs_greater_10k', 'total_read_pairs'),
                               'perc_pairs_greater_10k_on_contigs_greater_10k': ('pairs_greater_10k_on_contigs_greater_10k', 'pairs_on_contigs_greater_10k'),
                               'perc_pairs_on_same_strand': ('pairs_on_same_strand', 'total_pairs_on_same_contig'),
                               'perc_intercontig_pairs': ('intercontig_pairs', 'total_read_pairs'),
                               'perc_intercontig_pairs_hq': ('intercontig_pairs_hq', 'total_read_pairs_hq'),
                               'perc_split_reads': ('split_reads', 'total_reads'),
                               'perc_hq_rp': ('total_read_pairs_hq', 'total_read_pairs'),
                               'perc_different_ref_stub_pairs': ('different_ref_stub_pairs', 'total_read_pairs')
                               }
        self.to_round = set([
                             'proximo_usable_rp_per_ctg_gt_5k',
                             'proximo_usable_rp_hq_per_ctg_gt_5k'
                             ])

        self.convert_to_pairs = set(['unmapped_reads', 'split_reads', 'duplicate_reads', 'mapq0_reads'])

        self.paths = {'script_dir': os.path.dirname(__file__), 'outfile_prefix': outfile_prefix}
        self.paths['pg_logo'] = os.path.join(self.paths['script_dir'], 'collateral', 'PGBlueLogoHorSmall.png')

        self.N50 = None
        self.stats = Counter()
        self.dists = Counter()
        self.total_array = []
        self.non_dup_array = []

        if rp_stats is not None and mq_stats is not None and edist_stats is not None:
            rps = []
            mqs = []
            edists = []
            self.mapping_dict = {}
            rp_stats = sorted(map(lambda x: int(1000 * x), rp_stats))
            mq_stats = sorted(map(int, mq_stats))
            edist_stats = list(reversed(sorted(map(int, edist_stats))))

            for rp in rp_stats:
                self.mapping_dict[rp] = {}
                for mq in mq_stats:
                    self.mapping_dict[rp][mq] = {}
                    for ed in edist_stats:
                        self.mapping_dict[rp][mq][ed] = 0
            self.rp_stats = rp_stats
            self.mq_stats = mq_stats
            self.edist_stats = edist_stats

        else:
            self.mapping_dict = None
            self.rp_array = None

    def parse_bam(self, bamfile, max_read_pairs=-1):
        '''Extract QC metrics from a specified bam file. It requires a read name sorted bam file.
        By default, it will parse all reads in the bam file, but a limit can be specified by max_read_pairs.
        This method performs read pairing.
        '''

        if max_read_pairs != -1:
            self.logger.info('parsing the first {} read pairs in bam file {} '\
                             'to QC Hi-C library quality'.format(max_read_pairs, bamfile))
        else:
            self.logger.info('parsing all read pairs in bam file {} '\
                             'to QC Hi-C library quality'.format(max_read_pairs, bamfile))

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
                    a = None
                    b = None
                else:
                    a = read

                if i % 1000 == 0:
                    self.update_dup_stats()
                i += 1

        self.finalize_stats()

    def extract_header_info(self, header, xx=50):
        '''Extract reference names, calculate N50, get total assembly length, and get set of contigs > 10kbp from a pysam header.
        Also checks if input bamfile is labeled as coordinate sorted and throws a ValueError if True.

            Args:
                header (pysam.AlignmentFile.header): the header from which to extract reference sequence lengths
                xx (float): the XX of NXX. we assume it is N50 so the default is 50

            Uses:
                self.mapping_dict (dict(int-->dict(int-->dict(int-->int)))): Nested dict with min_size, mapq, and edist keys
                                                                             and read pair counts as the inner values.

            Sets:
                self.refs dict(int-->str): ref_id --> ref_name mappings for the assembly
                self.contig_len (int): the NXX (probably N50) of the assembly based on the header
                self.total (int): the total length of the assembly
                self.contigs_greater_10k (set(str)): The set of names of contigs with length > 10Kbp
                self.contigs_greater_5k (set(str)): The set of names of contigs with length > 5Kbp
                self.contigs_greater (dict(int-->set(str))): Dictionary with minimum lengths as keys and sets of contigs as values
                self.command_line(str): Full command-line argument used for alignment
                self.bwa_command(str): Subset of self.command_line containing only the BWA options used
                self.samblaster(str): Command used by samblaster

            Raises:
                ValueError if header labels bamfile as coordinate sorted
        '''

        if 'HD' in header and 'SO' in header['HD'] and header['HD']['SO'].lower().strip() == 'coordinate':
            raise ValueError('Error: hic_qc.py requires read name sorted input, but bamfile {} is coordinate sorted'.format(self.paths['bamfile']))

        self.refs = header.references
        self.stats['num_refs'] = len(self.refs)

        self.N50, self.total_length = calc_nxx(header)
        self.contigs_greater_10k = set([contig['SN'] for contig in header['SQ'] if contig['LN'] > 10000])
        self.contigs_greater_5k = set([contig['SN'] for contig in header['SQ'] if contig['LN'] > 5000])

        self.contigs_greater = {}

        if self.mapping_dict is not None:
            for min_size in self.mapping_dict.keys():
                self.contigs_greater[min_size] = set([contig['SN'] for contig in header['SQ'] if contig['LN'] > min_size])

        # TODO: add more robust logic to different BAM headers, and/or comment the assumptions made in this code
        if 'PG' in header and 'bwa' in header['PG'][0]['CL']:
            self.bwa_command_line = header['PG'][0]['CL']
            self.bwa_command = re.search(r'(bwa )[^//]*', self.bwa_command_line).group()
        else:
            self.bwa_command = 'BWA command not found'

        if 'PG' in header and 'samblaster ' in header['PG'][1]['CL']:
            self.samblaster = header['PG'][1]['CL']
        else:
             self.samblaster = 'samblaster command not found'

        if 'PG' in header:
            self.ref_assembly_path = re.search(r'(?<= /).*(?!\.fastq)(\.fasta|\.fna|\.fa)', header['PG'][0][
                'CL']).group(0)
            self.ref_assembly = os.path.basename(self.ref_assembly_path.strip())
        else:
            self.ref_assembly = "reference assembly not found"

    def process_pair(self, a, b):
        '''Extract stats from a pair of reads.
        Updates per read and per mapped pair stats separately.

        Args:
            a (pysam.AlignedSegment) One read in pair
            b (pysam.AlignedSegment) Other read in pair
        '''

        self.stats['total_read_pairs'] += 1

        self.update_read_stats(a)
        self.update_read_stats(b)

        if self.is_noninformative_read_pair(a, b):
            self.stats['noninformative_read_pairs'] += 1
        elif self.is_informative_pair(a, b):
            self.stats['informative_pairs'] += 1

        if not a.is_unmapped and not b.is_unmapped:
            self.update_mapped_pair_stats(a, b)

    def update_read_stats(self, read):
        '''Update per read stats based on given read.

        Args:
            read (pysam.AlignedSegment): read to extract stats from
        '''

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
        '''Update mapped pair stats given a pair of reads.

        Args:
            a (pysam.AlignedSegment): One read
            b (pysam.AlignedSegment): Other read
        '''

        is_high_qual_pair = self.is_high_qual_pair(a, b)

        if is_high_qual_pair:
            self.stats['total_read_pairs_hq'] += 1

        if a.reference_name != b.reference_name:
            self.stats['intercontig_pairs'] += 1

            if is_high_qual_pair:
                self.stats['intercontig_pairs_hq'] += 1
                if a.reference_name in self.contigs_greater_10k and b.reference_name in self.contigs_greater_10k:
                    self.stats['pairs_intercontig_hq_gt10kbp'] += 1

            refa_stub = a.reference_name.split('.')[0]
            refb_stub = b.reference_name.split('.')[0]

            if refa_stub != refb_stub:
                self.stats['different_ref_stub_pairs'] += 1

            if a.reference_name in self.contigs_greater_5k and b.reference_name in self.contigs_greater_5k:
                if min(a.mapping_quality, b.mapping_quality) > 0 and not any([a.is_duplicate, b.is_duplicate]):
                    self.stats['proximo_usable_rp'] += 1
                if is_high_qual_pair:
                    self.stats['proximo_usable_rp_hq'] += 1

        else:
            self.stats['total_pairs_on_same_contig'] += 1
            if (a.is_reverse and b.is_reverse) or (not a.is_reverse and not b.is_reverse):
                self.stats['pairs_on_same_strand'] += 1
            dist = abs(a.reference_start - b.reference_start)
            self.dists[dist] += 1

            if is_high_qual_pair:
                self.stats['pairs_intracontig_hq'] += 1
                if dist > 10000:
                    self.stats['pairs_intracontig_hq_gt10kbp'] += 1
                if (a.is_reverse and b.is_reverse) or (not a.is_reverse and not b.is_reverse):
                    self.stats['pairs_on_same_strand_hq'] += 1

            if dist > 10000:
                self.stats['pairs_greater_10k'] += 1
            if dist == 0:
                self.stats['zero_dist_pairs'] += 1

            ## Additional stats for dynamo table.
            if dist < 0:
                raise ValueError('Error: Distance between reads is less than zero.')
            elif dist < 1e3:
                self.stats['reads_spanning_up_to_1k'] += 1
            elif dist < 1e4:
                self.stats['reads_spanning_1k_to_10k'] += 1
            elif dist < 1e5:
                self.stats['reads_spanning_10k_to_100k'] += 1
            elif dist < 1e6:
                self.stats['reads_spanning_100k_to_1000k'] += 1
            else:
                self.stats['reads_spanning_greater_than_1000k'] += 1

            if a.reference_name in self.contigs_greater_10k:
                self.stats['pairs_on_contigs_greater_10k'] += 1
                if is_high_qual_pair:
                    self.stats['pairs_on_contigs_greater_10k_hq'] += 1
                    if dist > 10000:
                        self.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] += 1
                if dist > 10000:
                    self.stats['pairs_greater_10k_on_contigs_greater_10k'] += 1

            if self.mapping_dict is not None and \
               not (a.is_duplicate or b.is_duplicate):
                self.update_rp_array(a, b)

    def update_rp_array(self, a, b):
        '''Update nested dict of mapping stats based on current read pair a, b.

        Uses:
            self.rp_stats (list(int)): List of minimum insert sizes sorted from low to high
            self.mq_stats (list(int)): List of minimum mapq values sorted from low to high
            self.edist_stats (list(int)): List of maximum edit distances sorted from high to low

        Modifies:
            self.mapping_dict (dict(int-->dict(int-->dict(int-->int)))): Nested dict with min_size, mapq, and edist keys
                                                                         and read pair counts as the inner values.
        '''
        mq = min(a.mapping_quality, b.mapping_quality)
        isize = abs(a.reference_start - b.reference_start)
        edist = max(a.get_tag('NM'), b.get_tag('NM'))

        for min_size in self.rp_stats:
            if isize >= min_size:
                for mapq in self.mq_stats:
                    if mq >= mapq:
                        for ed in self.edist_stats:
                            if edist <= ed:
                                self.mapping_dict[min_size][mapq][ed] += 1

    def is_noninformative_read_pair(self, a, b):
        is_noninformative = False
        if a.is_unmapped or b.is_unmapped:
            is_noninformative = True
        elif a.is_duplicate or b.is_duplicate:
            is_noninformative = True
        elif a.mapping_quality == 0 or b.mapping_quality == 0:
            is_noninformative = True
        elif a.reference_name == b.reference_name:
            if abs(a.reference_start - b.reference_start) == 0:
                is_noninformative = True
        return is_noninformative

    def is_informative_pair(self, a, b):
        is_informative = True
        if a.is_unmapped or b.is_unmapped:
            is_informative = False
        elif a.is_duplicate or b.is_duplicate:
            is_informative = False
        elif a.mapping_quality == 0 or b.mapping_quality == 0:
            is_informative = False
        elif a.reference_name == b.reference_name and \
             abs(a.reference_start - b.reference_start) < 10000:
            is_informative = False
        return is_informative

    def is_high_qual_pair(self, a, b):
        return min(a.mapping_quality, b.mapping_quality) >= 20 and \
               max(a.get_tag('NM'), b.get_tag('NM') <= 5) and \
               not a.is_duplicate and not b.is_duplicate

    def update_dup_stats(self):
        '''Update lists of duplication statistics.
        '''

        self.total_array.append(self.stats['total_reads'])
        self.non_dup_array.append(self.stats['total_reads'] - self.stats['duplicate_reads'])

    def finalize_stats(self):
        '''Finalize stats from a bam file after reads are processed.
        '''
        self.total_array = np.array(self.total_array)
        self.non_dup_array = np.array(self.non_dup_array)
        if self.stats['pairs_on_contigs_greater_10k'] > 0:
            self.stats['proportion_pairs_greater_10k_on_contigs_greater_10k'] = self.stats['pairs_greater_10k_on_contigs_greater_10k'] / \
                                                                            self.stats['pairs_on_contigs_greater_10k']
        else:
            self.stats['proportion_pairs_greater_10k_on_contigs_greater_10k'] = 0

        if len(self.contigs_greater_5k) > 0:
            self.stats['proximo_usable_rp_per_ctg_gt_5k'] = self.stats['proximo_usable_rp'] / len(self.contigs_greater_5k)
            self.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = self.stats['proximo_usable_rp_hq'] / len(self.contigs_greater_5k)
        else:
            self.stats['proximo_usable_rp_per_ctg_gt_5k'] = 0
            self.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 0

        if self.mapping_dict is not None:
            self.write_mapping_stats()

        # We are stricter on wanting a low number of dupes when it looks like we are only looking at a QC amount of sequencing (<10M read pairs)
        self.allowed_dupe_percentage = 1.0 if self.stats['total_read_pairs'] > 1e7 else 0.5

    def write_mapping_stats(self):
        with open('{}.mapping_stats.tsv'.format(self.paths['outfile_prefix']), 'w') as outfile:
            print('edist', 'mapq', 'min_size', 'count', sep='\t', file=outfile)
            for min_size in self.rp_stats:
                for mapq in self.mq_stats:
                    for ed in self.edist_stats:
                        count = self.mapping_dict[min_size][mapq][ed]
                        print(ed, mapq, min_size, count, sep='\t', file=outfile)

        if self.mapping_dict is not None:
            self.write_mapping_stats()

    def write_mapping_stats(self):
        with open('{}.mapping_stats.tsv'.format(self.paths['outfile_prefix']), 'w') as outfile:
            print('edist', 'mapq', 'min_size', 'count', sep='\t', file=outfile)
            for min_size in self.rp_stats:
                for mapq in self.mq_stats:
                    for ed in self.edist_stats:
                        count = self.mapping_dict[min_size][mapq][ed]
                        print(ed, mapq, min_size, count, sep='\t', file=outfile)

    def plot_dup_saturation(self, target_x=100000000, min_sample=100000, target_y=None):
        '''Fit and plot a saturation curve from cumulative total and non-dup read counts.

            Args:
                target_x (int): Total reads to extrapolate to.
                target_y (int): Non-dup reads to extrapolate to. If specified with target_x, plots a point.

            Uses:
                self.paths['outfile_prefix'] (str): Path prefix for output files.
                self.total_array (np.array(int)): Numpy array of total reads, recorded every 1000 reads.
                self.non_dup_array (np.array(int)): Numpy array of non-duplicate read counts, recorded every 1000 reads.

            Sets:
                self.stats['observed_dup_rate'] (float): Observed rate of read duplication
                self.stats['extrapolated_dup_rate'] (float): Rate of duplication extrapolated to target_read_total
                self.stats['target_read_total'] (int): Read count to extrapolate to.
                self.stats['dup_sat_V'] (float): Optimization parameter for saturation curve.
                self.stats['dup_sat_K'] (float): Optimization parameter for saturation curve.
                self.paths['dup_sat_curve'] (str): Path to dup_sat_curve plot.
        '''

        self.stats['observed_dup_rate'] = -1
        self.stats['extrapolated_dup_rate'] = -1
        self.stats['target_read_total'] = target_x
        self.stats['dup_sat_V'] = -1
        self.stats['dup_sat_K'] = -1

        outfile = self.paths['outfile_prefix'] + '.dup_saturation.png'
        self.paths['dup_sat_curve'] = outfile

        if self.total_array.size == 0 or self.total_array[-1] < min_sample:
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
            self.logger.info('At {} reads, estimated {:.0f} non-dup reads'.format(target_x, saturation(target_x, *params)))
            self.logger.info('True non-dup reads: {}'.format(target_y))
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
                self.paths['bamname'], self.stats['observed_dup_rate'], 1-non_dup_rate))
        else:
            plt.title('{}\nproportion duplicated (sampled): {:.2f}'.format(
                self.paths['bamname'], self.stats['observed_dup_rate']))
        plt.xlabel('Total reads')
        plt.ylabel('Non-duplicate reads')
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()

        self.logger.info('Best V = {}, best K = {}'.format(self.stats['dup_sat_V'], self.stats['dup_sat_K']))

        return 0

    def plot_histograms(self):
        '''Make the read distance long, short, and log_log histograms using matplotlib and write them to disk.

        Args:
            self.dists (dict(int, int) of mate distances and counts): Distances to plot in histogram.
            self.num_pairs (int): number of read pairs analyzed

        Uses:
            self.paths['outfile_prefix'] (str):  Path prefix for output files.
        '''

        num_dists = sum(self.dists.values())
        num_pairs = self.stats['total_read_pairs']
        title_string = '\nMate distance distribution for first {} read pairs for sample\n{}'.format(num_pairs,
                                                                                                    self.paths['bamname'])
        key_len = len(self.dists)

        long_hist_path = self.paths['outfile_prefix'] + '_long.png'
        fig1, ax = plt.subplots(1)
        if key_len > 0:
            plt.hist(list(self.dists.keys()), weights=list(self.dists.values()), bins=50, edgecolor='black', color='red')

            ax.set_ylim(0.5, max(num_dists * 2, 1))
            plt.yscale('log', nonposy='clip')
            plt.title(title_string)
            plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig)')
            plt.ylabel('Number of reads')
        else:
            plt.title('Warning: No read pair distribution to plot')
        fig1.savefig(long_hist_path)
        plt.close(fig1)
        self.paths['long_hist'] = long_hist_path

        fig2, ax = plt.subplots(1)
        short_hist_path = self.paths['outfile_prefix'] + '_short.png'

        if key_len > 0:
            plt.hist(list(self.dists.keys()), weights=list(self.dists.values()), bins=range(0, 20000, 500), edgecolor='black', color='red')
            ax.set_xlim(0, 20000)
            ax.set_ylim(0.5, num_pairs * 2)
            plt.yscale('log', nonposy='clip')
            plt.title(title_string)
            plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig)')
            plt.ylabel('Number of reads')
        else:
            plt.title('Warning: No read pair distribution to plot')
        fig2.savefig(short_hist_path)
        plt.close(fig2)
        self.paths['short_hist'] = short_hist_path

        fig3, ax = plt.subplots(1)
        log_log_hist_path = self.paths['outfile_prefix'] + '_log_log.png'

        offset_dists = {}
        for key, value in self.dists.items():
            offset_dists[key+1] = value

        if key_len > 0:
            min_dist = min(self.dists.keys())
            max_dist = max(self.dists.keys())

            plt.hist(list(offset_dists.keys()),
                     weights=list(offset_dists.values()),
                     bins=np.logspace(np.log10(min_dist+1),
                                      np.log10(max_dist),
                                      50),
                     log=True, edgecolor='black', color='red')
            ax.set_ylim(0.5, max(num_pairs * 2, 1))
            plt.yscale('log', nonposy='clip')
            plt.xscale('log')
            plt.xlim(left=1)
            plt.title(title_string)
            plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig, log scale)')
            plt.ylabel('Number of reads (log scale)')
            plt.tight_layout()
        else:
            plt.title('Warning: No read pair distribution to plot')

        fig3.savefig(log_log_hist_path)
        plt.close(fig3)

        self.paths['log_log_hist'] = log_log_hist_path

        fig4, ax = plt.subplots(1)
        log_log_norm_hist_path = self.paths['outfile_prefix'] + '_log_log_norm.png'

        offset_dists = {}
        for key, value in self.dists.items():
            offset_dists[key+1] = value

        if key_len > 0:
            min_dist = min(self.dists.keys())
            max_dist = max(self.dists.keys())

            plt.hist(list(offset_dists.keys()),
                     weights=list(offset_dists.values()),
                     bins=np.logspace(np.log10(min_dist+1),
                                      np.log10(max_dist),
                                      50),
                     log=True, edgecolor='black', color='red', density=True)
            ax.set_ylim(0.0001, 1.0)
            plt.yscale('log', nonposy='clip')
            plt.xscale('log')
            plt.xlim(left=1)
            plt.title(title_string)
            plt.xlabel('Distance between read pair mates in Hi-C mapping (same contig, log scale)')
            plt.ylabel('Density of reads (density, log scale)')
            plt.tight_layout()
        else:
            plt.title('Warning: No read pair distribution to plot')

        fig4.savefig(log_log_norm_hist_path)
        plt.close(fig4)

        self.paths['log_log_norm_hist'] = log_log_norm_hist_path

    def html_from_judgement(self):
        '''Set a formatted HTML string based on two judgment bool values
            Uses:
                self.judge_good (bool): does the hi-c library show characteristics of 'goodness', e.g. many long-distance contacts etc.
                self.judge_bad (bool): does the hi-c library show 'bad' characteristics, e.g. zero-distance reads or too many duplicates.
            Sets:
                self.judge_html (str): an HTML string to be substituted into the report to subjectively grade the assembly. 4 possibilities.
            Raises:
                ValueError: if impossible logical situations occur given two bools.

        '''
        if self.judge_good and not self.judge_bad:
            self.judge_html = '<span class="pass">SUFFICIENT</span>'
        elif not self.judge_good and self.judge_bad:
            self.judge_html = '<span class="fail">INSUFFICIENT</span>'
        elif self.judge_good and self.judge_bad:
            self.judge_html = '<span class="mixed-results">MIXED RESULTS</span>'
        elif not self.judge_good and not self.judge_bad:
            self.judge_html = '<span class="low-signal">LOW SIGNAL</span>'
        else:
            raise ValueError('logical impossibility!')

        # driving metrics
        if self.good_same_strand:
            self.same_strand_hq_html = '<span class="pass">{0}</span>'
        else:
            self.same_strand_hq_html = '<span class="fail">{0}</span>'

        if self.good_informative_read_pairs:
            self.informative_read_pairs_html = '<span class="pass">{0}</span>'
        else:
            self.informative_read_pairs_html = '<span class="fail">{0}</span>'

        # other good metrics
        if self.good_long_contacts:
            self.long_contacts_html = '<span class="pass">{0}</span>'
        else:
            self.long_contacts_html = '<span class="fail">{0}</span>'

        if self.good_intercontig_contacts:
            self.intercontig_hq_contacts_html = '<span class="pass">{0}</span>'
        else:
            self.intercontig_hq_contacts_html = '<span class="fail">{0}</span>'

        if self.good_usable_reads:
            self.usable_hq_gt_5k_html = '<span class="pass">{0}</span>'
        else:
            self.usable_hq_gt_5k_html = '<span class="fail">{0}</span>'

        # noninformative reads breakdown
        if not self.bad_noninformative_read_pairs:
            self.noninformative_read_pairs_html = '<span class="pass">{0}</span>'
        else:
            self.noninformative_read_pairs_html = '<span class="fail">{0}</span>'

        if not self.high_dupe:
            self.high_dupe_html = '<span class="pass">{0}</span>'
        else:
            self.high_dupe_html = '<span class="fail">{0}</span>'

        if not self.many_zero_dist_pairs:
            self.many_zero_dist_pairs_html = '<span class="pass">{0}</span>'
        else:
            self.many_zero_dist_pairs_html = '<span class="fail">{0}</span>'

        if not self.many_unmapped_reads:
            self.many_unmapped_reads_html = '<span class="pass">{0}</span>'
        else:
            self.many_unmapped_reads_html = '<span class="fail">{0}</span>'

        if not self.many_mapq_zero_reads:
            self.many_zero_mapq_reads_html = '<span class="pass">{0}</span>'
        else:
            self.many_zero_mapq_reads_html = '<span class="fail">{0}</span>'

    def pass_judgement(self):
        '''Pass judgement on the library according to certain mostly subjective ideas about what is good

        Uses:
            self.stats ({str: float/str}): mapping of lib characteristics to their values

        Sets:
            self.judge_good (bool): does the hi-c library show characteristics of 'goodness', e.g. many same strand contacts.
            self.judge_bad (bool): does the hi-c library show 'bad' characteristics, e.g. zero-mapq reads or too many duplicates.
            self.judge_html (str): an HTML string to put into pass/fail box
        '''
        # these metrics drive the subjective quality judgement
        self.good_same_strand              = float(self.stats['pairs_on_same_strand_hq']) / \
                                             max(self.stats['pairs_intracontig_hq'], 1) > \
                                             self.min_same_strand_hq_percentage
        self.good_informative_read_pairs   = float(self.stats['informative_pairs']) / \
                                             max(self.stats['total_read_pairs'], 1) > \
                                             self.min_informative_read_pairs_percentage

        # other good metrics
        self.good_long_contacts        = float(self.stats['pairs_greater_10k_on_contigs_greater_10k_hq']) / \
                                         max(self.stats['pairs_on_contigs_greater_10k_hq'], 1) > \
                                         self.min_long_contact_percentage
        self.good_intercontig_contacts = float(self.stats['pairs_intercontig_hq_gt10kbp']) / \
                                         max(self.stats['total_read_pairs_hq'], 1) > \
                                         self.min_intercontig_contact_percentage
        self.good_usable_reads         = float(self.stats['proximo_usable_rp_hq_per_ctg_gt_5k']) > \
                                         self.min_usable_reads_per_contig

        # noninformative read breakdown
        # We are stricter on wanting a low number of dupes when it looks like we are only looking at a QC amount of sequencing (<10M read pairs)
        self.bad_noninformative_read_pairs = float(self.stats['noninformative_read_pairs']) / \
                                             max(self.stats['total_read_pairs'], 1) > \
                                             self.max_noninformative_read_pair_percentage
        self.high_dupe = float(self.stats['duplicate_reads']) / \
                         max(self.stats['total_reads'], 1) > \
                         self.max_dupe_percentage * \
                         self.allowed_dupe_percentage
        self.many_zero_dist_pairs = float(self.stats['zero_dist_pairs']) / \
                                    max(self.stats['total_read_pairs'], 1) > \
                                    self.max_zero_dist_percentage
        self.many_unmapped_reads  = float(self.stats['unmapped_reads']) / \
                                    max(self.stats['total_reads'], 1) > \
                                    self.max_unmapped_percentage
        self.many_mapq_zero_reads = float(self.stats['mapq0_reads']) / \
                                    max(self.stats['total_reads'], 1) > \
                                    self.max_zero_mapq0_percentage

        self.judge_good = self.good_same_strand
        self.judge_bad  = not self.good_informative_read_pairs

        self.html_from_judgement()

    def stringify_stats(self):
        '''Convert stats to output dictionary with pretty strings and percents.

        Uses:
            self.to_percents ({str: (int, int)}): Mapping of keys to numerator, denominator pair for conversion to percents.
            self.to_round (set(str)): Set of keys from stats dict to round.
            self.convert_to_pairs (set(str)): Set of keys from stats dict that represent per read statistics.
            self.per_pair_metrics (set(str)): Set of keys from stats dict that represent per read pair statistics.
            self.paths ({str: str}): Mapping of names to paths.

        Sets:
            self.other_stats ({str: (float, str)}): Mapping of keys to value, format pair.
            self.out_stats ({str: str}): Mapping of stat keys to formatted strings.
        '''


        if self.stats['extrapolated_dup_rate'] > 0:
            extrap_dup_rate = self.stats['extrapolated_dup_rate'] * 100
        else:
            extrap_dup_rate = self.stats['extrapolated_dup_rate']

        # Dict of key --> (value, fmt) pairs for items that aren't counts
        self.other_stats = {
                            'N50': (self.N50, '{:,}'),
                            'contigs': (len(self.refs), '{:,}'),
                            'contigs_greater_10k': (len(self.contigs_greater_10k), '{:,}'),
                            'contigs_greater_5k': (len(self.contigs_greater_5k), '{:,}'),
                            'total_length': (self.total_length, '{:,}'),
                            'total_reads': (self.stats['total_reads'], '{:,}'),
                            'target_read_total': (self.stats['target_read_total'], '{:,}'),
                            'extrapolated_dup_rate': (extrap_dup_rate, '{:.2f}%'),
                            'judgment': (self.judge_html, '{}'),
                            'qc_purpose': (self.qc_purpose, '{}'),
                            'same_strand_threshold': (100.0 * self.min_same_strand_hq_percentage, '{}'),
                            'informative_read_pairs_threshold': (100.0 * self.min_informative_read_pairs_percentage, '{}'),
                            'noninformative_read_pairs_threshold': (100.0 * self.max_noninformative_read_pair_percentage, '{}'),
                            'long_contacts_threshold': (100.0 * self.min_long_contact_percentage, '{}'),
                            'intercontig_hq_contacts_threshold': (100.0 * self.min_intercontig_contact_percentage, '{}'),
                            'usable_hq_gt_5k_threshold': (self.min_usable_reads_per_contig, '{}'),
                            'high_dupe_threshold': (100.0 * self.max_dupe_percentage * self.allowed_dupe_percentage, '{}'),
                            'many_zero_dist_threshold': (100.0 * self.max_zero_dist_percentage, '{}'),
                            'many_zero_mapq_threshold': (100.0 * self.max_zero_mapq0_percentage, '{}'),
                            'many_unmapped_threshold': (100.0 * self.max_unmapped_percentage, '{}'),
                            'alignment_command_line': (self.bwa_command, '{}'),
                            'samblaster': (self.samblaster, '{}'),
                            'lib_enzyme': (', '.join(self.lib_enzyme), '{}'),
                            'ref_assembly': (self.ref_assembly, '{}'),
                            }
        self.out_stats = {}
        for key, (num, denom) in self.to_percents.items():
            try:
                self.out_stats[key] = '{:.2f}%'.format((self.stats[num] / self.stats[denom]) * 100)
            except ZeroDivisionError as e:
                self.out_stats[key] = 'NaN'

        for key in self.to_round:
            self.out_stats[key] = '{:.2f}'.format(self.stats[key])

        for item in self.convert_to_pairs:
            self.out_stats[item] = '{:,}'.format(self.stats[item] // 2)

        for item in self.per_pair_metrics:
            self.out_stats[item] = '{:,}'.format(self.stats[item])

        for item in self.paths:
            if item == "bamname":
                self.out_stats[item] = self.paths[item]
            else:
                self.out_stats[item] = os.path.abspath(self.paths[item])

        for key, value in self.contigs_greater.items():
            key_str = 'contigs_greater_{:.0f}k'.format(key / 1000)
            self.out_stats[key_str] = '{:,}'.format(len(value))

        for key, (value, fmt) in self.other_stats.items():
            self.out_stats[key] = fmt.format(value)

        # driving stats
        self.out_stats['same_strand_hq_html'] = self.same_strand_hq_html.format(self.out_stats['perc_pairs_on_same_strand_hq'])
        self.out_stats['informative_read_pairs_html'] = self.informative_read_pairs_html.format(self.out_stats['perc_informative_read_pairs'])
        self.out_stats['noninformative_read_pairs_html'] = self.noninformative_read_pairs_html.format(self.out_stats['perc_noninformative_read_pairs'])
        # other good metrics
        self.out_stats['long_contacts_html'] = self.long_contacts_html.format(self.out_stats['perc_pairs_intra_hq_gt10kbp'])
        self.out_stats['intercontig_hq_contacts_html'] = self.intercontig_hq_contacts_html.format(self.out_stats['perc_intercontig_pairs_hq_gt10kbp'])
        self.out_stats['usable_hq_gt_5k_html'] = self.usable_hq_gt_5k_html.format(self.out_stats['proximo_usable_rp_hq_per_ctg_gt_5k'])
        # noninformative breakdown
        self.out_stats['high_dupe_html'] = self.high_dupe_html.format(self.out_stats['perc_duplicate_reads'])
        self.out_stats['many_zero_dist_pairs_html'] = self.many_zero_dist_pairs_html.format(self.out_stats['perc_zero_dist_pairs'])
        self.out_stats['many_unmapped_reads_html'] = self.many_unmapped_reads_html.format(self.out_stats['perc_unmapped_reads'])
        self.out_stats['many_zero_mapq_reads_html'] = self.many_zero_mapq_reads_html.format(self.out_stats['perc_mapq0_reads'])

        self.out_stats['version'] = __version__

    def log_stats(self, count_diff_refname_stub=False):
        '''Log statistical summary.

        Uses:
            self.paths ({str: str}): Mapping of names to paths.
            self.out_stats ({str: str}): Mapping of stat keys to formatted strings.
            count_diff_refname_stub (bool): Whether we are counting the contig name stub differences.
        '''

        ('Histograms written to:', self.paths['long_hist'], self.paths['short_hist'], self.paths['log_log_hist'])
        self.logger.info('Duplicate saturation curve written to: {}'.format(self.paths['dup_sat_curve']))

        self.logger.info('Number of contigs (more is harder):')
        self.logger.info(self.out_stats['contigs'])

        self.logger.info('Number of contigs greater than 10KB (longer contigs are better):')
        self.logger.info(self.out_stats['contigs_greater_10k'])

        self.logger.info('N50 of input assembly (longer contigs are better):')
        self.logger.info(self.out_stats['N50'])

        self.logger.info('Length of input assembly (bigger is harder):')
        self.logger.info(self.out_stats['total_length'])

        self.logger.info('Counts of zero distances (many is a sign of bad prep):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['zero_dist_pairs'],
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_zero_dist_pairs'])
              )

        self.logger.info('Count of same-contig read pairs with distance > 10KB (many is a sign of good prep):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['pairs_greater_10k'],
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_pairs_greater_10k'])
              )

        self.logger.info('Proportion of reads mapping to contigs > 10 Kbp with inserts > 10 Kbp:')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['pairs_greater_10k_on_contigs_greater_10k'],
              self.out_stats['pairs_on_contigs_greater_10k'],
              self.out_stats['perc_pairs_greater_10k_on_contigs_greater_10k'])
              )

        self.logger.info('Count of read pairs with mates mapping to different chromosomes/contigs ' \
                         '(sign of good prep IF same genome):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['intercontig_pairs'],
              self.out_stats['total_read_pairs'],
              self.out_stats['perc_intercontig_pairs'])
              )

        self.logger.info('Count of split reads (more is usually good, as indicates presence of Hi-C junction in read):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['split_reads'],
              self.stats['total_reads'],
              self.out_stats['perc_split_reads'])
              )

        self.logger.info('Count of MAPQ zero reads (bad, ambiguously mapped):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['mapq0_reads'],
              self.out_stats['total_reads'],
              self.out_stats['perc_mapq0_reads'])
              )

        self.logger.info('Count of duplicate reads (-1 if insufficient to estimate; duplicates are bad; ' \
                         'WILL ALWAYS BE ZERO UNLESS BAM FILE IS PREPROCESSED TO SET THE DUPLICATES FLAG):')
        self.logger.info('{} of total {} {}%'.format(
              self.out_stats['duplicate_reads'],
              self.out_stats['total_reads'],
              self.out_stats['perc_duplicate_reads'])
              )

        self.logger.info('Percent duplicated at {} reads: {} ' \
                         '(-1 if insufficient to estimate)'.format(self.out_stats['target_read_total'],
                                                                   self.out_stats['extrapolated_dup_rate']))

        if count_diff_refname_stub:
            self.logger.info('Count of read pairs with mates mapping to different reference groupings, ' \
                             'e.g. genomes (sign of bad prep potentially):')
            self.logger.info('{} of total {} {}%'.format(
                  self.out_stats['different_ref_stub_pairs'],
                  self.out_stats['total_read_pairs'],
                  self.out_stats['perc_different_ref_stub_pairs'])
                  )

    def write_stat_table(self):
        '''Writes the stats as a plain text file.

        Uses:
            self.out_stats ({str: str}): Mapping of stat keys to formatted strings.
            self.paths['outfile_prefix'] (str):  Path prefix for output files.
        '''

        tsv = self.paths['outfile_prefix'] + '.tsv'

        with open(tsv, "w") as outfile:
            for k, v in self.out_stats.items():
                if k == "refs" or k == "dists":
                    #skip long metadata fields we don't really need
                    continue
                print(k, v, sep="\t", file=outfile)

    def write_dists_file(self):
        '''Writes the dists as a plain text file.

        Uses:
            self.out_stats ({str: str}): Mapping of stat keys to formatted strings.
            self.paths['outfile_prefix'] (str):  Path prefix for output files.

        Args:
            stat_dict ({str:str/float}): dict mapping stat labels to their values and other info.
        '''
        if not self.paths['outfile_prefix'].endswith(".dists"):
            outfile_name = self.paths['outfile_prefix'] + ".dists"
        else:
            outfile_name = self.paths['outfile_prefix']

        with open(outfile_name, "w") as outfile:
            for k, v in self.dists.items():
                print(k, v, sep="\t", file=outfile)

    def write_pdf_report(self, quiet=False):
        '''Make the pdf report using the template.
        Requires markdown template in the collateral directory of hic_qc.

        Uses:
            self.paths['script_dir']: Path to hic_qc directory.
            self.out_stats ({str: str}): Mapping of stat keys to formatted strings.
            self.paths['outfile_prefix'] (str):  Path prefix for output files.
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

        if quiet:
            options["quiet"] = ""

        template_path = os.path.join(self.paths['script_dir'], "collateral", "HiC_QC_report_template.md")
        style_path = os.path.join(self.paths['script_dir'], "collateral", "style.css")

        if not os.path.exists(template_path):
            raise FileNotFoundError("Can't find markdown template at {}! Exitting...".format(
                template_path)
            )

        with open(template_path) as template_fh:
            template_string = template_fh.read()
            sub_str = template_string.format(**self.out_stats)  # splat the statistics and path into the markdown, render as html
            html = md.markdown(sub_str, extensions=['markdown.extensions.tables', 'markdown.extensions.nl2br'])

            # write out just html
            with open(self.paths['outfile_prefix'] + "_qc_report.html", 'w') as html_out:
                html_out.write(html)

            # print html
            pdfkit.from_string(html, self.paths['outfile_prefix'] + "_qc_report.pdf", options=options, css=style_path)

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
    parser.add_argument('--target_read_total', type=int, default=100000000,
                        help='Total read count for duplicate read extrapolation (Default: %(default)s)')
    parser.add_argument('--rp_stats', nargs='+', default=[0, 1, 2, 5, 10, 20, 50],
                        help='List of distances in Kbp to calculate RP stats for (Default: %(default)s)')
    parser.add_argument('--mq_stats', nargs='+', default=[0, 1, 10, 20, 30, 40],
                        help='List of min MQ scores to calculate RP stats for (Default: %(default)s)')
    parser.add_argument('--edist_stats', nargs='+', default=[100, 10, 5, 3, 1, 0],
                        help='List of max edist scores to calculate RP stats for (Default: %(default)s)')
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--thresholds', default='{0}/collateral/thresholds.json'.format(os.path.dirname(os.path.realpath(__file__))),
                        help='JSON file containing QC thresholds (Default: %(default)s)')
    parser.add_argument('--sample_type', default='genome', choices=['genome', 'metagenome'],
                        help='Use QC thresholds for the specified sample type (Default: %(default)s)')
    parser.add_argument('--lib_enzyme', default=['unspecified'], nargs='+', type=str,
                        help='Name of the enzyme(s) used for Hi-C library preparation.')

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

    QC = HiCQC(outfile_prefix=args.outfile_prefix,
               sample_type=args.sample_type,
               thresholds_file=args.thresholds,
               rp_stats=args.rp_stats,
               mq_stats=args.mq_stats,
               edist_stats=args.edist_stats,
               lib_enzyme=args.lib_enzyme
               )

    QC.parse_bam(args.bam_file, max_read_pairs=args.num_reads)
    QC.plot_dup_saturation()
    QC.pass_judgement()
    QC.html_from_judgement()
    QC.plot_histograms()
    QC.stringify_stats()
    QC.log_stats()
    QC.write_stat_table()
    QC.write_dists_file()
    QC.write_pdf_report()
