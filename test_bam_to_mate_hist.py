#!usr/bin/python
'''
Max Press
August 22, 2018
Phase Genomics

bam_to_mate_hist/test_bam_to_mate_hist.py

This file contains unit tests for functions of the bam_to_mate_hist.py script.

Copyright 2018, Phase Genomics Inc. All rights reserved.

The contents of this file are proprietary and private and are not intended for
distribution or use by any person or entity except Phase Genomics. You may not
use, modify, or distribute it in any fashion. You may not copy this file. You
may not describe the contents of this file to any other party.
'''

from __future__ import print_function
from __future__ import division

import unittest
import bam_to_mate_hist as b2mh
import pysam

class MyTestCase(unittest.TestCase):
    def setUp(self):
        num_reads = 1000
        bamfile = "collateral/abc_test.bam"
        count_diff_refname_stub = False

        QC = b2mh.HiCQC()
        QC.parse_bam(bamfile, max_read_pairs=num_reads)
        self.QC = QC
        self.stats = QC.stats
        # self.stat_dict, total_reads, num_dupes = b2mh.parse_bam_file(
        #     num_reads=num_reads, bamfile=bamfile, count_diff_refname_stub=count_diff_refname_stub)

        self.example_read = pysam.AlignedSegment()
        self.example_read.reference_start = 30
        self.example_read.query_name = 'read1'
        self.example_read.mapping_quality = 30
        self.example_read.query_sequence = "AAAAACAAAACAAAAT"
        self.example_read.query_qualities = [30] * 16
        self.example_read.cigarstring = '16M'
        self.example_read.set_tag("NM", 0)
        self.example_read.set_tag("MD", 100)
        self.example_read.set_tag("AS", 100)
        self.example_read.set_tag("XS", 0)

    def tearDown(self):
        pass

    # all manually measured in the BAM file...
    def test_count_diff_chr_pairs(self):
        self.assertEqual(self.stats['intercontig_pairs'], 5)

    def test_count_splits(self):
        self.assertEqual(self.stats['split_reads'], 6)

    def test_count_dupe_reads(self):
        self.assertEqual(self.stats['duplicate_reads'], 2)

    def test_refs_right(self):
        self.assertEqual(len(self.QC.refs), 1288)

    def test_greater_10k_contigs_right(self):
        self.assertEqual(len(self.QC.contigs_greater_10k), 229)

    def test_count_zero_dist_pairs(self):
        self.assertEqual(self.stats['zero_dist_pairs'], 38)

    def test_count_num_pairs(self):
        self.assertEqual(self.stats['total_read_pairs'], 107)

    def test_count_gt_10kbp(self):
        self.assertEqual(self.stats['pairs_greater_10k'], 1)

    def test_count_gt_10kbp_actual(self):
        self.assertEqual(self.stats['pairs_greater_10k_on_contigs_greater_10k'], 1)

    def test_count_gt_10kbp_possible(self):
        self.assertEqual(self.stats['pairs_on_contigs_greater_10k'], 65)

    def test_dists_right_len(self):
        self.assertEqual(sum(self.QC.dists.values()) + self.stats['intercontig_pairs'], self.stats['total_read_pairs'])

    def test_dists_right_num_zeros(self):
        num_zeros = self.QC.dists[0]
        self.assertEqual(num_zeros, self.stats['zero_dist_pairs'])

    def test_is_split_read_false(self):
        QCtmp = b2mh.HiCQC()
        QCtmp.update_read_stats(self.example_read)

        self.assertEqual(QCtmp.stats['split_reads'], 0)

    def test_is_split_read_true(self):
        QCtmp = b2mh.HiCQC()
        self.example_read.set_tag("SA", 1)

        QCtmp.update_read_stats(self.example_read)
        self.assertEqual(QCtmp.stats['split_reads'], 1)


if __name__ == '__main__':
    unittest.main()
