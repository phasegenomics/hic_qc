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

import unittest
import bam_to_mate_hist as b2mh
import pysam

class MyTestCase(unittest.TestCase):
    def setUp(self):
        num_reads = 1000
        bamfile = "collateral/abc_test.bam"
        count_diff_refname_stub = False
        self.stat_dict = b2mh.parse_bam_file(
            num_reads=num_reads, bamfile=bamfile, count_diff_refname_stub=count_diff_refname_stub)

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
        self.assertEqual(self.stat_dict["NUM_DIFF_CONTIG_PAIRS"], 5)

    def test_count_splits(self):
        self.assertEqual(self.stat_dict["NUM_SPLIT_READS"], 6)

    def test_count_dupe_reads(self):
        self.assertEqual(self.stat_dict["NUM_DUPE_READS"], 2)

    def test_refs_right(self):
        self.assertEqual(len(self.stat_dict["refs"]), 1288)

    def test_count_zero_dist_pairs(self):
        self.assertEqual(self.stat_dict["ZERO_DIST_PAIRS"], 38)

    def test_count_num_pairs(self):
        self.assertEqual(self.stat_dict["NUM_PAIRS"], 107)

    def test_dists_right_len(self):
        self.assertEqual(len(self.stat_dict["dists"]), self.stat_dict["NUM_PAIRS"])

    def test_dists_right_num_zeros(self):
        num_zeros = list(self.stat_dict["dists"]).count(0)
        self.assertEqual(num_zeros, self.stat_dict["ZERO_DIST_PAIRS"])

    def test_is_split_read_false(self):
        self.assertFalse(b2mh.is_split_read(self.example_read))

    def test_is_split_read_true(self):
        self.example_read.set_tag("SA", 1)
        self.assertTrue(b2mh.is_split_read(self.example_read))


if __name__ == '__main__':
    unittest.main()
