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
        self.stat_dict, total_reads, num_dupes = b2mh.parse_bam_file(
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

    def test_greater_10k_contigs_right(self):
        self.assertEqual(self.stat_dict["GREATER_10K_CONTIGS"], 229)

    def test_count_zero_dist_pairs(self):
        self.assertEqual(self.stat_dict["ZERO_DIST_PAIRS"], 38)

    def test_count_num_pairs(self):
        self.assertEqual(self.stat_dict["NUM_PAIRS"], 107)

    def test_count_gt_10kbp(self):
        self.assertEqual(self.stat_dict["NUM_10KB_PAIRS"], 1)

    def test_count_gt_10kbp_actual(self):
        self.assertEqual(self.stat_dict["LARGE_INSERT_ACTUAL"], 1)

    def test_count_gt_10kbp_possible(self):
        self.assertEqual(self.stat_dict["LARGE_INSERT_POSSIBLE"], 65)

    def test_dists_right_len(self):
        self.assertEqual(sum(self.stat_dict["dists"].values()) + self.stat_dict["NUM_DIFF_CONTIG_PAIRS"], self.stat_dict["NUM_PAIRS"])

    def test_dists_right_num_zeros(self):
        num_zeros = self.stat_dict["dists"][0]
        self.assertEqual(num_zeros, self.stat_dict["ZERO_DIST_PAIRS"])

    def test_is_split_read_false(self):
        self.assertFalse(b2mh.is_split_read(self.example_read))

    def test_is_split_read_true(self):
        self.example_read.set_tag("SA", 1)
        self.assertTrue(b2mh.is_split_read(self.example_read))
    
    def test_represents_number(self):
        self.assertTrue(b2mh.represents_number(5))
        self.assertTrue(b2mh.represents_number(0))
        self.assertTrue(b2mh.represents_number(-5))
        self.assertTrue(b2mh.represents_number('5'))
        self.assertTrue(b2mh.represents_number('0'))
        self.assertTrue(b2mh.represents_number('-5'))
        self.assertTrue(b2mh.represents_number('+5'))
        
        self.assertTrue(b2mh.represents_number(5.0))
        self.assertTrue(b2mh.represents_number(0.0))
        self.assertTrue(b2mh.represents_number(-5.0))
        self.assertTrue(b2mh.represents_number('5.0'))
        self.assertTrue(b2mh.represents_number('0.0'))
        self.assertTrue(b2mh.represents_number('-5.0'))
        self.assertTrue(b2mh.represents_number('+5.0'))
        
        self.assertTrue(b2mh.represents_number(5.2))
        self.assertTrue(b2mh.represents_number(0.2))
        self.assertTrue(b2mh.represents_number(-5.2))
        self.assertTrue(b2mh.represents_number('5.2'))
        self.assertTrue(b2mh.represents_number('0.2'))
        self.assertTrue(b2mh.represents_number('-5.2'))
        self.assertTrue(b2mh.represents_number('+5.2'))
        
        self.assertFalse(b2mh.represents_number('a'))
        self.assertFalse(b2mh.represents_number('1a'))
        self.assertFalse(b2mh.represents_number('a1'))
        self.assertFalse(b2mh.represents_number('1 + 1'))
        self.assertFalse(b2mh.represents_number('1+1'))
        self.assertFalse(b2mh.represents_number(''))
        self.assertFalse(b2mh.represents_number(' '))
        self.assertFalse(b2mh.represents_number('\t'))
        self.assertFalse(b2mh.represents_number('hello world'))
        self.assertFalse(b2mh.represents_number('.'))
        self.assertFalse(b2mh.represents_number('+'))
        self.assertFalse(b2mh.represents_number('-'))
        
        self.assertTrue(b2mh.represents_number('inf'))
        self.assertFalse(b2mh.represents_number(None))
        self.assertFalse(b2mh.represents_number(dict()))
        self.assertFalse(b2mh.represents_number(list()))
        self.assertFalse(b2mh.represents_number([5]))
        self.assertFalse(b2mh.represents_number([5, 6]))
        self.assertFalse(b2mh.represents_number(['5']))
        self.assertFalse(b2mh.represents_number({'a':5}))
        self.assertFalse(b2mh.represents_number({5:5}))
        self.assertFalse(b2mh.represents_number({5:'a'}))
    
    def test_represents_int(self):
        self.assertTrue(b2mh.represents_int(5))
        self.assertTrue(b2mh.represents_int(0))
        self.assertTrue(b2mh.represents_int(-5))
        self.assertTrue(b2mh.represents_int('5'))
        self.assertTrue(b2mh.represents_int('0'))
        self.assertTrue(b2mh.represents_int('-5'))
        self.assertTrue(b2mh.represents_int('+5'))
        
        self.assertFalse(b2mh.represents_int(5.0))
        self.assertFalse(b2mh.represents_int(0.0))
        self.assertFalse(b2mh.represents_int(-5.0))
        self.assertFalse(b2mh.represents_int('5.0'))
        self.assertFalse(b2mh.represents_int('0.0'))
        self.assertFalse(b2mh.represents_int('-5.0'))
        self.assertFalse(b2mh.represents_int('+5.0'))
        
        self.assertFalse(b2mh.represents_int(5.2))
        self.assertFalse(b2mh.represents_int(0.2))
        self.assertFalse(b2mh.represents_int(-5.2))
        self.assertFalse(b2mh.represents_int('5.2'))
        self.assertFalse(b2mh.represents_int('0.2'))
        self.assertFalse(b2mh.represents_int('-5.2'))
        self.assertFalse(b2mh.represents_int('+5.2'))
        
        self.assertFalse(b2mh.represents_int('a'))
        self.assertFalse(b2mh.represents_int('1a'))
        self.assertFalse(b2mh.represents_int('a1'))
        self.assertFalse(b2mh.represents_int('1 + 1'))
        self.assertFalse(b2mh.represents_int('1+1'))
        self.assertFalse(b2mh.represents_int(''))
        self.assertFalse(b2mh.represents_int(' '))
        self.assertFalse(b2mh.represents_int('\t'))
        self.assertFalse(b2mh.represents_int('hello world'))
        self.assertFalse(b2mh.represents_int('.'))
        self.assertFalse(b2mh.represents_int('+'))
        self.assertFalse(b2mh.represents_int('-'))
        
        self.assertFalse(b2mh.represents_int('inf'))
        self.assertFalse(b2mh.represents_int(None))
        self.assertFalse(b2mh.represents_int(dict()))
        self.assertFalse(b2mh.represents_int(list()))
        self.assertFalse(b2mh.represents_int([5]))
        self.assertFalse(b2mh.represents_int([5, 6]))
        self.assertFalse(b2mh.represents_int(['5']))
        self.assertFalse(b2mh.represents_int({'a':5}))
        self.assertFalse(b2mh.represents_int({5:5}))
        self.assertFalse(b2mh.represents_int({5:'a'}))
    
    def test_make_pretty_stat_dict(self):
        test_dict = {
                'a' : '1234567',
                'b' : '0.123456',
                'c' : '1',
                'd' : '0',
                'e' : '1.0',
                'f' : '0.0',
                'g' : 1234567,
                'h' : 0.123456,
                'i' : 1,
                'j' : 0,
                'k' : 1.0,
                'l' : 0.0,
                'm' : 'abcd',
                'n' : None,
                'o' : '',
                'p' : dict(),
                'q' : list(),
                'r' : '+12345',
                's' : '+1.23456',
                't' : -1,
                'u' : -1.0,
                'v' : '-1.0',
                'w' : '-0.12345',
                'x' : '-0.0',
                'y' : '-0',
                'z' : -0
            }
        answer_dict = {
                'a' : '1,234,567',
                'b' : '12.35%',
                'c' : '1',
                'd' : '0',
                'e' : '100.0%',
                'f' : '0.0%',
                'g' : '1,234,567',
                'h' : '12.35%',
                'i' : '1',
                'j' : '0',
                'k' : '100.0%',
                'l' : '0.0%',
                'm' : 'abcd',
                'n' : None,
                'o' : '',
                'p' : dict(),
                'q' : list(),
                'r' : '12,345',
                's' : '123.46%',
                't' : '-1',
                'u' : '-100.0%',
                'v' : '-100.0%',
                'w' : '-12.35%',
                'x' : '0.0%',
                'y' : '0',
                'z' : '0'
            }
        pretty_dict = b2mh.make_pretty_stat_dict(test_dict)
        self.assertEqual(sorted(test_dict.keys()), sorted(pretty_dict.keys()))
        self.assertEqual(sorted(pretty_dict.keys()), sorted(answer_dict.keys()))
        for k in pretty_dict:
            self.assertEqual(pretty_dict[k], answer_dict[k],
                             'Key: {0}, Got {1}, Answer {2}'.format(k, pretty_dict[k], answer_dict[k]))

if __name__ == '__main__':
    unittest.main()
