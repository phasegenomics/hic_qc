#!usr/bin/python
'''
Max Press
August 22, 2018
Phase Genomics

hic_qc/test_hic_qc.py

This file contains unit tests for functions of the hic_qc.py script.

Copyright 2020, Phase Genomics Inc. All rights reserved.

The contents of this file are proprietary and private and are not intended for
distribution or use by any person or entity except Phase Genomics. You may not
use, modify, or distribute it in any fashion. You may not copy this file. You
may not describe the contents of this file to any other party.
'''

from __future__ import print_function
from __future__ import division

import os
import sys
import unittest
import hic_qc
import pysam

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.dirname = os.path.dirname(__file__)
        self.collateral_dir = self.dirname + "/collateral/"
        self.input_dir = self.collateral_dir + "input/"
        self.output_dir = self.collateral_dir + "output/"
        self.output_prefix = self.output_dir + "Read_mate_dist"
        num_reads = 1000
        bamfile = self.input_dir + "abc_test.bam"
        count_diff_refname_stub = False

        QC = hic_qc.HiCQC()
        QC.logger.setLevel("ERROR")
        QC.output_prefix = self.output_prefix
        QC.parse_bam(bamfile, max_read_pairs=num_reads)
        self.QC = QC
        self.stats = QC.stats
        # self.stat_dict, total_reads, num_dupes = hic_qc.parse_bam_file(
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

    def setUp(self):
        self.QCtmp = hic_qc.HiCQC()
        self.QCtmp.allowed_dupe_percentage = 0.5

        # Default stats
        self.QCtmp.stats['total_reads'] = 4000
        self.QCtmp.stats['total_read_pairs'] = 2000
        self.QCtmp.stats['total_read_pairs_hq'] = 1000
        self.QCtmp.stats['pairs_intracontig_hq'] = 200
        self.QCtmp.stats['pairs_on_contigs_greater_10k_hq'] = 600
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 4 # 2%
        self.QCtmp.stats['proximo_usable_rp'] = 101 # 5.05%
        self.QCtmp.stats['informative_pairs'] = 101 # >5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1000 # 50%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 4 # 0.67%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 5 # 0.5%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 10
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 800 # 20%
        self.QCtmp.stats['zero_dist_pairs'] = 400 # 20%
        self.QCtmp.stats['unmapped_reads'] = 400 # 10%
        self.QCtmp.stats['mapq0_reads'] = 800 # 20%

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

    def test_reads_spanning_up_to_1k(self):
        self.assertEqual(self.stats['reads_spanning_up_to_1k'], 101)

    def test_reads_spanning_1k_to_10k(self):
        self.assertEqual(self.stats['reads_spanning_1k_to_10k'], 1)

    def test_reads_spanning_10k_to_100k(self):
        self.assertEqual(self.stats['reads_spanning_10k_to_100k'], 1)

    def test_reads_spanning_100k_to_1000k(self):
        self.assertEqual(self.stats['reads_spanning_100k_to_1000k'], 2)

    def test_reads_spanning_greater_than_1000k(self):
        self.assertEqual(self.stats['reads_spanning_greater_than_1000k'], 1)

    def test_all_reads(self):
        self.assertEqual(self.stats['reads_spanning_greater_than_1000k'] + self.stats['reads_spanning_100k_to_1000k']
                         + self.stats['reads_spanning_10k_to_100k'] + self.stats['reads_spanning_1k_to_10k'] +
                         self.stats['reads_spanning_up_to_1k'] + self.stats['intercontig_pairs'], self.stats[
            'total_reads'] / 2)

    def test_count_num_pairs(self):
        self.assertEqual(self.stats['total_read_pairs'], 111)

    def test_count_gt_10kbp(self):
        self.assertEqual(self.stats['pairs_greater_10k'], 4)

    def test_count_gt_10kbp_actual(self):
        self.assertEqual(self.stats['pairs_greater_10k_on_contigs_greater_10k'], 4)

    def test_count_gt_10kbp_possible(self):
        self.assertEqual(self.stats['pairs_on_contigs_greater_10k'], 69)

    def test_dists_right_len(self):
        self.assertEqual(sum(self.QC.dists.values()) + self.stats['intercontig_pairs'], self.stats['total_read_pairs'])

    def test_dists_right_num_zeros(self):
        num_zeros = self.QC.dists[0]
        self.assertEqual(num_zeros, self.stats['zero_dist_pairs'])

    def test_is_split_read_false(self):
        self.QCtmp.update_read_stats(self.example_read)
        self.assertEqual(self.QCtmp.stats['split_reads'], 0)

    def test_is_split_read_true(self):
        self.example_read.set_tag("SA", 1)

        self.QCtmp.update_read_stats(self.example_read)
        self.assertEqual(self.QCtmp.stats['split_reads'], 1)

    def test_python_version(self):
        '''Confirms that PYTHON version in Travis CI env matches expectation'''
        if "PYTHON" in os.environ:
            expected_python = os.environ["PYTHON"] if os.environ["PYTHON"] != "default" else "3.6"
            version_string = "{}.{}".format(*sys.version_info)
            self.assertEqual(version_string, expected_python)
        else:
            return True

    def test_pass_judgement_sufficient(self):
        '''should pass'''
        # major denominators
        # driving metrics
        self.QCtmp.stats['noninformative_read_pairs'] = 100 # 10%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 400 # 66.7%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 500 # 50%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 1000
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 40 # 1%
        self.QCtmp.stats['zero_dist_pairs'] = 40 # 2%
        self.QCtmp.stats['unmapped_reads'] = 20 # 0.5%
        self.QCtmp.stats['mapq0_reads'] = 0 #0%
        self.QCtmp.pass_judgement()
        self.assertTrue(self.QCtmp.good_same_strand)
        self.assertTrue(self.QCtmp.good_informative_read_pairs)
        self.assertFalse(self.QCtmp.bad_noninformative_read_pairs)
        self.assertTrue(self.QCtmp.good_long_contacts)
        self.assertTrue(self.QCtmp.good_intercontig_contacts)
        self.assertTrue(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertTrue(self.QCtmp.judge_good)
        self.assertFalse(self.QCtmp.judge_bad)

    def test_pass_judgement_insufficient(self):
        '''should fail'''
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 1 # 0.5%
        self.QCtmp.stats['proximo_usable_rp'] = 10 # 0.25%
        self.QCtmp.stats['informative_pairs'] = 100 # >5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1500 # 75%
        # other good metrics
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 2000 # 50%
        self.QCtmp.stats['zero_dist_pairs'] = 1000 # 50%
        self.QCtmp.stats['unmapped_reads'] = 800 # 20%
        self.QCtmp.stats['mapq0_reads'] = 1200 # 30%
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertFalse(self.QCtmp.good_informative_read_pairs)
        self.assertTrue(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertTrue(self.QCtmp.high_dupe)
        self.assertTrue(self.QCtmp.many_zero_dist_pairs)
        self.assertTrue(self.QCtmp.many_unmapped_reads)
        self.assertTrue(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertTrue(self.QCtmp.judge_bad)

    def test_pass_judgement_mixed(self):
        # should be mixed results
        # major denominators
         # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 100 # 50%
        self.QCtmp.stats['proximo_usable_rp'] = 500 # 25%
        self.QCtmp.stats['informative_pairs'] = 100 # 5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1500 # 75%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 400 # 66.7%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 500 # 50%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 1000
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 2000 # 50%
        self.QCtmp.stats['zero_dist_pairs'] = 1000 # 50%
        self.QCtmp.stats['unmapped_reads'] = 800 # 20%
        self.QCtmp.stats['mapq0_reads'] = 1200 # 30%
        self.QCtmp.pass_judgement()
        self.assertTrue(self.QCtmp.good_same_strand)
        self.assertFalse(self.QCtmp.good_informative_read_pairs)
        self.assertTrue(self.QCtmp.bad_noninformative_read_pairs)
        self.assertTrue(self.QCtmp.good_long_contacts)
        self.assertTrue(self.QCtmp.good_intercontig_contacts)
        self.assertTrue(self.QCtmp.good_usable_reads)
        self.assertTrue(self.QCtmp.high_dupe)
        self.assertTrue(self.QCtmp.many_zero_dist_pairs)
        self.assertTrue(self.QCtmp.many_unmapped_reads)
        self.assertTrue(self.QCtmp.many_mapq_zero_reads)
        self.assertTrue(self.QCtmp.judge_good)
        self.assertTrue(self.QCtmp.judge_bad)

    def test_pass_judgement_low_signal(self):
        # should low signal
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 1 # 0.5%
        self.QCtmp.stats['proximo_usable_rp'] = 10 # 0.25%
        self.QCtmp.stats['informative_pairs'] = 101 # >5%
        self.QCtmp.stats['noninformative_read_pairs'] = 100 # 10%
        # other good metrics
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 40 # 1%
        self.QCtmp.stats['zero_dist_pairs'] = 40 # 2%
        self.QCtmp.stats['unmapped_reads'] = 20 # 0.5%
        self.QCtmp.stats['mapq0_reads'] = 0 #0%
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertTrue(self.QCtmp.good_informative_read_pairs)
        self.assertFalse(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertFalse(self.QCtmp.judge_bad)

    def test_pass_judgement_close_not_sufficient(self):
        # should barely not pass
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 3 # 1.5% - exactly on the threshold fails
        self.QCtmp.stats['proximo_usable_rp'] = 500 # 25%
        self.QCtmp.stats['informative_pairs'] = 101 # >5%
        self.QCtmp.stats['noninformative_read_pairs'] = 100 # 10%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 400 # 66.7%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 500 # 50%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 1000
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 40 # 1%
        self.QCtmp.stats['zero_dist_pairs'] = 40 # 2%
        self.QCtmp.stats['unmapped_reads'] = 20 # 0.5%
        self.QCtmp.stats['mapq0_reads'] = 0 #0%
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertTrue(self.QCtmp.good_informative_read_pairs)
        self.assertFalse(self.QCtmp.bad_noninformative_read_pairs)
        self.assertTrue(self.QCtmp.good_long_contacts)
        self.assertTrue(self.QCtmp.good_intercontig_contacts)
        self.assertTrue(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertFalse(self.QCtmp.judge_bad)

    def test_pass_judgement_close_still_insufficient(self):
        # should still be insufficient just over the line
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 1 # 0.5%
        self.QCtmp.stats['proximo_usable_rp'] = 10 # 0.25%
        self.QCtmp.stats['informative_pairs'] = 100 # >5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1001 # 50.05%
        # other good metrics
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 2000 # 50%
        self.QCtmp.stats['zero_dist_pairs'] = 1000 # 50%
        self.QCtmp.stats['unmapped_reads'] = 800 # 20%
        self.QCtmp.stats['mapq0_reads'] = 1200 # 30%
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertFalse(self.QCtmp.good_informative_read_pairs)
        self.assertTrue(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertTrue(self.QCtmp.high_dupe)
        self.assertTrue(self.QCtmp.many_zero_dist_pairs)
        self.assertTrue(self.QCtmp.many_unmapped_reads)
        self.assertTrue(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertTrue(self.QCtmp.judge_bad)

    def test_pass_judgement_on_boundaries_should_be_low_signal(self):
        # should be low signal
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 3 # 1.5%
        self.QCtmp.stats['proximo_usable_rp'] = 100 # 5%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 15 # 2.5%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 25 # 2.5%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 600
        # noninformative breakdown
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertTrue(self.QCtmp.good_informative_read_pairs)
        self.assertFalse(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertFalse(self.QCtmp.judge_bad)

    def test_pass_judgement_just_past_boundaries_should_be_mixed_results(self):
        # should be low signal
        # major denominators
        # driving metrics
        self.QCtmp.stats['informative_pairs'] = 100 # 5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1001 # 50.05%
        # other good metrics
        self.QCtmp.stats['pairs_greater_10k_on_contigs_greater_10k_hq'] = 19 # 3.17%
        self.QCtmp.stats['pairs_intercontig_hq_gt10kbp'] = 26 # 2.6%
        self.QCtmp.stats['proximo_usable_rp_hq_per_ctg_gt_5k'] = 601
        # noninformative breakdown
        self.QCtmp.stats['duplicate_reads'] = 801 # 20.025%
        self.QCtmp.stats['zero_dist_pairs'] = 401 # 20.05%
        self.QCtmp.stats['unmapped_reads'] = 401 # 10.025%
        self.QCtmp.stats['mapq0_reads'] = 801 # 20.025%
        self.QCtmp.pass_judgement()
        self.assertTrue(self.QCtmp.good_same_strand)
        self.assertFalse(self.QCtmp.good_informative_read_pairs)
        self.assertTrue(self.QCtmp.bad_noninformative_read_pairs)
        self.assertTrue(self.QCtmp.good_long_contacts)
        self.assertTrue(self.QCtmp.good_intercontig_contacts)
        self.assertTrue(self.QCtmp.good_usable_reads)
        self.assertTrue(self.QCtmp.high_dupe)
        self.assertTrue(self.QCtmp.many_zero_dist_pairs)
        self.assertTrue(self.QCtmp.many_unmapped_reads)
        self.assertTrue(self.QCtmp.many_mapq_zero_reads)
        self.assertTrue(self.QCtmp.judge_good)
        self.assertTrue(self.QCtmp.judge_bad)

    def test_pass_judgement_just_across_sufficient_boundary(self):
        # should pass
        # major denominators
        # driving metrics
        # other good metrics
        # noninformative breakdown
        self.QCtmp.pass_judgement()
        self.assertTrue(self.QCtmp.good_same_strand)
        self.assertTrue(self.QCtmp.good_informative_read_pairs)
        self.assertFalse(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertTrue(self.QCtmp.judge_good)
        self.assertFalse(self.QCtmp.judge_bad)

    def test_pass_judgement_just_across_insufficient_boundary(self):
        # should pass
        # major denominators
        # driving metrics
        self.QCtmp.stats['pairs_on_same_strand_hq'] = 3 # 1.5%
        self.QCtmp.stats['proximo_usable_rp'] = 100 # 5%
        self.QCtmp.stats['informative_pairs'] = 100 # 5%
        self.QCtmp.stats['noninformative_read_pairs'] = 1001 # 50.05%
        # other good metrics
        # noninformative breakdown
        self.QCtmp.pass_judgement()
        self.assertFalse(self.QCtmp.good_same_strand)
        self.assertFalse(self.QCtmp.good_informative_read_pairs)
        self.assertTrue(self.QCtmp.bad_noninformative_read_pairs)
        self.assertFalse(self.QCtmp.good_long_contacts)
        self.assertFalse(self.QCtmp.good_intercontig_contacts)
        self.assertFalse(self.QCtmp.good_usable_reads)
        self.assertFalse(self.QCtmp.high_dupe)
        self.assertFalse(self.QCtmp.many_zero_dist_pairs)
        self.assertFalse(self.QCtmp.many_unmapped_reads)
        self.assertFalse(self.QCtmp.many_mapq_zero_reads)
        self.assertFalse(self.QCtmp.judge_good)
        self.assertTrue(self.QCtmp.judge_bad)

    def test_empty_bam(self):
        bamfile = self.input_dir + "abc_test.empty.bam"
        self.QCtmp.parse_bam(bamfile, max_read_pairs=1000)
        self.QCtmp.plot_dup_saturation()
        self.QCtmp.pass_judgement()
        self.QCtmp.html_from_judgement()
        self.QCtmp.plot_histograms()
        self.QCtmp.stringify_stats()
        self.QCtmp.log_stats()
        self.QCtmp.write_stat_table()
        self.QCtmp.write_dists_file()
        self.QCtmp.write_pdf_report(quiet=True)

    def test_bad_report_path(self):
        self.QCtmp.paths["script_dir"] = "/not/a/path"
        with self.assertRaises(FileNotFoundError):
            self.QCtmp.write_pdf_report(quiet=True)

if __name__ == '__main__':
    unittest.main()
