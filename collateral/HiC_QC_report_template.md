# Hi-C library QC report

&copy; 2018 Phase Genomics Inc.

## Library statistics
<center>

| Label                                                           | Your library          | Expected values                               |
|-----------------------------------------------------------------|:---------------------:|----------------------------------------------:|
| BAM file                                                        | {bamname}         | N/A                                           |
| Assembly size                                                   | {total_length}             | N/A                                           |
| Contig N50                                                      | {N50}                   | N/A                                           |
| Number of contigs                                               | {contigs}           | N/A                                           |
| Number of contigs greater than 10KB                             | {contigs_greater_10k}   | N/A                                           |
| Number of read pairs analyzed                                   | {total_read_pairs}             | N/A                                           |
| Percent of read pairs >10KB apart                              | {perc_pairs_greater_10k}        | 1-15%                                    |
| Percent of read pairs >10KB apart mapping to contigs >10Kbp    | {perc_pairs_greater_10k_on_contigs_greater_10k} | 1-15%                                   |
| Percent of read pairs mapping to different contigs/chromosomes | {perc_intercontig_pairs} | 1-6% (contigs)<br>1-20% (chromosomes)      |
| Percent of split reads                                         | {perc_split_reads}       | 1-10% (PG libraries) 30%+ (other libraries) |
| Percent of zero-distance pairs                                 | {perc_zero_dist_pairs}       | 0-20%                                        |
| Percent of reads with zero map quality                         | {perc_mapq0_reads}           | 0-10%                                        |
| Percent of duplicate reads*                                    | {perc_duplicate_reads}        | 0-10%                                        |
| Percent of duplicate reads (extrapolated to {target_read_total} reads)* | {extrapolated_dup_rate} | 0-50%                               |
| Percent of unmapped reads                                      | {perc_unmapped_reads} | 0-10%                               |
| (SUBJECTIVE!!) Hi-C library judgment                            | {judgment}             | pass/fail/mixed-results/low-signal           |
</center>

*If this quantity is zero, see duplicate read section below. If quantity is negative, there are too few reads sampled to estimate duplicates.

See below for information on differences between Phase Genomics Hi-C libraries and traditional Hi-C libraries.


## Aligned mate distance histograms

!["Long range interaction histogram"]({long_hist})
!["Short range interaction histogram"]({short_hist})
!["Log-log interaction histogram"]({log_log_hist})

## Duplicate read saturation curve

!["Duplicate read saturation curve"]({dup_sat_curve})

## Alignment distance statistics and plots
We briefly describe some of the statistics we compute below to aid interpretation of this report.

### Percent of read pairs > 10KB apart
More is better.

### Percent of read pairs > 10KB apart mapping to contigs >10Kbp
More is better. Attempts to correct for assembly contiguity.

### Percent of read pairs mapping to different contigs or chromosomes
More is better.

### Percent of zero-distance pairs
Less is better. Seems to be due to very short library inserts or sometimes mapping artifacts; can be an issue on some Illumina instruments such as iSeq and HiSeqX. Still, even libraries with up to 70% zero-distance pairs can yield usable results sometimes.

### Percent of split reads
Traditionally, split reads have been a favored measure of Hi-C library quality. This is because traditional Hi-C library preparations are expected to produce many reads reading through junctions.

Phase Genomics libraries, whether produced in our laboratory or by means of ProxiMeta &copy;, Animal, Plant, or Human kits, will have a generally lower Percent of split reads. This is because we have optimized our Hi-C protocol to enrich for slightly longer fragments around Hi-C junctions, such that each read is less likely to read through a junction even when a junction is present. This innovation improves mappability.

We therefore rely more heavily on metrics that directly relate to the usefulness of Hi-C reads for proximity analysis, such as the percentage of read pairs with mates mapping far away, or mapping to different contigs.

### Duplicate reads
**IMPORTANT NOTE: THE DUPLICATE FLAG IS NOT SET BY DEFAULT IN A BAM FILE. YOU NEED TO EXPLICITLY SET IT BY E.G. RUNNING SAMBLASTER ON YOUR BAM FILE. IF THE PERCENT OF DUPLICATES IS EXACTLY ZERO, IT PROBABLY MEANS THAT THE FLAG HAS NOT BEEN SET.**

Sequencing libraries frequently contain duplicate reads, due to PCR or optical issues. These are generally considered to be non-informative, and are thus bad. Higher proportions of duplicate reads are also correlated with low library complexity and poor library performance generally.

### Percent of unmapped reads
A high percent of unmapped reads may indicate sequence is missing from the reference, the reads are mapped to the wrong reference, or the sample is contaminated.

#### REPORT VERSION: COMMIT_VERSION
