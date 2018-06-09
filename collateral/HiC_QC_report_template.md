# Hi-C library QC report

&copy; 2018 Phase Genomics Inc.

## Library statistics
<center>

| Label                                                           | Your library          | Expected values                               |
|-----------------------------------------------------------------|:---------------------:|----------------------------------------------:|
| BAM file                                                        | {BAM_FILE_PATH}                      | N/A                                           |
| Number of read pairs                                            | {NUM_PAIRS}             | N/A                                           |
| Fraction of read pairs >10KB apart                              | {NUM_10KB_PAIRS}        | 0.005-0.05                                    |
| Fraction of read pairs mapping to different contigs/chromosomes | {NUM_DIFF_CONTIG_PAIRS} | 0.1-0.5 (contigs)<br>0.01-0.1 (chromosomes)      |
| Fraction of split reads                                         | {NUM_SPLIT_READS}       | 0.1-0.4 (PG libraries) 0.3+ (other libraries) |
| Fraction of zero-distance pairs                                 | {ZERO_DIST_PAIRS}       | 0-0.15                                        |


</center>
See below for information on differences between Phase Genomics Hi-C libraries and traditional Hi-C libraries.

## Aligned mate distance histograms
!["Long range interaction histogram"]({PATH_TO_LONG_HIST})
!["Short range interaction histogram"]({PATH_TO_SHORT_HIST})


## Alignment distance statistics and plots
We briefly describe some of the statistics we compute below to aid interpretation of this report.
### Fraction of read pairs > 10KB apart
### Fraction of read pairs mapping to different contigs or chromosomes
### Fraction of zero-distance pairs

## Split reads
Traditionally, split reads have been a favored measure of Hi-C library quality. This is because traditional Hi-C library preparations are expected to produce many reads reading through junctions. 

Phase Genomics libraries, whether produced in our laboratory or by means of ProxiMeta &copy;, Animal, Plant, or Human kits, will have a generally lower fraction of split reads. This is because we have optimized our Hi-C protocol to 
