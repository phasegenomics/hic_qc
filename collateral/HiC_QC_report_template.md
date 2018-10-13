# Hi-C library QC report

&copy; 2018 Phase Genomics Inc.

## Library statistics
<center>

| Label                                                           | Your library          | Expected values                               |
|-----------------------------------------------------------------|:---------------------:|----------------------------------------------:|
| BAM file                                                        | {BAM_FILE_PATH}         | N/A                                           |
| Assembly size                                                   | {TOTAL_LEN}             | N/A                                           |
| Contig N50                                                      | {N50}                   | N/A                                           |
| Contigs                                               | {NUM_CONTIGS}           | N/A                                           |
| Contigs greater than 10kbp                             | {GREATER_10K_CONTIGS}   | N/A                                           |
| Total read pairs analyzed                                   | {NUM_PAIRS}             | N/A                                           |
| Read pairs >10kbp apart                              | {NUM_10KB_PAIRS}        | 1-15%                                    |
| Read pairs >10kbp apart on contigs >10kbp    | {LARGE_INSERT_PROPORTION} | 1-15%                                   |
| Read pairs mapping to different contigs/chromosomes | {NUM_DIFF_CONTIG_PAIRS} | 10-60% (contigs)<br>1-20% (chromosomes)      |
| Split reads                                         | {NUM_SPLIT_READS}       | 1-10% (PG libraries)<br> 30%+ (other libraries) |
| Zero-distance pairs                                 | {ZERO_DIST_PAIRS}       | 0-20%                                        |
| Zero map quality reads                         | {MAPQ0_READS}           | 0-10%                                        |
| Duplicate reads*                                    | {NUM_DUPE_READS}        | 0-10%                                        |
| Duplicate reads (extrapolated to {TARGET_READ_TOTAL} reads)* | {NUM_DUPE_READS_EXTRAP} | 0-50%                               |
| Subjective Hi-C library judgement                           | {JUDGEMENT}             | see below           |
</center>

*If this quantity is zero, see duplicate read section below. If negative, there are too few reads sampled to estimate duplicates.
See below for information on differences between Phase Genomics Hi-C libraries and traditional Hi-C libraries.

## Aligned mate distance histograms

!["Long range interaction histogram"]({PATH_TO_LONG_HIST})
!["Short range interaction histogram"]({PATH_TO_SHORT_HIST})

!["Log-log interaction histogram"]({PATH_TO_LOG_LOG_HIST})

## Duplicate read saturation curve

!["Duplicate read saturation curve"]({PATH_TO_DUP_SAT})

## Alignment distance statistics and plots
We briefly describe some of the statistics we compute below to aid interpretation of this report.

### Subjective Hi-C Library Judgement
While Hi-C data is nuanced and some analyses are more sensitive to data quality than others,  a basic quality assessment can usually be made by examining the mapping characteristics of the Hi-C library. Based on our experience working with Hi-C data, we classify libraries into one of four QC categories:

- **Pass** means that from everything we can tell, the library looks to be in great shape. Proceed to full sequencing or analysis with confidence.

- **Mixed Results** means that the library is good in some ways, but not in others. Perhaps it has a good amount of long range data, but there are also an elevated number of read pairs with MAPQ 0. Usually, data generated from Mixed Results libraries works out just fine (a high MAPQ 0 number can be due to repetitiveness in the assembly, for example), but it is good to know there may have been a few hiccups in the library prep in case troubleshooting is needed down the line.

- **Low Signal** means that the library contains good Hi-C signal, but it's in lower proportion than usual. These libraries are generally good for generating useful Hi-C data, but you may need to sequence a little deeper than normal to get enough of it. You might consider size selecting the library to discard reads outside the 300-700bp range, as these are unlikely to be good Hi-C junctions. Alternatively, you might just want to prep a new library.

- **Fail** means that the library, or perhaps the library in combination with a low-contiguity or error-prone assembly, does not look useful. Sometimes size selection can rescue such libraries, but sometimes a new prep is the only way forward. (Contact us)[mailto:support@phasegenomics.com] if you get a fail and we will help you out.

### Read pairs > 10kbp apart
These are the proportion of read pairs which map to the same contig, with at least 10kbp separating them. More is always better, but because this number is affected by assembly contiguity, there is not a specifc target threshold. Note that for some analyses, such as scaffolding or metagenomic deconvolution, read pairs that map to the same contig are not useful because they do not provide information that the assembly doesn't already contain. This statistic is more useful for these projects because it correlates with library prep success. These reads are useful for analyses like structural variant analysis or assembly misjoin detection, because they provide detailed structural information about existing assembled sequences.

### Read pairs > 10kbp apart mapping to contigs >10kbp
These are the proportion of read pairs which map to the same contig, with at least 10kbp separating them, but only considering read pairs mapping to contigs that are at least 10kbp long. This attempts to corrects for assembly contiguity differences. More is always better, but typically at least 5% is desired. Note that for some analyses, such as scaffolding or metagenomic deconvolution, read pairs that map to the same contig are not useful because they do not provide information that the assembly doesn't already contain. This statistic is more useful for these projects because it correlates with library prep success. These reads are useful for analyses like structural variant analysis or assembly misjoin detection, because they provide detailed structural information about existing assembled sequences.


### Read pairs mapping to different contigs or chromosomes
These are the proportion of read pairs which map to different contigs, which is particularly important . More is always better, but because this number is affected by assembly quality, there is not a specifc target threshold, although at least 20% on *de novo* assembly projects is helpful. These reads are the primary source of information for Hi-C scaffolding or metagenomic deconvolution analyses. This statistic useful on most *de novo* proects because it also correlates with library prep success. These reads may be useful for analyses like structural variant analysis or assembly misjoin detection if those are performed on lower contiguity assemblies, because they provide detailed structural information about sequences which were not assembled together into contigs.


### Zero-distance pairs
A failure mode for Hi-C libraries is a large number of read pairs which map as if there were no distance between them. This can occur for many reasons, but common ones include very short library inserts (possible overdigestion) or sometimes mapping artifacts. These can also be an issue on some Illumina instruments such as iSeq and HiSeqX. These essentially represent read pairs which provide no long-range interaction information. Less is always better for these kinds of read pairs, although even libraries with up to 70-80% zero-distance pairs can yield usable results.


### Split reads
Traditionally, split reads have been a favored measure of Hi-C library quality because they directly exhibit Hi-C junctions. Most traditional Hi-C library preparations produce many reads that sequence through junctions because their Hi-C junctions tend to occur randomly on the proximity ligated chimeric molecules. However, Phase Genomics has introduced steps into our protocol which tend to make junctions occur closer to the center of the proximity ligated chimeric molecules.  

Phase Genomics libraries, whether produced in our laboratory or by means of our Plant, Animal, Human, or Microbe Hi-C kits, will have a generally lower fraction of split reads. This is because we have optimized our Hi-C protocol to enrich for slightly longer fragments around Hi-C junctions, such that each read is less likely to read through a junction even when a junction is present.

This innovation improves mappability and increases the amount of useful data, and reduces the utility of split read measurements to assess library quality. We therefore rely more heavily on metrics that directly relate to the usefulness of Hi-C reads for proximity analysis, such as the percentage of read pairs with mates mapping far away, or mapping to different contigs.

### Duplicate reads
**IMPORTANT NOTE: THE DUPLICATE FLAG IS NOT SET BY DEFAULT IN A BAM FILE. YOU NEED TO EXPLICITLY SET IT BY E.G. RUNNING SAMBLASTER ON YOUR BAM FILE. IF THE FRACTION OF DUPLICATES IS EXACTLY ZERO, IT PROBABLY MEANS THAT THE FLAG HAS NOT BEEN SET.**

Sequencing libraries frequently contain duplicate reads due to PCR or optical issues. These are generally considered to be non-informative because they are chemical artifacts rather than biological signal, and are thus typically excluded from further analysis. Higher proportions of duplicate reads are also correlated with low library complexity and poor library performance, making the proportion of duplicate reads a useful quality control measure.

Also, because we often recommend sequencing a few million read pairs for QC prior to a full sequencing run, it is helpful to project the amount of duplicate reads a full run might generate. To do this, we attempt to fit a curve to the rate at which duplicates are observed in a library, and then project the proportion of duplicates that would be expected in a deeper sequencing run. This projection is a reasonably useful heuristic, but it is not a guarantee that the true duplicate rate will be near a specific value. QC sequencing data with a very low number of reads or which mapped well at a very low rate can distort this calculation. We attempt to identify low confidence calculations and report that we could not extrapolate the expected duplicate frequency when it is possible to do so. 
<br>
<h4 class="centered">
REPORT VERSION: COMMIT_VERSION
</h4>
