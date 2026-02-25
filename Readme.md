[![Build Status](https://travis-ci.com/phasegenomics/hic_qc.svg?branch=master)](https://travis-ci.com/phasegenomics/hic_qc)

# `hic_qc.py` readme.

This script is intended as a simple QC method for Hi-C libraries, based on reads in a BAM file aligned to some genome/assembly. For our full recommendations on aligning and QCing Hi-C data, please [see here](https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html).

The most informative Hi-C reads are the ones that are long-distance contacts, or contacts between contigs of an assembly. This tool quantifies such contacts and makes plots of contact distance distributions. The most successful Hi-C libraries have many long-distance and multi-contig contacts. Note that for chromosome-scale assemblies, long-distance contacts may not be between seperate contigs.

Hi-C connectivity drops off in approximately a power-law with increasing linear sequence distance. Consequently, one expects Hi-C reads to follow a characteristic distribution, wherein there is a spike of many read pairs at distances close to zero, which drops off smoothly (in log space) with increasing distance. If there are odd spikes or discontinuities, or if there are few long-distance contacts, there may be a problem either with the library or the assembly. You may also see such abberations in cases where the data you are aligning comes from a different individual than the assembly and structural variants are expected.

## Dependencies
* python >=3.8
* numpy
* pysam
* matplotlib
* scipy
* markdown
* pdfkit
* wkhtmltopdf (must be installed manually if using PDF report generation)

## Installation

### Conda installation

These steps assume you already have **conda** and **git** installed.

```bash
git clone https://github.com/phasegenomics/hic_qc.git
cd hic_qc
```

# Create and activate the conda environment
conda env create -n hic_qc --file env.minimal.yml
conda activate hic_qc

# Install dependencies from pip that are not available/reliable in conda
pip install pdfkit

# Install hic_qc into the environment
pip install .

## PDF Report Generation (wkhtmltopdf Required)

`hic_qc` can generate a PDF QC report in addition to an html.  
PDF generation requires the external `wkhtmltopdf` binary to be installed on your system.

If `wkhtmltopdf` is not installed, you can still run `hic_qc` normally by using:

```bash
python3 hic_qc.py --disable_report -b sample.bam
```

## Installing wkhtmltopdf

### macOS (Homebrew)

```bash
brew install wkhtmltopdf
```

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y wkhtmltopdf
```

### RHEL / CentOS / Rocky Linux

```bash
sudo yum install -y wkhtmltopdf
```

### Verify wkhtmltopdf

```bash
wkhtmltopdf --version
```

## Usage
In the most basic use-case, you can run the script in a terminal

`python3 hic_qc.py -b input.bam -n num_reads_to_use`

where `input.bam` is your BAM file from aligning Hi-C reads to your reference, and `num_reads_to_use` is just the number of read pairs you want to sample from the BAM file (default 1 million read pairs; assuming there are this many reads in the file).

The script will write plots in PNG format (long, short, and log-log mate distance histograms) with the prefix `Read_mate_dist` in the working directory, unless `-o` or `--outfile_prefix` are set as described below.

The script will also quantify some basic QC metrics and print those to the screen.

The script will generate an HTML report of those metrics with the plots embedded, and will attempt to generate a PDF report as well (PDF generation requires `wkhtmltopdf`). To disable PDF generation, use the `--disable_report` flag.

Coverage metrics are computed for autosome-style references when available. For non-human/non-autosome references, coverage metrics are automatically skipped, or you can explicitly disable coverage with `-c` / `--disable_coverage`.

To set the name of the files written out, such as the PNG figures and the report files, set the `-o /path/to/outfile` or `--outfile_prefix /path/to/outfile` parameters.

QC is performed using a set of thresholds in JSON format. By default, the file `hic_qc/collateral/thresholds.json` is used. The chosen file may be changed with the `--thresholds` flag. Note that the thresholds in the default file are informed by Phase Genomics' analysis of thousands of Hi-C libraries, and reflect what we ourselves use for QC.

Different QC thresholds may be present in a thresholds file. The default file includes thresholds for genome scaffolding projects and metagenome deconvolution projects. The `--sample_type` argument is used to specify which set of thresholds in the thresholds file should be used in the run, and is also noted at the top of the report. By default, the `genome` sample type is used. Additional sample types may be added to a thresholds JSON file by making them keys in the file.

## Library judgement and thresholds

### Judgement categories
The report includes a judgement about the library and the assembly it was mapped to based on the observed statistics, shown at the top of the report. Libraries are given one of four classifications:

* **SUFFICIENT** - the library and assembly appear to be sufficient for the purposes shown in the report.
* **INSUFFICIENT** - the library and assembly appear to be insufficient for the purposes shown in the report.
* **MIXED RESULTS** - the library and assembly are probably sufficient for the purposes shown in the report, but there is some additional noise or other unexpected properties in the report.
* **LOW SIGNAL** - the library and assembly do not appear to be actively bad, but there is not very much observable long-range Hi-C signal. It is likely they are insufficient for the purposes shown in the report.

IMPORTANT NOTE: because the input assembly is a significant contributor to the ability to perform a given analysis, a good library can still generate a failed result when mapped to a poor assembly. This can particularly occur with difficult-to-align assemblies, such as polyploid, highly repetitive, or extremely fragmented assemblies. Using a related species instead of an exact organism match can also negatively impact mapping percentages and therefore the utility of the QC report.

### "Good" properties
Several statistics are used to determine whether a set of alignments show strong long-range Hi-C signal:

* **HQ RPs >10KB apart (CTGs >10KB)**: the percentage of read pairs that map with high quality (MAPQ ≥ 20, max edit distance ≤ 5, not duplicates) in which both mates align to the same contig, the contig is at least 10kbp long, and the mates are at least 10kbp apart.
* **Intercontig HQ RPs (CTGs >10KB)**: the percentage of read pairs that map with high quality in which each mate aligns to a different contig, and each of those contigs is at least 10kbp. Highly contiguous assemblies, including chromosome-scale assemblies, will naturally show a reduced percentage due to fewer opportunities for read pairs to map to distinct contigs.
* **Same strand HQ RPs**: the percentage of read pairs that map with high quality to the same strand. Such reads are strong indicators of true Hi-C junctions.

These metrics are compared to thresholds defined in the thresholds JSON file to determine whether the dataset demonstrates sufficient long-range Hi-C signal.

### "Bad" properties
Several statistics are used to determine if a set of alignments shows problematic properties:

* **Duplicate reads**: the percentage of reads flagged as PCR duplicates (by a prior tool such as Picard or, our recommendation, SAMBLASTER). Duplicate rates can increase with deeper sequencing, but unusually high values (e.g., >5% at low depth) can indicate library preparation issues.
* **Zero map quality reads**: the percentage of reads aligning with MAPQ = 0. These reads align non-uniquely and are common in repetitive or low-complexity assemblies.
* **Unmapped reads**: the percentage of reads that could not be mapped. High unmapped rates can result from a low-quality draft assembly, contamination, or a mismatch between the Hi-C data and the reference assembly.

These metrics are also compared against thresholds defined in the thresholds JSON file.

### Generating the final judgement
Thresholds for the "good" and "bad" aspects are defined in the specified thresholds JSON file and are displayed in the report. Fields in the report are highlighted to indicate whether each metric passed or failed its threshold. These comparisons are used to generate the final judgement:

* **SUFFICIENT** - strong Hi-C signal and no major problematic properties
* **INSUFFICIENT** - weak Hi-C signal and problematic properties
* **MIXED RESULTS** - strong Hi-C signal but also problematic properties
* **LOW SIGNAL** - weak Hi-C signal without major problematic properties

## Histogram plot characteristics to look for
Histogram plots should show some characteristic features:

* Substantial long-range contacts (note that contact distance is bounded by the assembly). You will want to see at least some contacts approximately as long as your longest contig. In the log-log histogram, the appearance of a second hump or positive slope after the initial dropoff is a very good qualitative sign.
* Gradual drop-off in signal with increasing distance (in log space). Choppiness or spikes in the distribution may indicate problems such as collapsed repeats or chimerisms, unless it can be attributed to sampling error due to (very) small numbers of reads. Periodicity in the distribution with distance often indicates problems.
* The leftmost spike of mates mapping very close is always the most prominent feature in the plot. However, it should not be too much larger than the rest of the distribution, or you are not having enough long-distance contacts. A dropoff of 3-4 orders of magnitude in the 0-20KB plot is the most you want to see in that sudden dropoff. Ideally it would be only 1-2 orders of magnitude dropoff in the 0-20KB range.

## Example histograms
The collateral folder includes several histograms which serve as examples of what is expected for a good Hi-C library.

## Statistics reported
* Number of read pairs with mates mapping to exactly the same position (distance == 0). These are bad. We observe these reads at some rate all the time, but they are especially abundant when there is a problem. This proportion should be small, no more than 10% and ideally much smaller. That said, if other measures look ok it might be worth trying a library even if there are many distance == 0 pairs.
* Number of read pairs with mates mapping >10KB apart. These are good. We would ideally like to see lots of very long-distance contacts between mates, as that is a sign of strong Hi-C signal. On the order of 1-20% is reasonable, though it depends on the assembly. For scaffolding best results are obtained when this is higher than 5%.
* Number of read pairs with mates mapping to different contigs/chromosomes. These are good if they represent contacts within a cell, but bad if they represent noise or contacts between cells (e.g. for metagenomic data). On the order of 10-40% seems standard, again it depends on particulars of the assembly.
* Number of split reads. These are good, usually, as they hopefully represent Hi-C junctions and thus successful Hi-C. There are of course other reasons why a read might be split.

## What do I do if there is a problem?
Problems observed in QC may indicate an issue either with the Hi-C reads or the assembly used for alignments. If the assembly is bad or comes from a distantly related organism or set of organisms, you should expect to see artifacts in the alignment of Hi-C reads.

* If the issue is the reads, you can try filtering your read alignments, either removing bad contigs or low-confidence reads. Our tool `matlock` has utilities for doing this.
* If the issue is the assembly, you can either get a new/more appropriate assembly somehow, or you can attempt to fix your existing assembly.
* The best ways to fix assemblies in our experience is to break up chimerically assembled contigs and/or to purge haplotigs. This can be achieved by either breaking on gaps if they exist in your assembly (e.g. runs of Ns) under the assumption that most chimerae span such gaps, or by directly inferring and breaking misjoins in your assembly. Breaking on gaps is fairly trivial, and our tool `polar_star` can help you infer and break misjoins using long read data. There are several open-source haplotig purgers worth trying and many modern assemblers come with parameters or tools to purge them as part of the assembly process.
* If the issue is that you simply don't have enough long-distance Hi-C contacts, **unavoidably you will sometimes have to remake the Hi-C library**.
