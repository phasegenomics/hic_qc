[![Build Status](https://travis-ci.com/phasegenomics/bam_to_mate_hist.svg?branch=master)](https://travis-ci.com/phasegenomics/bam_to_mate_hist)

# `hic_qc.py` readme.

This script is intended as a simple QC method for Hi-C libraries, based on reads in a BAM file aligned to some genome/assembly. For our full recommendations on aligning and QC'ing Hi-C data, please [see here](https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html).

The most informative Hi-C reads are the ones that are long-distance contacts, or contacts between contigs of an assembly. This tool quantifies such contacts and makes plots of contact distance distributions. The most successful Hi-C libraries have many long-distance and among-contig contacts.

Hi-C connectivity drops off in approximately a power-law with increasing linear sequence distance. Consequently, one expects Hi-C reads to follow a characteristic distribution, wherein there is a spike of many read pairs at distances close to zero, which drops off smoothly (in log space) with increasing distance. If there are odd spikes or discontinuities, or if there are few long-distance contacts, there may be a problem either with the library or the assembly.

## Dependencies
* python 2.7, 3.5, or 3.6
* numpy
* pysam
* matplotlib
* pdfkit
* markdown
* wkhtmltopdf (needs to be installed manually)

For installation, run this statement in a terminal:

`git clone https://github.com/phasegenomics/bam_to_mate_hist.git && cd bam_to_mate_hist && pip install --user -r requirements.txt && python setup.py install --user`

We include a `requirements.txt` file with dependencies, which should be installed if you use the above command. However, if you want to use the PDF report feature of this tool, you will need to install `wkhtmltopdf` externally, as we cannot install this readily.

This script has been verified on MacOSX, Ubuntu Linux, and Amazon Linux.

Some dependencies such as matplotlib don't play nicely with all pythons, such that some pythons in e.g. virtualenvs may not work. In that specific case you can just deactivate the virtualenv.

### Conda installation
We have also successfully installed requirements using the following conda command (replace X and Y with desired python version):

`conda env create -n bam_to_mate_hist --file env.yml`

This will create a python 3.6 environment and automatically install the `wkhtmltopdf` dependency.

To use a different python version, run:
```
conda --add channels bioconda --add channels conda-forge
conda create -q -n test-environment python=$PYTHON_VERSION pysam numpy scipy matplotlib wkhtmltopdf pdfkit markdown
```

## Usage
In the most basic use-case, you can run the script in a terminal

`python hic_qc.py -b input.bam -n num_reads_to_use`

where `input.bam` is your BAM file from aligning Hi-C reads to your reference, and `num_reads_to_use` is just the number of reads you want to sample from the BAM file (default 1 million; assuming there are this many reads in the file).

The script will write plots in PDF format to the file `Read_mate_dist.pdf` in the working directory, unless `-o` or `--outfile` are set as described below.

The script will also quantify some basic QC metrics and print those to the screen.

The script can also make a full-on PDF report of those metrics with the plots embedded if you set the `-r`/`--make_report` flag.

To set the name of the files written out, such as the PNG figures and the report PDF, set the `-o /path/to/outfile` or `--outfile_name /path/to/outfile` parameters.

QC is performed using a set of thresholds in JSON format. By default, the file `bam_to_mate_hist/collateral/thresholds.json` is used. The chosen file may be changed with the `--thresholds` flag. Note that the thresholds in the default file are informed by Phase Genomics' analysis of thousands of Hi-C libraries, and reflect what we ourselves use for QC.

Different QC thresholds may be present in a thresholds file. The default file includes thresholds for genome scaffolding projects and metagenome deconvolution projects. The `--sample_type` argument is used to specify which set of thresholds in the thresholds file should be used in the run, and is also noted at the top of the report. By default, the `genome` sample type is used. Additional sample types may be added to a thresholds JSON file by making them keys in the file.

## Library judgement and thresholds

### Judgement categories
The report includes a judgement about the library and the assembly it was mapped to based on the observed statistics, shown at the top of the report. Libraries are given one of four classifications:
* PASS - the library and assembly appear to be sufficient for the purposes shown in the report.
* FAIL - the library and assembly appear to be insufficient for the purposes shown in the report.
* MIXED RESULTS - the library and assembly are probably sufficient for the purposes shown in the report, but there is some additional noise or other unexpected properties in the report as well.
* LOW SIGNAL - the library and assembly don't appear to be actively bad, but there is not very much observable long-range Hi-C signal. It's likely they are insufficient for the purposes shown in the report.
IMPORTANT NOTE: because input assembly is a significant contributor to the ability to perform a given analysis, a good library can still generate a failed result with a bad assembly. This can particularly occur with difficult-to-align-to assemblies, such as polyploids or highly heterozygous/repetitive assemblies.

### "Good" properties
Three statistics are used to determine if a set of alignments has "good" aspects to it:
* HQ RPs >10KB apart (CTGs >10KB): the percentage of read pairs that map with high quality (MAPQ >=20, max edit distance <=5, not dupes) in which both mates align to the same contig, the contig is at least 10kbp long, and the mates are at least 10kbp apart.
* Intercontig HQ RPs (CTGs >10KB): the percentage of read pairs that map with high quality in which each mates aligns to a different contig, and each of those contigs are at least 10kbp.
* Same strand HQ RPs: the percentage of read pairs that map with high quality to the same strand (such reads are almost guaranteed to be Hi-C junctions).
The entire set of alignments is deemed to have "good" properties if ALL THREE good criteria pass the threshold.

### "Bad" properties
Three statistics are used to determine if a set of alignments has "bad" aspects to it:
* Duplicate reads: the percentage of reads that are flagged as PCR duplicates (by a prior tool such as Picard or SAMBLASTER).
* Zero map quality reads: the percentage of reads aligning with zero MAPQ.
* Unmapped reads: the percentage of reads that could not be mapped.
The entire set of alignments is deemed to have "bad" properties if ANY OF THE THREE bad criteria fail the threshold.

### Generating the final judgement
Thresholds for these "good" and "bad" aspects come from the specified thresholds JSON file and are shown in the report. The fields in the report are highlighted for convenience to show whether a specific metric passed or failed the threshold. These are used to generate the judgement calls:
* PASS - "good" and "not bad"
* FAIL - "not good" and "bad"
* MIXED RESULTS - "good" and "bad"
* LOW SIGNAL - "not good" and "not bad"

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
Problems observed in QC may indicate an issue either with the Hi-C reads or the assembly used for alignments. If the assembly is bad or e.g. comes from a distantly related organism or set of organisms, you should expect to see artifacts in the alignment of Hi-C reads.

* If the issue is the reads, you can try filtering your read alignments, either removing bad contigs or low-confidence reads. Our tool `matlock` has utilities for doing this.
* If the issue is the assembly, you can either get a new/more appropriate assembly somehow, or you can attempt to fix your existing assembly.
* The best way to fix assemblies in our experience is to break up chimerically assembled contigs. This can be achieved by either breaking on gaps if they exist in your assembly (e.g. runs of Ns) under the assumption that most chimerae span such gaps, or by directly inferring and breaking misjoins in your assembly. Breaking on gaps is fairly trivial, and our tool `polar_star` can help you infer and break misjoins using long read data.
* If the issue is that you simply don't have enough long-distance Hi-C contacts, **unavoidably you will sometimes have to remake the Hi-C library**.
