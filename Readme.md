# `bam_to_mate_hist.py` readme.

This script is intended as a simple QC method for Hi-C libraries, based on reads in a BAM file aligned to some genome/assembly. 

Hi-C connectivity drops off in approximately a power-law with increasing linear sequence distance. Consequently, one expects Hi-C reads to follow a characteristic distribution, wherein there is a spike of many read pairs at distances close to zero, which drops off smoothly (in log space) with increasing distance. If there are odd spikes or discontinuities, there may be a problem either with the library or the assembly.

## Dependencies
* python2
* numpy
* pysam
* matplotlib

These can be installed with pip, e.g.:

`pip install pysam`

I've tested this script on MacOSX/Linux. 

Some dependencies such as matplotlib don't play nicely with all pythons, such that some pythons in e.g. virtualenvs may not work. In this case you can just deactivate the virtualenv. 

## Usage
No installation should be required, simply clone the repo:
`git clone [repo_name] && cd [repo_name]`

You can run the script in a terminal
`$ python bam_to_mate_hist.py -b input.bam -n num_reads_to_use`
where `input.bam` is your BAM file from aligning Hi-C reads to your reference, and `num_reads_to_use` is just the number of reads you want to sample from the BAM file (default 1 million; assuming there are this many reads in the file). 

The script will write plots in PDF format to the file `Read_mate_dist.pdf` in the working directory.