# aln_stats
Generate alignment statistics from a paf file

## Dependencies

* [minimap2](https://github.com/lh3/minimap2) 
* [numpy](https://numpy.org/)

## Installation

Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

```
git clone https://github.com/microbemarsh/aln_stats.git

cd aln_stats

conda env create -n aln_stats -f aln_stats.yml

conda activate aln_stats
```

## Usage

```
# Activate the conda environment
conda activate aln_stats

# Run minimap2 and output the paf file of interest
minimap2 -x map-ont -t $nproc -c --cs=long --secondary=no tests/emu_reference.fasta tests/full_length.fa > tests/test.paf

# Run aln_stats (make sure to use the correct order of files)
aln_stats tests/test.paf tests/emu_reference.fasta tests/full_length.fa

# Bonus! Extract the aligned sequences
paftools.js longcs2seq tests/test.paf tests/full_length.fa
```
