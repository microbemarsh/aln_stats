# aln_stats
Generate alignment statistics from a paf file

## Dependencies

[minimap2](https://github.com/lh3/minimap2)
[numpy](https://numpy.org/)

## Installation

Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

```
git clone https://github.com/microbemarsh/aln_stats.git

cd aln_stats

conda create env -n aln_stats -f aln_stats.yml

conda activate aln_stats
```

## Usage

```
# Activate the conda environment
conda activate aln_stats

# Run minimap2 and output the paf file of interest
minimap2 -x map-ont -t $nproc -c --cs=long --secondary=no test/emu_reference.fasta test/full_length.fa > test/test.paf


# Bonus! Extract the aligned sequences we can then run
paftools.js longcs2seq test/test.paf test/full_length.fa
```