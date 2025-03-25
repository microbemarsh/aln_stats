# aln_stats
Generate alignment statistics from alignment files (sam, bam or paf files)

## Dependencies

* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://www.htslib.org/)
* [numpy](https://numpy.org/)
* [pysam](http://pysam.readthedocs.io/en/latest/)

## Installation

Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

```
# Clone the respository
git clone https://github.com/microbemarsh/aln_stats.git

# Go to repo
cd aln_stats

# Create conda environment
conda env create -n aln_stats -f aln_stats.yml

# Activate conda environment
conda activate aln_stats
```

## Usage

### SAMs ###
```
# Activate the conda environment
conda activate aln_stats

# Run minimap2 and output the paf file of interest
minimap2 -ax map-ont -t $nproc --cs=long --secondary=no tests/emu_reference.fasta tests/full_length.fa > tests/test.sam

# Run aln_stats (make sure to use the correct order of files)
aln_stats tests/test.sam tests/emu_reference.fasta tests/full_length.fa

# Extract the aligned sequences
paftools.js longcs2seq tests/test.paf tests/full_length.fasamtools view -F 4 tests/test.sam | samtools fasta -o tests/aligned_sequences.fa
```

### BAMs ###
```
# Activate the conda environment
conda activate aln_stats

# Run minimap2 and output the paf file of interest
minimap2 -ax map-ont -t $nproc --cs=long --secondary=no tests/emu_reference.fasta tests/full_length.fa | samtools view -b -o tests/test.bam

# Run aln_stats (make sure to use the correct order of files)
aln_stats tests/test.bam tests/emu_reference.fasta tests/full_length.fa

# Extract the aligned sequences
samtools view -F 4 tests/test.bam | samtools fasta -o tests/aligned_sequences.fa
```

### PAFs ###
```
# Activate the conda environment
conda activate aln_stats

# Run minimap2 and output the paf file of interest
minimap2 -x map-ont -t $nproc -c --cs=long --secondary=no tests/emu_reference.fasta tests/full_length.fa > tests/test.paf

# Run aln_stats (make sure to use the correct order of files)
aln_stats tests/test.paf tests/emu_reference.fasta tests/full_length.fa

# Extract the aligned sequences
paftools.js longcs2seq tests/test.paf tests/full_length.fa
```
