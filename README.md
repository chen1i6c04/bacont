# Installation
You will need to install the external dependencies separately with `conda`, which can be found in `requirement.yaml`
```commandline
git clone https://github.com/chen1i6c04/bacont.git
cd bacont
conda env create -f requirement.yaml
conda activate bacont
```
# Quick usage
```commandline
bacont.py -o output_dir long_reads.fastq.gz
```
* Estimate genome completeness with `BUSCO`
```commandline
bacont.py -o output_dir long_reads.fastq.gz --busco_db busco_database
```
* Filter raw reads with length and Q-score
```commandline
bacont.py -o output_dir long_reads.fastq.gz -l 1000 -q 10
```