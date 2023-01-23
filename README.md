# qtl2-founder-haps.py

Wrapper for [R/qtl2](https://kbroman.org/qtl2/) to calculate founder haplotype probabilities

The qtl2 R package can be used to infer founder haplotype probabilities (or "genotype probabilities") at a set of loci for a multi-parent population using a hidden Markov model. However, if the genotypes are in VCF format, preprocessing is required. This Python script can be run from the command line to get a 3D array of founder haplotype probabilities for a set of individuals using a VCF file for the individuals and a VCF file for the founders. It prepares input files, runs qtl2, and saves the result.

By default, the found haplotype inference does not consider phase, and the probabilities pertain to the combination of the two genome haplotypes. That is, at each locus in each individual, probabilities are calculated over all possible unordered pairs of founder haplotypes (including homozygous pairs), and then collapsed into one probability per founder. This script provides an option to get two sets of probabilities, one for each haplotype in the phased chromosomes, using a sort of hack. Assuming the VCF file contains phased genotypes, pass `--haplotype 1` or `--haplotype 2`, and each genotype will be converted to a homozygous genotype with either the first or second allele. The collapsed output probabilities should then reflect only the first or second haplotype of all chromosomes. See `run.sh` for a usage example.

Currently it works only for HS rats at around generation 90. If you would like expanded options, please let me know.

## Dependencies

- R
- qtl2 R package
- Python 3.6+
- pandas
- pysam

## Usage

```
usage: qtl2-founder-haps.py [-h] [--snps SNPS] [--gmap-dir GMAP_DIR] [--working-dir WORKING_DIR]
                            [--founder-pairs] [--haplotype HAPLOTYPE] [--cores CORES]
                            individuals founders out

Wrapper for R/qtl2 to calculate founder haplotype probabilities

positional arguments:
  individuals           path to VCF file for individuals
  founders              path to VCF file for founder strains
  out                   Name of 3D array output file (*.rds, serialized R object)

optional arguments:
  -h, --help            show this help message and exit
  --snps SNPS, -s SNPS  File of SNP IDs to subset VCFs (e.g. to include observed and not imputed
                        SNPs)
  --gmap-dir GMAP_DIR   Directory containing genetic mapping files
  --working-dir WORKING_DIR
                        Name of directory to write qtl2 input files
  --founder-pairs       Output probabilities per founder pair instead of collapsing to per-founder
  --haplotype HAPLOTYPE
                        If set to 1 or 2, VCFs are assumed to be phased and the output will reflect
                        only the first or second haplotypes in the VCF
  --cores CORES         Number of cores to use when calculating probabilities
```

## TODO

- Option to save output as text file
- Option to save or delete intermediate files
- Handle additional complex cross designs and generation numbers
- Check for and fix reference allele mismatches between VCFs
- Option to get probabilities for regularly spaced pseudomarkers
