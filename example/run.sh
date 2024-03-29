##############
## Unphased ##
##############

python3 ../qtl2-founder-haps.py individuals.vcf.gz founders.vcf.gz probs.rds --snps snps.txt

# Inspect output:
R -e 'x <- readRDS("probs.rds"); names(x); dim(x[[1]]); x[[1]][, , 1:5]'

############
## Phased ##
############

python3 ../qtl2-founder-haps.py individuals.vcf.gz founders.vcf.gz probs1.rds --snps snps.txt --working-dir tmp-qtl2-founder-haps-1 --haplotype 1
python3 ../qtl2-founder-haps.py individuals.vcf.gz founders.vcf.gz probs2.rds --snps snps.txt --working-dir tmp-qtl2-founder-haps-2 --haplotype 2

# Inspect outputs:
R -e 'x <- readRDS("probs1.rds"); names(x); dim(x[[1]]); x[[1]][, , 1:5]'
R -e 'x <- readRDS("probs2.rds"); names(x); dim(x[[1]]); x[[1]][, , 1:5]'
