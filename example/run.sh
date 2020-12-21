python3 ../qtl2-founder-haps.py individuals.vcf.gz founders.vcf.gz probs.rds --snps snps.txt

# Inspect output:
R -e 'x <- readRDS("probs.rds"); names(x); dim(x[[1]]); x[[1]][, , 1:5]'
