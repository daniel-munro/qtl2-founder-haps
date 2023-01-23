import argparse
import os
import pysam
import pandas as pd
from bisect import bisect_left
import subprocess


def genotype_code(gt: tuple, founder: bool = False) -> str:
    if gt == (None, None):
        return "-"
    elif gt == (0, 0):
        return "A"
    elif (gt[0] == 0 and gt[1] > 0) or (gt[0] > 0 and gt[1] == 0):
        return "-" if founder else "H"
    elif gt[0] > 0 and gt[1] > 0:
        return "B"
    else:
        raise ValueError(f"GT not recognized: {gt}")


def genetic_pos(chrmap: pd.DataFrame, pos: int) -> float:
    r = bisect_left(chrmap["pos"], pos)
    if r == len(chrmap["pos"]):
        return chrmap["cm"][r - 1]
    elif chrmap["pos"][r] == pos or r == 0:
        return chrmap["cm"][r]
    else:
        # Interpolate the genetic position.
        p_lo = chrmap["pos"][r - 1]
        p_hi = chrmap["pos"][r]
        g_lo = chrmap["cm"][r - 1]
        g_hi = chrmap["cm"][r]
        rel = (pos - p_lo) / (p_hi - p_lo)
        return g_lo + rel * (g_hi - g_lo)


def make_qtl_inputs(args):
    if args.snps is None:
        IDs = None
    else:
        IDs = set(open(args.snps, "r").read().splitlines())

    maps = {}
    for chrom in range(1, 21):
        filename = os.path.join(args.gmap_dir, f"MAP4chr{chrom}.txt.gz")
        maps[chrom] = pd.read_table(filename, sep=" ", names=["pos", "ratio", "cm"])

    vcf = pysam.VariantFile(args.individuals)
    samples = list(vcf.header.samples)
    genos = {}
    refs = {}
    ID_list = []
    for rec in vcf.fetch():
        ID = rec.id if rec.id is not None else f"{rec.contig}:{rec.pos}"
        if IDs is None or ID in IDs:
            gt = [rec.samples[sample]["GT"] for sample in samples]
            if args.haplotype in {1, 2}:
                gt = [(g[args.haplotype - 1], g[args.haplotype - 1]) for g in gt]
            labels = [genotype_code(g, founder=False) for g in gt]
            genos[ID] = labels
            refs[ID] = rec.ref
            ID_list.append(ID)

    # ID_list = [x for x in ID_list if x in genos.keys()]
    IDs = set(ID_list)

    vcf = pysam.VariantFile(args.founders)
    strains = list(vcf.header.samples)
    founder_genos = {}
    ref_mismatch = 0
    # remove = set()
    ID_list = []
    for rec in vcf.fetch():
        ID = rec.id if rec.id is not None else f"{rec.contig}:{rec.pos}"
        if ID in IDs:
            gt = [rec.samples[strain]["GT"] for strain in strains]
            labels = [genotype_code(g, founder=True) for g in gt]
            # assert rec.ref == refs[ID]
            if rec.ref != refs[ID]:
                ref_mismatch += 1
                # remove.add(ID)
                # del genos[ID]
                # del founder_genos[ID]
            else:
                founder_genos[ID] = labels
                ID_list.append(ID)


    if ref_mismatch > 0:
        print(f"{ref_mismatch} SNPs removed due to reference mismatch.")
        # ID_list = [ID for ID in ID_list if ID not in remove]

    if not os.path.exists(args.working_dir):
        os.makedirs(args.working_dir)
    with open(os.path.join(args.working_dir, "geno.csv"), "w") as out:
        out.write(f"id,{','.join(samples)}\n")
        for ID in ID_list:
            out.write(f"{ID},{','.join(genos[ID])}\n")

    with open(os.path.join(args.working_dir, "founder_geno.csv"), "w") as out:
        out.write(f"id,{','.join(strains)}\n")
        for ID in ID_list:
            out.write(f"{ID},{','.join(founder_genos[ID])}\n")

    with open(os.path.join(args.working_dir, "pmap.csv"), "w") as out:
        out.write("marker,chr,pos\n")
        for ID in ID_list:
            chrom, pos = tuple(ID.replace("chr", "").split(":"))
            pos = int(pos) / 1e6  # Units are Mbp.
            out.write(f"{ID},{chrom},{pos}\n")

    with open(os.path.join(args.working_dir, "gmap.csv"), "w") as out:
        out.write("marker,chr,pos\n")
        for ID in ID_list:
            chrom, pos = tuple(ID.replace("chr", "").split(":"))
            gpos = genetic_pos(maps[int(chrom)], int(pos))
            out.write(f"{ID},{chrom},{round(gpos, 6)}\n")

    with open(os.path.join(args.working_dir, "covar.csv"), "w") as out:
        out.write("id,generations\n")
        for sample in samples:
            out.write(f"{sample},90\n")

    strain_str = ", ".join([f'"{strain}"' for strain in strains])
    cntrl_command = (
        'qtl2::write_control_file('
        'output_file = "control.yaml", '
        'overwrite = TRUE, '
        'crosstype = "hs", '
        'geno_file = "geno.csv", '
        'founder_geno_file = "founder_geno.csv", '
        'gmap_file = "gmap.csv", '
        'pmap_file = "pmap.csv", '
        'covar_file = "covar.csv", '
        'crossinfo_covar = "generations", '
        'geno_codes = c(A = 1L, H = 2L, B = 3L), '
        f'alleles = c({strain_str}), '
        'na.strings = "-", '
        'geno_transposed = TRUE, '
        'founder_geno_transposed = TRUE'
        ')'
    )
    subprocess.run(f"cd {args.working_dir} && R -e '{cntrl_command}'", shell=True)

gmaps = os.path.join(os.path.dirname(__file__), "genetic_map")

p = argparse.ArgumentParser(description="Wrapper for R/qtl2 to calculate founder haplotype probabilities")
p.add_argument("individuals", help="path to VCF file for individuals")
p.add_argument("founders", help="path to VCF file for founder strains")
p.add_argument("out", help="Name of 3D array output file (*.rds, serialized R object)")
p.add_argument("--snps", "-s", help="File of SNP IDs to subset VCFs (e.g. to include observed and not imputed SNPs)")
p.add_argument("--gmap-dir", default=gmaps, help="Directory containing genetic mapping files")
p.add_argument("--working-dir", default="tmp-qtl2-founder-haps", help="Name of directory to write qtl2 input files")
# p.add_argument("--save-interm", action="store_true", help="With this flag, qtl2 input files will be saved")
p.add_argument("--founder-pairs", action="store_true", help="Output probabilities per founder pair instead of collapsing to per-founder")
p.add_argument("--haplotype", type=int, default=0, help="If set to 1 or 2, VCFs are assumed to be phased and the output will reflect only the first or second haplotypes in the VCF")
p.add_argument("--cores", type=int, default=1, help="Number of cores to use when calculating probabilities")
args = p.parse_args()

make_qtl_inputs(args)
qtl_command = (
    'library(qtl2); '
    f'cross <- read_cross2("{os.path.join(args.working_dir, "control.yaml")}"); '
    f'pr <- calc_genoprob(cross, error_prob = 0.01, cores = {args.cores}); '
)
if not args.founder_pairs:
    qtl_command += f'pr <- genoprob_to_alleleprob(pr); '
qtl_command += f'saveRDS(pr, "{args.out}")'
subprocess.run(f"R -e '{qtl_command}'", shell=True)
