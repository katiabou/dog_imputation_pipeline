import pandas as pd

##### load config file #####
configfile: 'config.yaml'

##### Define total chromosome numbers ##### 
CHROM = [f"chr{i}" for i in range(1,int(config['chromosome_number'])+1)]

##### set wildcard constraints ##### 
wildcard_constraints:
    chrom="chr\d+",
    info_cutoff="\d+.\d+",
    maf_cutoff="\d+.\d+"

##### Bam files which will be imputed #####
bams_df = pd.read_table(config['bam_imputation'], dtype=str, delimiter="\t").set_index("Sample", drop=False)
BAM = list(bams_df['Sample'])

##### Paths to programmes #####
glimpse_chunk = config['glimpse_chunk']
glimpse_impute = config['glimpse_impute']
glimpse_ligate = config['glimpse_ligate']
glimpse_sample = config['glimpse_sample']
glimpse_concordance = config['glimpse_concordance']

##### Rules to include #####
#include: "rules/ref_panel.smk"
include: "rules/GLIMPSE_impute.smk"
include: "rules/filter_MAF_INFO_imputation.smk"

##### Target rules #####
rule all:
    input:
        expand('output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info_cutoff}.vcf.gz', chrom=CHROM, maf_cutoff=config['maf_cutoff'], info_cutoff=config['info_cutoff']),
        expand('output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz', chrom=CHROM, maf_cutoff=config['maf_cutoff'], info_cutoff=config['info_cutoff'])
