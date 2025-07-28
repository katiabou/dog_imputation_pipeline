#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

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

#Bams which will be imputed (replace names and bam file paths to your samples)
samples_df = pd.read_table(
    "sample_lists/bams.tsv", dtype=str,delimiter="\t",
).set_index("Sample", drop=False)
SAMPLE = list(samples_df["Sample"])
BAM = list(samples_df["Bam"])

##### Rules to include #####
include: "rules/genetic_map.smk"
include: "rules/fetch_data.smk"
include: "rules/GLIMPSE_impute.smk"
include: "rules/filter_MAF_INFO_imputation.smk"

##### Target rules #####
rule all:
    input:
        expand(
            ["output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz"
            ],
            chrom=CHROM, 
            maf_cutoff=config["maf_cutoff"],
            info_cutoff=config["info_cutoff"]
        )