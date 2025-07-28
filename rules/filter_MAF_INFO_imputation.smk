#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

##########################################################################
#  Filter imputed dataset based on MAF of reference panel and INFO score #
##########################################################################


rule merge_phased_bcfs:
    """ 
    Merge all phased samples 
    """
    input: 
        phased_vcf_annotate = expand(
            "output/GLIMPSE_phased/{sample}_phased_annotated.{chrom}.vcf.gz", 
            sample = SAMPLE, 
            allow_missing = True,
        ),
    output:
        merged_phased_vcf = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz"
    log:
        "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz"
    threads: 4
    shell:
        """
        bcftools merge \
        {input.phased_vcf_annotate} \
        -Oz -o {output.merged_phased_vcf} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_phased_vcf}
        """

rule MAF_sites_ref_pan_phased:
    """
    Extract sites from the reference panel based on selected MAF filter
    """
    input:
        ref_panel_phased="data/reference_panel_vcf/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
    output:
        ref_maf_vcf = "output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.vcf.gz",
        ref_maf_tsv = "output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.tsv.gz"
    params:
        maf=config["maf_cutoff"]
    log:
        "output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.vcf.log"
    threads: 4
    shell:
        """
        bcftools view \
        -q {params.maf}:minor \
        {input.ref_panel_phased} \
        -Oz -o {output.ref_maf_vcf} 2> {log}
        
        bcftools index -f {output.ref_maf_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_maf_vcf} | \
        bgzip -c > {output.ref_maf_tsv}

        tabix -s1 -b2 -e2 {output.ref_maf_tsv}
        """

rule maf_sites_phased_vcf:
    """
    Extract MAF sites from merged phased VCF
    """
    input:
        ref_maf_tsv = "output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.tsv.gz",
        merged_phased_vcf = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz"
    output:
        merged_phased_vcf_maf = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz"
    log:
        "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz.log"
    threads: 4
    shell:
        """
        bcftools view {input.merged_phased_vcf} \
        --regions-file {input.ref_maf_tsv} \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf}
        """

rule recalibrate_info_phased_vcf:
    """
    Re-calibrate INFO scores based on all samples present in the merged VCF
    """
    input:
        merged_phased_vcf_maf = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz"
    output:
        merged_phased_vcf_maf_recalibrated = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz"
    log:
        "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz.log"
    threads: 4
    shell:
        """
        bcftools plugin impute-info \
        {input.merged_phased_vcf_maf} \
        -Ob -o {output.merged_phased_vcf_maf_recalibrated} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf_recalibrated}
        """

rule filter_recalibrated_INFO_phased_vcf:
    """
    Filter for INFO after re-calibrating INFO score for all samples together
    """
    input:
        merged_phased_vcf_maf_recalibrated = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz"
    output:
        merged_phased_vcf_maf_recalibrated_info = "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz"
    params:
        info=config["info_cutoff"]
    log:
        "output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz.log"
    threads: 4
    shell:
        """
        bcftools view \
        {input.merged_phased_vcf_maf_recalibrated} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf_recalibrated_info} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf_recalibrated_info}
        """