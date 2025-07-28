#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

###########################################
#  Download reference panel and BAM files #
###########################################

rule get_filt_reference_panel:
    """
    Download filtered reference panel, ready to use
    """
    output:
        vcf = "data/reference_panel_vcf/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        tbi = "data/reference_panel_vcf/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tbi",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/d1p5Gd2PaB/ref-panel_{wildcards.chrom}_sample-snp_filltags_filter.phased.vcf.gz > {output.vcf}
        wget --quiet -O - https://sid.erda.dk/share_redirect/d1p5Gd2PaB/ref-panel_{wildcards.chrom}_sample-snp_filltags_filter.phased.vcf.gz.tbi > {output.tbi}
        """

rule get_ref_fasta:
    """
    Download canfam3.1 reference fasta file
    """
    output:
        fa = "data/reference_fasta/canFam3_withY.fa",
        fai = "data/reference_fasta/canFam3_withY.fa.fai",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_fasta/canFam3_withY.fa > {output.fa}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_fasta/canFam3_withY.fa.fai > {output.fai}
        """