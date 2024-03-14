##########################################################################
#  Filter imputed dataset based on MAF of reference panel and INFO score #
##########################################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

rule merge_phased_bcfs:
    """ 
    Merge all phased samples 
    """
    input: 
        phased_vcf_annotate = expand('output/GLIMPSE_phased/{bam_imputation}_phased_annotated.{chrom}.vcf.gz', bam_imputation = BAM, allow_missing = True)
    output:
        merged_phased_vcf = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz'
    log:
        'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz'
    threads: 10
    shell:
        '''
        bcftools merge \
        {input.phased_vcf_annotate} \
        -Oz -o {output.merged_phased_vcf} \
        --threads {threads}

        tabix -p vcf {output.merged_phased_vcf}
        '''

rule MAF_sites_ref_pan_phased:
    """
    Extract sites using a MAF filter from the reference panel 
    """
    input:
        ref_panel_phased = config['reference_panel']
    output:
        ref_maf_vcf = 'output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.vcf.gz',
        ref_maf_tsv = 'output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.tsv.gz'
    params:
        maf=config['maf_cutoff']
    log:
        'output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.vcf.log'
    shell:
        '''
        bcftools view \
        -q {params.maf}:minor \
        {input.ref_panel_phased} \
        -Oz -o {output.ref_maf_vcf} 2> {log}
        
        bcftools index -f {output.ref_maf_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_maf_vcf} | \
        bgzip -c > {output.ref_maf_tsv}

        tabix -s1 -b2 -e2 {output.ref_maf_tsv}
        '''

rule maf_sites_phased_vcf:
    """
    Extract MAF sites from merged phased VCF
    """
    input:
        ref_maf_tsv = 'output/reference_panel/ref_panel_{chrom}_MAF_{maf_cutoff}.phased.tsv.gz',
        merged_phased_vcf = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}.vcf.gz'
    output:
        merged_phased_vcf_maf = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    log:
        'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz.log'
    threads: 10
    shell:
        '''
        bcftools view {input.merged_phased_vcf} \
        --regions-file {input.ref_maf_tsv} \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf}
        '''


rule filter_INFO_phased_vcf:
    """
    Filter for INFO WITHOUT re-calibrating INFO score for all samples together
    """
    input:
        merged_phased_vcf_maf = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    output:
        merged_phased_vcf_maf_info = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info_cutoff}.vcf.gz'
    params:
        info=config['info_cutoff']
    log:
        'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info_cutoff}.vcf.gz.log'
    threads: 10
    shell:
        '''
        bcftools view \
        {input.merged_phased_vcf_maf} \
        --trim-alt-alleles -Ou | \
        bcftools view \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf_info}

        tabix -p vcf {output.merged_phased_vcf_maf_info}
        '''

rule recalibrate_info_phased_vcf:
    """
    Re-calibrate INFO scores based on all samples present in the merged VCF
    """
    input:
        merged_phased_vcf_maf = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    output:
        merged_phased_vcf_maf_recalibrated = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz'
    log:
        'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz.log'
    threads: 10
    shell:
        '''
        bcftools plugin impute-info \
        {input.merged_phased_vcf_maf} \
        -Ob -o {output.merged_phased_vcf_maf_recalibrated} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf_recalibrated}
        '''

rule filter_recalibrated_INFO_phased_vcf:
    """
    Filter for INFO after re-calibrating INFO score for all samples together
    """
    input:
        merged_phased_vcf_maf_recalibrated = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz'
    output:
        merged_phased_vcf_maf_recalibrated_info = 'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz'
    params:
        info=config['info_cutoff']
    log:
        'output/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info_cutoff}.vcf.gz.log'
    threads: 10
    shell:
        '''
        bcftools view \
        {input.merged_phased_vcf_maf_recalibrated} \
        --trim-alt-alleles -Ou | \
        bcftools view \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf_recalibrated_info}

        tabix -p vcf {output.merged_phased_vcf_maf_recalibrated_info}
        '''