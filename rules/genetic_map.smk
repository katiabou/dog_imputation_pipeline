#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

#######################################################
# Rules for creating a GLIMPSE compatible genetic map #
#######################################################

global CHROM


rule download_canine_genetic_map:
    """
    Download the canFam3.1 genetic map
    """
    output:
        average=temp(expand("data/gen_map/{chrom}_average_canFam3.1.txt", chrom=CHROM)),
        male=temp(expand("data/gen_map/{chrom}_male_canFam3.1.txt", chrom=CHROM)),
        female=temp(expand("data/gen_map/{chrom}_female_canFam3.1.txt", chrom=CHROM)),
        tar=temp("data/gen_map/dog_genetic_maps.b38.tar.gz"),
    params:
        path="data/gen_map/",
    shell:
        "wget --quiet -O {output.tar} -o /dev/null https://github.com/cflerin/dog_recombination/raw/refs/heads/master/dog_genetic_maps.tar.gz && "
        "tar -xzf {output.tar} -C {params.path}"


rule glimpse_genetic_map:
    """
    Format the canFam3.1 genetic map for GLIMPSE
    """
    input:
        "data/gen_map/{chrom}_average_canFam3.1.txt",
    output:
        "data/gen_map/{chrom}_average_canFam3.1_modified.txt",
    shell:
        """
        awk '{{ print $2" "$1" "$4 }}' {input} > {output}
        """

