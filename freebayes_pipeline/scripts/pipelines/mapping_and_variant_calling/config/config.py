#!/usr/bin/env python
"""
all variables without "_dir" name ending (except freebayes_dir) or "__" name beginning will 
be saved directly in {base_name}.config.yaml
"""

import os
from socket import gethostname
from sacgf.util.server import is_on_server
from sacgf.util.file_utils import bioinformatics_dirname
from os.path import join
from os.path import dirname
import sys

user_bioinfo_dir = bioinformatics_dirname()

#Nedd to set gatk_bundle_dir and tools_dir below if NOT running on SACGF, Bigmem or Tango
# __hostname = gethostname()

if is_on_server():
    tools_dir = '/data/sacgf/reference/map_n_call_pipeline/'
    # snakemake = '/data/sacgf/admin/python_virtual_environments/python3/bin/snakemake'
    # snakemake = '/data/sacgf/admin/miniconda3/bin/snakemake'
    # snakemake = '/data/sacgf/tools/miniconda36/bin/snakemake'
    snakemake = '/data/sacgf/tools/snakemake/miniconda3/bin/snakemake'
    sacgf_env_runner = join(user_bioinfo_dir, 'scripts/bin/tango_env_runner.sh')
    # sacgf_py3_env_runner = join(bioinfo_dir, 'scripts/bin/sacgf_py3_env_runner.sh')
    #Current Nimblegen exome capturer
    nimblegen_dir = '/data/sacgf/reference/hg19/ExomeCaptureRegions/Nimblegen/Nimblegen_SeqCap_EQ_Exome_v3/'
    nimblegen_call_region_bed = nimblegen_dir + \
            'SeqCap_EZ_Exome_v3_capture.1000bp_padded_each_side.merged.Ensembl_style.cut_columns.bed'
    nimblegen_stat_region_bed = nimblegen_dir + \
            'SeqCap_EZ_Exome_v3_capture.b37.for_qualimap.bed'
    gatk_bundle_dir = '/data/sacgf/reference/b37/GATK_bundles/GATK_2.8_bundle_b37/'
# elif __is_on_tau:
else:
    tools_dir = '/tau/references/for_franks_map_n_call_pipeline/map_n_call_pipeline/'
    snakemake = '/tau/code/python_virtual_environments/python3/bin/snakemake'
    sacgf_env_runner = join(user_bioinfo_dir, 'scripts/bin/sacgf_env_runner.sh')
    # sacgf_py3_env_runner = join(bioinfo_dir, 'scripts/bin/sacgf_py3_env_runner.sh')
    gatk_bundle_dir = '/tau/references/for_franks_map_n_call_pipeline/gatk_bundle/2.8/b37/'
    #FIXME these Nimblegen files are not yet copied to TAU
    nimblegen_dir = '/data/sacgf/reference/hg19/ExomeCaptureRegions/Nimblegen/Nimblegen_SeqCap_EQ_Exome_v3/'
    nimblegen_call_region_bed = nimblegen_dir + \
            'SeqCap_EZ_Exome_v3_capture.1000bp_padded_each_side.merged.Ensembl_style.cut_columns.bed'
    nimblegen_stat_region_bed = nimblegen_dir + \
            'SeqCap_EZ_Exome_v3_capture.b37.for_qualimap.bed'
walltime = '60:00:00' #use 35 for --threads_per_sample, normally should be fine for one 40X WGS each time on Tizard

config_file_dir = dirname(os.path.abspath(__file__))
#uncomment this if passing the file to --custom_config
#config_file_dir = os.path.join(user_bioinfo_dir, 'scripts/pipelines/mapping_and_variant_calling/config') #can't end with "/"

project_dir = dirname(config_file_dir) # the dir to mapping_and_variant_calling
ruledir = join(project_dir, 'rules/')
scripts_dir = join(project_dir, 'scripts/')
bioinfo_dir = dirname(dirname(dirname(project_dir)))
dir_to_wgs_chr_beds = join(project_dir, 'config/wgs_chr_beds/') # this will be saved to {base_name}.config.yaml

cutadapt = '/data/sacgf/tools/miniconda3_tango_py3_env/envs/tango_py3_env/bin/cutadapt'

bwa = tools_dir + 'bwa/0.7.12-r1044/bwa'
bwakit_run_bwamem = tools_dir + 'bwakit/0.7.12-r1044/run-bwamem'
sambamba = tools_dir + 'sambamba/sambamba_v0.6.5'
gatk = tools_dir + 'gatk/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar'
freebayes_dir = tools_dir + 'freebayes/v1.2.0/'
snpsift = tools_dir + 'snpeff/4.1L/SnpSift.jar'
bcftools = tools_dir + 'bcftools/1.5-32/bcftools'
igvtools = tools_dir + 'igvtools/2.3.81/igvtools'
vt = tools_dir + 'vt/v0.5772-60f436c3/vt'
bedtools = tools_dir + 'bedtools/2.17/bin/bedtools'
bgzip = tools_dir + 'htslib/1.3-8/bgzip'
tabix = tools_dir + 'htslib/1.3-8/tabix'
qualimap = tools_dir + 'qualimap/2.2.1/qualimap'
java8= tools_dir + 'java/oracle_jre1.8.0_144/bin/java'

remove_star_non_variants_in_vcf_script  = join(scripts_dir, 'remove_star_non_variants_in_vcf.sh')
combine_regional_depth_script  = join(scripts_dir, 'combine_regional_depth.py')
combine_qualimap_script  = join(scripts_dir, 'combine_qualimap_genome_results.py')
trim_amp_sam_script = join(scripts_dir, 'trim_amp_sam.py')
gatk_dna_hard_filters_script = join(scripts_dir, 'gatk_dna_hard_filters.py')
combine_flagstat_script = join(scripts_dir, 'combine_flagstat.py')
combine_depth_of_coverage_script = join(scripts_dir, 'combine_depth_of_coverage.py')
combine_key_bam_stats_script = join(scripts_dir, 'combine_key_bam_stats.py')
freebayes_bed_generate_regions_script = join(scripts_dir, 'freebayes_bed_generate_regions.py')
remove_an_intermediate_file_script = join(scripts_dir, 'remove_an_intermediate_file.sh')
sex_prediction_script = join(scripts_dir, 'determine_sex_with_x_variants.py')
calc_homozygous_variant_rate_script = join(scripts_dir, 'homozygosity_rate_each_autosome.py')
b37_region_for_freebayes = config_file_dir + '/b37.freebayes_parallel.fasta.fai'
freebayes_vcf_header = config_file_dir + '/freebayes_vcf_header.v1.2.0.txt'

#always use the annotation scripts in user's BIOINFO_DIR
var_annot_script_dir = join(user_bioinfo_dir, 'scripts/pipelines/variant_annotation_and_filtration/')
variant_level_annotate_script = join(var_annot_script_dir, 'variant_level_annotation.py')
gene_and_sample_level_annotationy_script = join(var_annot_script_dir, 'gene_and_sample_level_annotation.py')
standard_filter_script = join(var_annot_script_dir, 'standard_filters.py')

#gatk bundle files
ref = gatk_bundle_dir + 'human_g1k_v37_decoy.fasta'
hapmap = gatk_bundle_dir + 'hapmap_3.3.b37.vcf'
dbsnp = gatk_bundle_dir + 'dbsnp_138.b37.vcf'
omni = gatk_bundle_dir + '1000G_omni2.5.b37.vcf'
g1k_snps = gatk_bundle_dir + '1000G_phase1.snps.high_confidence.b37.vcf'
g1k_indels = gatk_bundle_dir + '1000G_phase1.indels.b37.vcf'
mills_indels = gatk_bundle_dir + 'Mills_and_1000G_gold_standard.indels.b37.vcf'

#the gvcf of the selected exomes from the 1KG project
extra_gvcf_dir = tools_dir + 'extra_gvcfs/v0/'

#default params for running commands and for the workflow
sample_names = []
sample_to_units = {}
sample_to_bam = {}
sample_to_analysis_ready_bam = {}
sample_to_analysis_ready_bam_index = {}
passed_gvcfs = []
platform = None
stop_after = None
skips = []
add_chris_gene_list = False
species = 'HUMAN'
call_region_bed = None
stat_region_bed = None


# #may be use for depth <200X
# freebayes_call_param = (
            # '--pooled-continuous '
            # '--haplotype-length 0 '
            # '--no-population-priors '
            # '--hwe-priors-off '
            # '--binomial-obs-priors-off --allele-balance-priors-off '
            # '--min-alternate-fraction 0.005 --min-alternate-count 2 '
            # '--min-mapping-quality 20 --min-base-quality 20 '
            # '--read-max-mismatch-fraction 0.05 '
            # '--reference-quality 60,40 '
            # '--use-best-n-alleles 6 '
            # '--min-alternate-qsum 50 '
            # '--min-supporting-mapping-qsum 50 '
            # '--no-complex '
            # '--strict-vcf '
            # )

freebayes_call_param = (
            '--pooled-continuous '
            '--haplotype-length 0 '
            '--no-population-priors '
            '--hwe-priors-off '
            '--binomial-obs-priors-off --allele-balance-priors-off '
            '--min-alternate-fraction 0.005 --min-alternate-count 3 '
            '--min-mapping-quality 20 --min-base-quality 20 '
            '--read-max-mismatch-fraction 0.03 '  #Was 0.05, before 14/9/2018
            '--reference-quality 60,40 '
            '--use-best-n-alleles 4 ' #Was 6, before 14/9/2018
            '--min-alternate-qsum 70 '
            '--min-supporting-mapping-qsum 70 '
            '--no-complex '
            '--strict-vcf '
            )

freebayes_filter_param = 'SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0' #capture-based data; In most of the case, it is better not excluding RPR > 0 & RPL > 0
freebayes_filter_param_for_amplicon = None #For amplicon-based NGS data

gatk_hc_param = ''
gatk_hc_param_for_amplicon = '--min_mapping_quality_score 30 --min_base_quality_score 30'
gatk_gtGVCF_param = '' #'--standard_min_confidence_threshold_for_emitting was deprecated
gatk_gtGVCF_param_for_amplicon = '' #'--standard_min_confidence_threshold_for_calling was deprecated
# double quotes must be used inside the filter statement below
gatk_snp_filters = (
            '--filterExpression "QD < 2.0" --filterName "QDFilter" '
            '--filterExpression "MQ < 40.0" --filterName "MQFilter" '
            '--filterExpression "FS > 60.0" --filterName "FSFilter" '
            '--filterExpression "SOR > 4.0" --filterName "SORFilter" '
            '--filterExpression "MQRankSum < -12.5" --filterName "MQRSFilter" '
            '--filterExpression "ReadPosRankSum < -8.0" --filterName "RPRSFilter" ')
gatk_indel_filters = (
            '--filterExpression "QD < 2.0" --filterName "QDFilter" '
            '--filterExpression "ReadPosRankSum < -20.0" --filterName "RPRSFilter" '
            '--filterExpression "FS > 200.0" --filterName "FSFilter" '
            '--filterExpression "SOR > 10.0" --filterName "SORFilter" ' )
gatk_snp_filters_for_amplicon = (
            '--filterExpression "QD < 6.0" --filterName "QDFilter" '
            '--filterExpression "MQ < 40.0" --filterName "MQFilter" '
            '--filterExpression "MQRankSum < -12.5" --filterName "MQRSFilter" '
            )
gatk_indel_filters_for_amplicon = (
            '--filterExpression "QD < 6.0" --filterName "QDFilter" '
            '--filterExpression "MQ < 30.0" --filterName "MQFilter" '
            )
