#!/usr/bin/env python

"""
When running on local machines, this script relies on the shell environment variable $SACGF_DIR
"""

import sys
sys.dont_write_bytecode = True
import os
from sacgf.util.server import get_sacgf_dir
from sacgf.util.file_utils import bioinformatics_dirname

# locate SACGF_DIR and BIOINFO_DIR
# Shell environment variable SACGF_DIR is needed on local machines
BIOINFO_DIR = bioinformatics_dirname()
SACGF_DIR = get_sacgf_dir()
if SACGF_DIR is None:
    sys.exit('Error: environment variable $SACGF_DIR is needed, when not running on SACGF, Tizard, or Bigmem.')

DATE_VERSION = '1/Oct/2018'
VERSIONS = {
    'DATE_VERSION': DATE_VERSION,
    'SnpEff': '4.1L with data file GRCh37.75',
    'COSMIC': '79',
    'CancerGeneCensus': 'COSMIC V79',
    'dbSNP': '149',
    'SIFT': 'dbNSFP V2.6',
    'Polyphen2': 'dbNSFP V2.6',
    'MutationTaster': 'dbNSFP V2.6',
    'FATHMM': 'dbNSFP V2.6',
    'miRNA_bind_site': 'hg19 version from the UCSC Genome Browser; downloaded January 2014',
    'phyloP': 'hg19 version from the UCSC Genome Browser; downloaded August 2014',
    'GERP': 'hg19 version from the UCSC Genome Browser; downloaded August 2013',
    'The 1000 Genome project': 'Phase 3 V5b',
    'UK10K': 'version 20160215',
    'ExAC': 'version 0.3.1',
    'gnomAD': 'r2.0.2',
    'Exome Sequencing Project': 'ESP6500SI-V2',
    'TGI Tier': 'hg19 version; downloaded August 2014 from the official website',
    'CADD': 'Version 1.3',
    'Segmental duplication': 'hg19/b37 version; downloaded December 2014 from the official website',
    'PFAM': 'hg19 version from the UCSC Genome Browser; downloaded September 2014',
    'INTERPRO': 'hg19 version from the UCSC Genome Browser; downloaded September 2014',
    }

VAR_ANNOT_FILE_DIR = os.path.join(SACGF_DIR, 'reference/variant_annotation/')
##################################################################
# TOOLS
VAR_ANNOT_TOOL_DIR = os.path.join(VAR_ANNOT_FILE_DIR, 'tools/')

# SNPEFF and its data file
SNPEFF_DIR = os.path.join(VAR_ANNOT_TOOL_DIR, 'snpEff_4.1L/') #TODO test and upgrade it to version 4.3i (all files are already available on SACGF)
SNPEFF_JAR  = os.path.join(SNPEFF_DIR, 'snpEff.jar')
SNPEFF_DATA = 'GRCh37.75'
SNPSIFT_JAR  = os.path.join(SNPEFF_DIR, 'SnpSift.jar')
SNPEFF_HEADERS_FILE = os.path.join(VAR_ANNOT_FILE_DIR, 'snpeff_headers.v3.txt')
SNPEFF_ALL_PREDICTS_HEADER_FILE = os.path.join(VAR_ANNOT_FILE_DIR, 'snpeff_headers.all_predicts.v3.txt')
COSMIC_SNPEFF_HEADERS_FILE = os.path.join(VAR_ANNOT_FILE_DIR, 'cosmic_snpeff_headers.v2.txt')

BCFTOOLS = os.path.join(VAR_ANNOT_TOOL_DIR, 'bcftools_1.3.1-207/bcftools') #upgraded 10/1/2017

VT = os.path.join(VAR_ANNOT_TOOL_DIR, 'vt_0.5772-60f436c3/vt')

GATK_JAR = os.path.join(VAR_ANNOT_TOOL_DIR, 'GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar')

BGZIP_CMD = os.path.join(VAR_ANNOT_TOOL_DIR, 'htslib-1.3.2-180/bgzip')

TABIX_CMD = os.path.join(VAR_ANNOT_TOOL_DIR, 'htslib-1.3.2-180/tabix')

IGVTOOLS_CMD = os.path.join(VAR_ANNOT_TOOL_DIR, 'igvtools_2.3.40/igvtools')

### scripts

# use a script from bioinformatics repository to perform MARK_LISTED_GENES_IN_CSV with Chris' Haematopoietic gene list
# TODO use pandas 
MARK_LISTED_GENES_IN_CSV_CMD =  os.path.join(BIOINFO_DIR,  
    'scripts/project_specific/finding_causative_variant/mark_listed_genes_in_csv.py')

###################################################################
###########################################
# External databases and data files
###########################################
VAR_ANNOT_DB_DIR = os.path.join(VAR_ANNOT_FILE_DIR, 'b37/')
GENE_ANNOT_DB_DIR = os.path.join(SACGF_DIR, 'reference/gene_annotation/current/')

B37_REF = os.path.join(VAR_ANNOT_DB_DIR, 'ref_genome/human_g1k_v37_decoy.fasta')
B37_GENOME_SIZE_FILE = os.path.join(VAR_ANNOT_FILE_DIR, 'genome_info_data/b37.chrom.sizes')

# miRNA binding site predition DB file. It is modified from the one downloaded by Paul
MIRNA_BINDING_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'mirna/targetScan.Jan_2014.unoverlapping.bed.gz')
MIRNA_BINDING_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'mirna/header.txt')

# the combined GERP and phyloP BED file
GERP_n_PHYLOP_BED_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'gerp_n_phyloP/gerp_atLeast2.phyloP_atLeast1.bed.gz')

# allele frequencies and others from the public population variant datasets
POPULATION_FRE_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'population_variant_datasets/merged.1kg_uk10k_exac_esp.tsv.gz')
POPULATION_FRE_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'population_variant_datasets/header.1kg_uk10k_exac_esp.txt')
GNOMAD_FRE_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'gnomAD/r2.0.2/gnomAD_combined.tsv.gz')
GNOMAD_FRE_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'gnomAD/r2.0.2/header.gnomAD_combined.tsv.gz.txt')

#COSMIC
cosmic_dir = os.path.join(VAR_ANNOT_DB_DIR, 'cosmic/')
COSMIC_FILE = os.path.join(cosmic_dir, 'CosmicCodingMuts.v79.tsv.gz')
COSMIC_HEADER = os.path.join(cosmic_dir, 'cosmic_header.txt')

NON_SPECIFIC_COSMIC_DB = os.path.join(cosmic_dir,  'non_specific_records.v79.pydict') # version is needed here
COSMIC_TRANSCRIPTS_FILE = os.path.join(cosmic_dir, 'transcripts_in_COSMIC_v75.Ensembl_v75.txt')

#clinVar
CLINVAR_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'clinvar/current_tsv') 
CLINVAR_FILE = os.path.realpath(CLINVAR_FILE) # convert symbolic link to real file path
CLINVAR_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'clinvar/clinvar_header.txt') 

#TGI Tier
TGI_TIER_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'tgi_tier/TGI_Tiers_combined.bed.gz')
TGI_TIER_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'tgi_tier/tgi_tier.header')

#CADD
CADD_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'cadd/1.3/at_least_10.tsv.gz')
CADD_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'cadd/header.txt')

#AA-Effects
AA_EFF_PREDICT_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'dbNSFP/2.6/extracted_prediction.v2.tsv.gz')
AA_EFF_PREDICT_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'dbNSFP/2.6/header.v1.txt')

#LCR low complexity region
LCR_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'lcr/LCR.edited.bed.gz')
LCR_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'lcr/header.txt')

#Segmental duplications
SEGMENT_DUP_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'segment_dup/sd.edited.sorted.merged.tsv.gz')
SEGMENT_DUP_HEADER = os.path.join(VAR_ANNOT_DB_DIR, 'segment_dup/header.txt')

#dbSNP, Only for rs-ID annotation
DBSNP_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'dbsnp/dbsnp_149.normalized.uniq.tsv.gz')

# protein domains, downloaded by Paul
PFAM_BED_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'pfam/pfam.Sep_2014.bed')
INTERPRO_BED_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'interpro/interpro.Sep_2014.bed')

#branch point bed file
BRANCH_POINT_BED_FILE = os.path.join(VAR_ANNOT_DB_DIR, 'branch_points/branch_points.bed')

#gene level annotation files
ENSEMBL_KEY_FULL_ANNOTATION_FILE = os.path.join(GENE_ANNOT_DB_DIR, 'ensembl_key.v75.full_annotation.tsv')
ENSEMBL_TO_HGNC_MAPPING_FILE = os.path.join(GENE_ANNOT_DB_DIR, 'ensem_id_to_hgnc_mapping.pydict')
ENSEMBL_CANONICAL_TRANSCRIPTS_FILE = os.path.join(GENE_ANNOT_DB_DIR, 
                                                '../materials/human_canonical_transcripts_v75.pyset')

CHRIS_ONE_PLUS_LIST  = os.path.join(SACGF_DIR, 'reference/gene_lists/for_programming/chris_HM_one_or_more.06-01-2016.tsv')
CHRIS_TWO_PLUS_LIST  = os.path.join(SACGF_DIR, 'reference/gene_lists/for_programming/chris_HM_two_or_more.06-01-2016.tsv')

