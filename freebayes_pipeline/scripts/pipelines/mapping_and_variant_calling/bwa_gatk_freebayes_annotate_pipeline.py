#!/usr/bin/env python

"""

"""

import argparse
from socket import gethostname
import imp
import os
from os.path import join
import pysam
import sys
import types
import yaml
import pandas as pd
from numpy import nan as Nan
from sys import exit
from sacgf.util.file_utils import absolute_file_path, random_file_name, bioinformatics_dirname
from sacgf.util.sam_utils import get_sample_list, get_platform_list
from sacgf.util.subprocess_utils import run_and_handle_returncode
from sacgf.util.vcf_utils import get_sample_list as get_sample_list_in_vcf

sys.dont_write_bytecode = True

ILLUMINA = 'ILLUMINA'
Data_type = ['WGS', 'WES', 'Panel', 'Amplicon-Panel']
WGS, WES, Panel, Amplicon = [t.lower() for t in Data_type]

def handle_cl_options():
    parser = argparse.ArgumentParser(description=('Pipeline script for running '
        'BWA-GATK-Freebayes-annotate as described in the Wiki-page.'))
    parser.add_argument('--base_name', required=True, help='The project name, batch name, or '
            'any string you want to name the final output file')
    parser.add_argument('--data_type', choices=Data_type, required=True,
            help=('The data type. WGS: whole genome sequencing; WES: whole exome sequencing; '
                'Panel: capture-based gene panel sequencing; Amplicon-Panel: amplicon-based gene panel sequencing'))
    parser.add_argument('--out_dir', dest='workdir', required=True,
            help='Output diretory. It will be created if not existing ')
    parser.add_argument('--fastq_dirs', nargs='*', default=[], metavar='FASTQ_DIR',
            help=('Space separated list of directories containing FASTQ files. All FASTQ files '
                'in the listed directories will be processed unless '
                '--include_samples or --exclude_samples is used'))
    parser.add_argument('--fastq_sample_sheets', nargs='*', default=[],
            metavar='SAMPLE_SHEET',
            help="Space separated list of SACGF's sample "
                "sheets of Illumina machine, in the corresponding order to the --fastq_dirs argument")
    parser.add_argument('--use_random_read_group_info',
            #metavar='PLATFORM',
            choices='ILLUMINA IONTORRENT SOLID HELICOS ONT PACBIO CAPILLARY LS454'.split(),
            help="Use this option to write random and proper read group info into BAM files, "
            "when fastq_info_spreadsheet is unavailable. "
            "You need to specify the sequencing platform")
    parser.add_argument('--adapter_fasta',
            help="Pass the adapter fasta file to activate adapter trimming. In this file, "
                "the first read is for read 1 and the second read is for read 2. "
                "Currently, adapter trimming is only for paired FASTQ files")
    parser.add_argument('--inbam_dirs', nargs='*', default=[], metavar='BAM_DIR',
            help=('Space separated list of directories containing BAM files to be passed to mapping. '
                'All BAM files in the listed directories will be processed unless '
                '--include_samples or --exclude_samples is used'))
    parser.add_argument('--inbams', nargs='*', default=[], metavar='BAM',
            help=('Space separated list of BAM files to be passed to mapping'))
    parser.add_argument('--analysis_ready_bam_dirs', nargs='*', default=[], metavar='ANALYSIS_READY_BAM_DIR',
            help=('Space separated list of directories containing analysis-ready BAM files to '
                'be passed to variant calling. All BAM files in the listed directories will be '
                'processed unless --include_samples or --exclude_samples is used'))
    parser.add_argument('--analysis_ready_bams', nargs='*', default=[], metavar='ANALYSIS_READY_BAM',
            help=('Space separated list of analysis-ready BAM files to '
                'be passed to variant calling.'))
    parser.add_argument('--include_samples', nargs='*', default=[], metavar='SAMPLE',
            help=('Space separated list of samples IDs to process. Case sensitive.'))
    parser.add_argument('--exclude_samples', nargs='*', default=[], metavar='SAMPLE',
            help=('Space separated list of samples IDs to exclude. Case sensitive. '))
    parser.add_argument('--gvcf_dirs', nargs='*', default=[], metavar='GVCF_DIR',
            help=("Space separated list of directories containing GATK's .g.vcf.gz files to be "
                "passed to GATK's genotyping. All .g.vcf.gz files in the listed directories will be "
                "processed"))
    parser.add_argument('--gvcfs', nargs='*', default=[], metavar='GVCF',
            help=("Space separated list of GATK's .g.vcf.gz files to be "
                "passed to GATK's genotyping "))
    parser.add_argument('--caller', dest='callers', choices=['GATK', 'Freebayes', 'both'],
            default='GATK',
            help=('The variant caller to be used. "both" means using both GATK and Freebayes. '
                '(default: %(default)s)'))
    parser.add_argument('--stop_after', choices=['GATK_BQSR', 'GATK_HC'],
            help=('The point after which to stop in the workflow' ))
    parser.add_argument('--skip', dest='skips', action='append',
            choices=['Marking_Duplicates', 'GATK_BQSR', 'BAM_Stats', 'Qualimap', 'Variant_Stats'],
            help='A step to skip. Can be specified multiple times to skip multiple steps.')
    parser.add_argument('--call_region_bed', help=('The BED file to feed to the variant caller, to '
        'limit the regions for variants calling '))
    parser.add_argument('--depthofcoverage_bed', help=("The BED file to feed to GATK's depthofcoverage"))
    parser.add_argument('--qualimap_bed', help=("The BED file to feed to Qualimap"))
    parser.add_argument('--amp_inserts_bed', help=('The amplicon insert regions bed. '
        'This bed file is also used as CALL_REGION_BED and DEPTHOFCOVERAGE_BED if the later two not specified'))
    parser.add_argument('--amp_regions_bed', help=('The amplicon regions bed. '
        'This bed file is also used as QUALIMAP_BED if the later not specified'))
    parser.add_argument('--max_threads', help=('The max total number of cpu cores/threads allowed to '
        'be required by the running and waiting jobs. (default: %(default)s)'), type=int, default=120)
    parser.add_argument('--threads_per_sample', help=('The number of cpu cores/threads used to run '
        'the command (e.g. BWA) on a single sample. (default: %(default)s)'), type=int, default=16)
    parser.add_argument('--threads_on_aggregate_samples', help=("The number of cpu cores/threads used "
        "to run the command (e.g. Freebayes' variant calling) with data of all samples as input. "
        "(default: %(default)s)"),
        type=int, default=24)
    parser.add_argument('--max_io_heavy_jobs', help='The max number of running and waiting jobs '
            '(e.g. merging BAMs), which can cause heavy disk I/O. '
            '(default: %(default)s)',
            type=int, default=48)
    parser.add_argument('--email', help='Specify your email if you want an email notification after '
            'the pipeline is finished or failed')
    parser.add_argument('--not_annotate', action='store_true',
            help="Do not perform variant or gene level annotation.")
    parser.add_argument('--not_gene_annotate', action='store_true',
            help="Do variant level but not gene level annotation. No variant spreadsheets will be generated.")
    parser.add_argument('--add_chris_gene_list', action='store_true',
            help="Use Chris Hahn's leukaemia gene list during the gene level annotation")
    parser.add_argument('--add_VAF_column', action='store_true',
            help="Force to add VAF column into the variant spreadsheet for GATK variants.")
    parser.add_argument('--spreadsheet_format', choices=['xlsx', 'tsv'], default='xlsx',
            help="The format of spreadsheet file of variants. (default: %(default)s)")
    parser.add_argument('--trio',
            help=('A comma seperate list of the trio information for inferring inheritance model '
                'that a (germline) variant is fitting. Use this option only when '
                'the samples are a trio. '
                'The format is "sex-of-child,child,mother,father" '
                '(e.g. "male,HSS333,HSS222,HSS111"), where sex is either "female" or "male". '
                'It only works with variant calls by GATK.'))
    parser.add_argument('--just_get_bam_stats', action='store_true',
            help=('Use this option to produce BAM stats, with --analysis_ready_bams and/or --analysis_ready_bam_dirs'))
    #to pass custom parameters to GATK and Freebayes
    parser.add_argument('--force_GATK_hard_filter', action='store_true',
            help=('Use this option to force to perform GATK hard filter.'))
    parser.add_argument('--force_extra_joint_gt', action='store_true',
            help=('When the input files are gVCF created on ILLUMINA WES data '
                'and the number of samples is small (e.g. <30), '
                'use this option to force joint-genotyping with the extra gVCF files.'))
    parser.add_argument('--freebayes_call_param',
            help=('Paramaters passed to Freebayes to call variants. Using this option will overwrite '
            'the default params used in the pipeline. See FAQ in the Wiki-page for the default params'))
    parser.add_argument('--freebayes_filter_param',
            help=('Paramaters passed to the vcffilter of Freebayes to filter variants. Using this '
                'option will overwrite the default params used in the pipeline. See FAQ in the '
                'Wiki-page for the default params'))
    parser.add_argument('--gatk_hc_param',
            help=('Paramaters passed to GATK HaplotypeCaller. Using this option will overwrite '
            'the default params used in the pipeline. See FAQ in the Wiki-page for the default params'))
    parser.add_argument('--gatk_gtGVCF_param',
            help=('Paramaters passed to GATK genotypeGVCFs. Using this option will overwrite '
            'the default params used in the pipeline. See FAQ in the Wiki-page for the default params'))
    parser.add_argument('--gatk_snp_filters',
            help=('Paramaters passed to GATK hard filtering SNPs. Using this option will overwrite '
            'the default params used in the pipeline. See FAQ in the Wiki-page for the default params'))
    parser.add_argument('--gatk_indel_filters',
            help=('Paramaters passed to GATK hard filtering INDELs. Using this option will overwrite '
            'the default params used in the pipeline. See FAQ in the Wiki-page for the default params'))
    parser.add_argument('--custom_config',
            help="For debug purpose. Pass a custom config file to execute the workflow.")
    parser.add_argument('--ccb_queue', action='store_true',
            help=('Use this option to force to submit jobs onto the CCB queue'))
    parser.add_argument('--non_scratch_dir_to_write_log_files',
            help=('Writes snakemake config and SLUM log files to this dir, '
                'when --ccb_queue is set and the --out_dir is at /scratch '
                'It will be created if not existing '))
    parser.add_argument('--dryrun', action='store_true',
            help="For debug purpose. Create Snakemake files in the OUT_DIR, "
                "print the Snakemake command, but not execute the workflow.")
    parser.add_argument('--keep_temp_files', action='store_true',
            help="For debug purpose. "
                "Do not remove intermediate files.")
    args = parser.parse_args()

    ##validate args
    in_type_number = 0
    if args.fastq_dirs:
        in_type_number += 1
    if args.inbam_dirs or args.inbams:
        in_type_number += 1
    if (args.analysis_ready_bam_dirs or args.analysis_ready_bams) and \
            not (args.gvcf_dirs or args.gvcfs):
        in_type_number += 1
    if args.gvcf_dirs or args.gvcfs:
        in_type_number += 1
    if in_type_number == 0:
        parser.error('No input files are specified. Please specify FASTQ, BAM, ANALYSIS_READY_BAM or GVCF files')
    elif in_type_number > 1:
        parser.error('Only one type of input file can be specified: FASTQ, BAM, ANALYSIS_READY_BAM or GVCF')

    if args.fastq_dirs:
        if args.fastq_sample_sheets and args.use_random_read_group_info:
            parser.error('These two options can not be used simultaneously: '
                    '--fastq_sample_sheets and --use_random_read_group_info')
        elif not (args.fastq_sample_sheets or args.use_random_read_group_info):
            parser.error('One of these two options must be used, when FASTQ files are the input: '
                    '--fastq_sample_sheets or --use_random_read_group_info')
        #check if Marking_Duplicates is in skips for Amplicon-Panel data
        if args.data_type == 'Amplicon-Panel' and 'Marking_Duplicates' not in args.skips:
            parser.error('For Amplicon-Panel data, please remember to skip Marking_Duplicates via '
                    'the --skips option. Contact Frank for more details.')
    if args.include_samples and args.exclude_samples:
        parser.error('These two options can not be used simultaneously: '
                '--include_samples and --exclude_samples')
    # if args.gvcf_dirs or args.gvcfs:
        # if not (args.analysis_ready_bam_dirs or args.analysis_ready_bams):
            # parser.error('--analysis_ready_bam_dirs and/or --analysis_ready_bams are needed '
                    # 'for GATK ReadBackedPhasing when '
                    # 'user-spedified .g.vcf.gz files are passed. Note '
                    # 'only BAM files for samples listed in .g.vcf will be used, and '
                    # '--include_samples and --exclude_samples will be ignored if specified.')

    if args.fastq_dirs and args.fastq_sample_sheets and \
            len(args.fastq_dirs) != len(args.fastq_sample_sheets):
        parser.error('The number of listed items are not the same between --fastq_dirs and --fastq_sample_sheets')
    if args.fastq_sample_sheets:
        not_ending_with_csv = []
        for f in args.fastq_sample_sheets:
            if not f.endswith('.csv'):
                not_ending_with_csv.append(f)
        if not_ending_with_csv:
            print 'The file name of fastq_info_spreadsheet does not end with ".csv":'
            for f in not_ending_with_csv:
                print '\t%s' % f
            exit(1)

    max_threads, threads_per_sample, threads_on_aggregate_samples = args.max_threads, \
            args.threads_per_sample, args.threads_on_aggregate_samples
    # if not (max_threads >= threads_on_aggregate_samples >= threads_per_sample):
        # parser.error('It requires max_threads >= threads_on_aggregate_samples >= threads_per_sample')

    #lower these strings
    args.data_type = args.data_type.lower()
    args.callers = args.callers.lower()
    args.stop_after = args.stop_after.lower() if args.stop_after else None
    args.skips = [askip.lower() for askip in args.skips] if args.skips else []

    args.callers = [args.callers] if args.callers != 'both' else ['gatk', 'freebayes']

    if args.data_type == WGS:
        args.call_region_bed = bioinformatics_dirname() + '/scripts/pipelines/mapping_and_variant_calling/config/wgs_chr_beds/b37.chrom.bed'

    #change files or dirs with absolute_file_path()
    args.fastq_dirs = [absolute_file_path(p) for p in args.fastq_dirs]
    args.fastq_sample_sheets = [absolute_file_path(p) for p in args.fastq_sample_sheets]
    args.inbam_dirs = [absolute_file_path(p) for p in args.inbam_dirs]
    args.inbams = [absolute_file_path(p) for p in args.inbams]
    args.analysis_ready_bam_dirs = [absolute_file_path(p) for p in args.analysis_ready_bam_dirs]
    args.analysis_ready_bams = [absolute_file_path(p) for p in args.analysis_ready_bams]
    args.adapter_fasta = absolute_file_path(args.adapter_fasta)
    args.gvcf_dirs = [absolute_file_path(p) for p in args.gvcf_dirs]
    args.gvcfs = [absolute_file_path(p) for p in args.gvcfs]
    args.workdir = absolute_file_path(args.workdir)
    args.call_region_bed = absolute_file_path(args.call_region_bed)
    # args.stat_region_bed = absolute_file_path(args.stat_region_bed)
    args.depthofcoverage_bed = absolute_file_path(args.depthofcoverage_bed)
    args.qualimap_bed = absolute_file_path(args.qualimap_bed)
    args.amp_inserts_bed = absolute_file_path(args.amp_inserts_bed)
    args.amp_regions_bed = absolute_file_path(args.amp_regions_bed)
    args.custom_config = absolute_file_path(args.custom_config)
    args.non_scratch_dir_to_write_log_files = absolute_file_path(args.non_scratch_dir_to_write_log_files)

    # depthofcoverage_bed = qualimap_bed

    if args.freebayes_call_param is None:
        del args.freebayes_call_param
    if args.freebayes_filter_param is None:
        del args.freebayes_filter_param
    if args.gatk_hc_param is None:
        del args.gatk_hc_param
    if args.gatk_gtGVCF_param is None:
        del args.gatk_gtGVCF_param
    if args.gatk_snp_filters is None:
        del args.gatk_snp_filters
    if args.gatk_indel_filters is None:
        del args.gatk_indel_filters

    if args.fastq_dirs and args.data_type == Amplicon and args.amp_inserts_bed is None:
            parser.error('--amp_inserts_bed is needed for %s data type' % Amplicon)

    return args

def parse_adapter_fasta_file(adapter_fasta):
    r1_adapter, r2_adapter = None, None
    if adapter_fasta is not None:
        lines = open(adapter_fasta).readlines()
        for i, line in enumerate(lines):
            if line.startswith('>'):
                if r1_adapter is None:
                    r1_adapter = lines[i+1].strip()
                elif r2_adapter is None:
                    r2_adapter = lines[i+1].strip()
                if r2_adapter is not None:
                    break
    return r1_adapter, r2_adapter

def is_illumina_convention(fastq):
    if fastq.endswith('_R1_001.fastq.gz') or fastq.endswith('_R2_001.fastq.gz'):
        return True
    return False

def is_simplified_single_end_convention(fastq):
    if not is_illumina_convention(fastq):
        if not (fastq.endswith('_R1.fastq.gz') or fastq.endswith('_R2.fastq.gz')):
            return True
    return False

def extract_sm_from_fastq(fastq):
    """
    It tries to fit the name to Illumina convensions first, and then simplified conventions.
    """

    if is_illumina_convention(fastq):
        # <SampleName>_S1_L001_R2_001.fastq.gz
        return '_'.join(fastq.split('_')[:-4])
    else:
        if is_simplified_single_end_convention(fastq):
            return fastq[:-len('.fastq.gz')]
        else:
            return fastq[:-len('_R1.fastq.gz')]

def is_new_sheet(sheet):
    first_line = open(sheet).readline()
    return first_line.startswith( '[Header]' )

def  count_num_of_lines_to_skip(sheet):
    with open(sheet) as source:
        c = 0
        for line in source:
            if not line.startswith('[Data]'):
                c += 1
            else:
                break
    assert c > 0, '[Data] was not found in sample sheet %s' % sheet
    return c + 1

def convert_sheet_to_df(sheet):
    """
    columns in returned df:
        sample_id, sample_name, flowcell_id, lane, barcode, date, platform
    """

    platform = ILLUMINA #for now, it is always ILLUMINA

    true_new_sheet = True if is_new_sheet(sheet) else False

    if true_new_sheet:
        # 160106_NB501008_0008_AHT33JBGXX.csv
        number_of_lines_to_skip = count_num_of_lines_to_skip(sheet)
        sample_id = 'Sample_ID'
    else: # old version sheet e.g.
        # ${BIOINFO}/scripts/pipelines/mapping_and_variant_calling/test/sample_sheets/150521_SN1101_0170_AC69B8ACXX.csv
        number_of_lines_to_skip = 0
        sample_id = 'SampleID'

    file_name_parts = os.path.basename(sheet)[: -len('.csv')].split('_')

    #extract date from sheet file name
    date_on_file = file_name_parts[0]
    if len(date_on_file) != 6:
        print 'The date (%s) indicated in the file name is not in format YYMMDD: ' % date_on_file
        print '\t%s' % sheet
        exit(1)
    yy, mm, dd = date_on_file[:2], date_on_file[2:4], date_on_file[4:]
    if int(mm) > 12 or int(dd) > 31:
        print 'The date indicated in the file name not in format YYMMDD: '
        print '\t%s' % sheet
        exit(1)
    date = '20{yy}-{mm}-{dd}T00:00:00+0930'.format(yy=yy, mm=mm, dd=dd)

    flow_cell_id = extract_flow_cell_id_from_sheet(sheet)

    orig_df = pd.read_csv(sheet, skiprows=number_of_lines_to_skip, comment='#',
            skip_blank_lines=True, skipinitialspace=True, dtype=str)
    orig_columns = orig_df.columns

    # Because Dave keep changing it, either Index or index can be in the sheet. It is hard to predit what he used
    if 'Index' in orig_columns:
        index = 'Index'
    elif 'index' in orig_columns:
        index = 'index'
    else:
        exit('Index or index column is not in the sample sheet')

    # sample_id, sample_name, flowcell_id, lane, barcode, date, platform
    df = pd.DataFrame()
    df['sample_id'] = orig_df[sample_id]
    df['flowcell_id'] = flow_cell_id
    df['barcode'] = orig_df[index]
    df['date'] = date
    df['platform'] = platform
    if true_new_sheet:
        df['lane'] = Nan #lane is not in the new sheet; will extract it from the fastq file name
        if 'Sample_Name' in orig_columns:
            df['sample_name'] = orig_df['Sample_Name']
            columns_to_strip = 'sample_id sample_name flowcell_id barcode date platform'.split()
        else:
            df['sample_name'] = Nan #column Sample_Name sometimes is not there, if Dave removes it
            columns_to_strip = 'sample_id flowcell_id barcode date platform'.split()
    else:
        df['sample_name'] = Nan #no sample_name colum in the old sheet
        df['lane'] = orig_df.Lane
        columns_to_strip = 'sample_id flowcell_id lane barcode date platform'.split()

    for column in columns_to_strip:
        df[column] = df[column].str.strip()

    return df

def extract_flow_cell_id_from_sheet(sheet):
    if is_new_sheet(sheet):
        #extract flow_cell_id from the Description line
        Description = 'Description'
        with open(sheet) as source:
            for line in source:
                if Description in line:
                    flow_cell_id = line.strip().replace(Description, '').replace(',', '')
                    if not flow_cell_id:
                        exit('flow_cell_id is missing in the "Description" line in sample sheet %s' % sheet)
                    break
    else: #old sheet
        orig_df = pd.read_csv(sheet, comment='#', skip_blank_lines=True, skipinitialspace=True, dtype=str)
        flow_cell_id_list = orig_df.FCID.tolist()
        flow_cell_id_list = [_id.strip() for _id in flow_cell_id_list]
        flow_cell_id_list = list(set(flow_cell_id_list))
        if len(flow_cell_id_list) > 1:
            exit('Error: more than one flow_cell_id are in sample sheet %s' % sheet)
        elif len(flow_cell_id_list) == 0:
            exit('Error: no flow_cell_id is in sample sheet %s' % sheet)
        else:
            flow_cell_id = flow_cell_id_list[0]
    return flow_cell_id

def extract_lane_from_fastq_name(fastq):
    """
    returns the lane number if ILLUMINA convention, or returns None
    """
    if is_illumina_convention(fastq):
        # It always has lane number, if Illumina convention
        tail_length = len('001_R1_001.fastq.gz')
        if not fastq[-(tail_length+1)].startswith('L'):
            exit('Lane number is missing the the FASTQ file with ILLUMINA naming convension: %s' % fastq)
        return str(int(fastq[-tail_length: -(tail_length-3)]))
    else:
        return None

def get_fastq_input_config(fastq_dirs, fastq_sample_sheets, include_samples, exclude_samples,
        use_random_read_group_info):
    """
    FASTQ files must be named in either Illumina or simplified convension.
    Samples in fastq_dirs must be a subset of those in fastq_sample_sheets if
    fastq_sample_sheets is used

    Illumina convention (bcl2fastq V2.17):
        paired-end:
            <SampleName>_S1_L001_R1_001.fastq.gz
            <SampleName>_S1_L001_R2_001.fastq.gz
        single-end:
            <SampleName>_S1_L001_R1_001.fastq.gz

    Simplified convention:
        paired-end:
            <SampleName>_R1.fastq.gz
            <SampleName>_R2.fastq.gz
        single-end:
            <SampleName>.fastq.gz
    """

    sample_names = set()
    dir_to_fastqs = {}
    for dir_ in fastq_dirs:
        dir_to_fastqs[dir_] = []
        files = os.listdir(dir_)
        for f in files:
            if f.endswith('.fastq.gz'):
                sm = extract_sm_from_fastq(f)
                if (include_samples and sm in include_samples) or \
                        (exclude_samples and sm not in exclude_samples) or \
                        (not (include_samples or exclude_samples)):
                    sample_names.add(sm)
                    dir_to_fastqs[dir_].append(f)
    sample_names = sorted(list(sample_names))

    dir_to_sample_df = {}
    if fastq_sample_sheets:
        for dir_, sheet in zip(fastq_dirs, fastq_sample_sheets):
            dir_to_sample_df[dir_] = convert_sheet_to_df(sheet)
    else:
        for dir_ in fastq_dirs:
            dir_to_sample_df[dir_] = None

    sample_to_units = {}
    unit_to_fqs = {}
    unit_to_read_group = {}
    for dir_ in fastq_dirs:
        fastqs = dir_to_fastqs[dir_]
        if fastqs:
            while fastqs:
                for fastq in fastqs:
                    sm = extract_sm_from_fastq(fastq)
                    if sm not in sample_to_units:
                        sample_to_units[sm] = []
                    unit = '{sm}_U{num}'.format(sm=sm, num=len(sample_to_units[sm]))
                    sample_to_units[sm].append(unit)

                    if is_illumina_convention(fastq):
                        r1_end = '_R1_001.fastq.gz'
                        r2_end = '_R2_001.fastq.gz'
                        end, mate_end = (r1_end, r2_end) if fastq.endswith(r1_end) else (r2_end, r1_end)
                        mate_fastq = fastq.replace(end, mate_end)
                        if mate_fastq in fastqs:
                            fq_for_the_unit = [fastq, mate_fastq] if fastq.endswith(r1_end) else \
                                    [mate_fastq, fastq]
                        else:
                            fq_for_the_unit = [fastq]
                    else:
                        if is_simplified_single_end_convention(fastq):
                            fq_for_the_unit = [fastq]
                        else:
                            r1_end = '_R1.fastq.gz'
                            r2_end = '_R2.fastq.gz'
                            end, mate_end = (r1_end, r2_end) if fastq.endswith(r1_end) else (r2_end, r1_end)
                            mate_fastq = fastq.replace(end, mate_end)
                            fq_for_the_unit = [fastq, mate_fastq] if fastq.endswith(r1_end) else \
                                    [mate_fastq, fastq]
                    unit_to_fqs[unit] = [os.path.join(dir_, f) for f in fq_for_the_unit]
                    for f in fq_for_the_unit: #remove files which are already assigned to a unit
                        fastqs.remove(f)

                    sample_df = dir_to_sample_df[dir_]
                    if sample_df is not None: #write complete RG info
                        # sample_id, sample_name, flowcell_id, lane, barcode, date, platform
                        # use the sample_name column rather than sample_id, when both exising, as the bcl2fastq script does
                        if not sample_df.sample_name.dropna().empty:
                            df_for_sm = sample_df[sample_df.sample_name == sm]
                        else:
                            df_for_sm = sample_df[sample_df.sample_id == sm]
                        #check if sm is in the sample name or sample ID column
                        if df_for_sm.empty:
                            sys.exit('Error: %s is not in the sample name or sample ID column in the sample sheet.' % sm )
                        #It always has lane with Illumina fastq file naming convention, or it is None
                        lane_in_file = extract_lane_from_fastq_name(fastq)
                        if lane_in_file and df_for_sm.lane.any(): #ILLUMINA naming fastq with old sheet
                            df_for_sm = df_for_sm[df_for_sm.lane == lane_in_file]
                        flowcell = '&'.join(df_for_sm.flowcell_id.drop_duplicates()) #Just in case; Only one flowcell_id can be in one sheet
                        if lane_in_file:
                            lane = lane_in_file
                        elif not df_for_sm.lane.dropna().empty: #lane in old sheet
                            lane = '&'.join(df_for_sm.lane.drop_duplicates())
                        else:
                            lane = None
                        barcode =  '&'.join(df_for_sm.barcode.drop_duplicates())
                        date = df_for_sm.date.tolist()[0]
                        platform = df_for_sm.platform.tolist()[0]
                        lane = 'LN' + lane if lane else ''
                        id_ = '.'.join([n for n in [flowcell, lane, sm] if n])
                        pu = '.'.join([n for n in [flowcell, lane, barcode] if n])
                        unit_to_read_group[unit] = \
                                ('@RG\\tID:{id_}\\tSM:{sample}\\tPL:{platform}'
                                '\\tPU:{pu}\\tDT:{date}\\tCN:SACGF').format(
                                id_ = id_, pu = pu, platform = platform,
                                sample = sm, date = date)
                        # # test
                        # import ipdb
                        # ipdb.set_trace()
                    else: #write random RG info
                        unit_to_read_group[unit] = \
                                '@RG\\tID:random_{rdm}\\tSM:{sample}\\tPL:{platform}'.format(
                                rdm=random_file_name()[:6], sample=sm, platform=use_random_read_group_info)
    return {'sample_names': sample_names,
            'unit_names': sorted(unit_to_fqs.keys()),
            'sample_to_units': sample_to_units,
            'units': unit_to_fqs,
            'read_groups': unit_to_read_group,
            'platform': ILLUMINA if fastq_sample_sheets else use_random_read_group_info
            }

def get_the_sample_id_from_a_bam(abam):
    samples = get_sample_list(abam)
    if len(samples) == 0:
        exit('\nRead group or SM tag is missing in BAM file: %s' % abam)
    if len(set(samples)) > 1:
        exit('\nMultiple and different samples defined in the read group in BAM file: %s' % abam)
    return samples[0]

def make_rg_from_a_bam(abam):
    sm = get_the_sample_id_from_a_bam(abam)

    header = pysam.AlignmentFile(abam, 'rb').header #@UndefinedVariable
    rgs = header.get('RG', None)
    if rgs is None:
        exit('\nRead group is missing in BAM file: %s' % abam)
    else:
        pls = []
        ids = []
        for rg in rgs:
            pl = rg.get('PL', None)
            if pl:
                pls.append(pl)
            id_ = rg.get('ID', None)
            if id_:
                ids.append(id_)
        pl_set = set(pls)
        if len(pl_set) == 0:
            exit('\nPL is mssing in the read group in BAM file: %s' % abam)
        elif len(pl_set) > 1:
            exit('\nMultiple and different PLs defined in the read group in BAM file: %s' % abam)
        else:
            out_pl = list(pl_set)[0]
        id_set = set(ids)
        out_id = sm + '-' + '&'.join(sorted(list(id_set)))
        #same tags as with use_random_read_group_info, except for the DS tag
        return '@RG\\tID:{id_}\\tSM:{sm}\\tPL:{pl}\\tDS:Remapped by SACGF with BAM input'.format(
                    id_=out_id, sm=sm, pl=out_pl)

def get_inbam_input_config(inbams, inbam_dirs, include_samples, exclude_samples):
    all_bams_in_dirs = []
    for dir_ in inbam_dirs:
        all_bams_in_dirs.extend([os.path.join(dir_, f) for f in os.listdir(dir_) if f.endswith('.bam')])

    #apply include_samples and exclude_samples
    if include_samples:
        needed_bams_from_dirs = [f for f in all_bams_in_dirs if get_the_sample_id_from_a_bam(f) in include_samples]
    elif exclude_samples:
        needed_bams_from_dirs = [f for f in all_bams_in_dirs if get_the_sample_id_from_a_bam(f) not in exclude_samples]
    else:
        needed_bams_from_dirs = all_bams_in_dirs

    needed_bams = needed_bams_from_dirs + inbams

    sample_to_bam = {}
    sample_to_read_group = {}
    platforms = []
    sms = []
    for abam in needed_bams:
        sm = get_the_sample_id_from_a_bam(abam)
        sms.append(sm)
        sample_to_bam[sm] = abam
        sample_to_read_group[sm] = make_rg_from_a_bam(abam)
        platforms.extend(get_platform_list(abam))

    #check if a sample occurs in more than one BAM file
    for sm in set(sms):
        if sms.count(sm) > 1:
            sys.exit('Error: Sample %s occurs in more than one BAM files' % sm)

    sample_names = sample_to_bam.keys()
    platforms = list(set(platforms))
    platform = platforms[0] if len(platforms) == 1 else None

    return {'sample_names': sample_names,
            'sample_to_bam': sample_to_bam,
            'read_groups': sample_to_read_group,
            'platform': platform
            }

def find_bam_index(sample_to_bam):
    """find either the samtools or the picard style index for the user-specified bam files"""
    sample_to_index = {}
    for sm, abam in sample_to_bam.iteritems():
        samtools_index = abam + '.bai'
        if os.path.exists(samtools_index):
            sample_to_index[sm] = samtools_index
            continue
        else:
            picard_index = abam[:-1] + 'i'
            if os.path.exists(picard_index):
                sample_to_index[sm] = picard_index
            continue
        exit('\nBAM index file is missing for the user-specified bam %s' % abam)
    return sample_to_index

def bed_checker(args):
    """
    simple checker to try to make sure bed file passed to the script satisfying the requirement
    checking bed for qualimap in the snake file
    """
    #TODO re-construct a bed for internal use, rather than throwing a erro message
    #QUALIMAP_BED specific checking (six columns requiement) is in the snake file

    all_auto_chrs = set([n for n in range(1,23)])

    for bed_file in args.call_region_bed, \
            args.depthofcoverage_bed, args.qualimap_bed, args.amp_inserts_bed, args.amp_regions_bed:
        if bed_file:
            prev_chr = 1
            with open(bed_file) as inhandle:
                for line in inhandle:
                    infos = line.rstrip().split('\t')
                    if len(infos) == 1:
                        if len(infos[0]) == 0:
                            exit('\nError: No empty line (the last line?) is allowed in bed file %s' % bed_file)
                        else:
                            #TODO to re-construct an bed for actual use so that it will allow the track info lines in the input bed file
                            exit(('\nError with bed file: At least three columns are required in each line. '
                                'Bed file %s has the line shown below (track info line?) violating the '
                                'requirement: \n\n%s' % (bed_file, line)))
                    elif len(infos) >= 3:
                        chrom_name, start, end = infos[:3]
                        if chrom_name.startswith('chr'):
                            #TODO to re-construct an bed with hg19 to b37 conversion for actual
                            exit('\nError with bed file: hg19 convention (contig name starting with "chr") is NOT allowed in %s. See FAQ in our Wiki page for help' % bed_file)
                        else:
                            try:
                                chr_ = int(chrom_name)
                                if chr_ in all_auto_chrs and  chr_ < prev_chr:
                                    #TODO to re-construct an karyotypic sorted bed for actual use
                                    exit(('\nError with bed file: '
                                            'Chromosomes must be listed in the karyotypic order '
                                            'in %s. See FAQ in our Wiki page for help' % bed_file ))
                                prev_chr = chr_
                            except ValueError:
                                pass


if __name__ == '__main__':
    args = handle_cl_options()

    bed_checker(args)

    config_out = {}

    #get stuff from default_config
    # all variables without "_dir" name ending (except freebayes_dir) or "__" name beginning will 
    # be saved directly in {base_name}.config.yaml
    if not args.custom_config:
        from config import config as default_config
    else:
        default_config = imp.load_source('default_config', args.custom_config)
    for k, v in vars(default_config).iteritems():
        if k.startswith('__') or \
                (k.endswith('_dir') and k not in ['freebayes_dir', 'extra_gvcf_dir']) or \
                type(v) in [types.ModuleType, types.FunctionType]:
            pass
        else:
            config_out[k] = v

    #handle Amplicon data
    if args.data_type == Amplicon:
        config_out['freebayes_filter_param'] = config_out['freebayes_filter_param_for_amplicon']
        config_out['gatk_hc_param'] = config_out['gatk_hc_param_for_amplicon']
        config_out['gatk_gtGVCF_param'] = config_out['gatk_gtGVCF_param_for_amplicon']
        config_out['gatk_snp_filters'] = config_out['gatk_snp_filters_for_amplicon']
        config_out['gatk_indel_filters'] = config_out['gatk_indel_filters_for_amplicon']

    #keys to remove
    keys_to_remove = ('gethostname '
                    'freebayes_filter_param_for_amplicon gatk_hc_param_for_amplicon '
                    'gatk_gtGVCF_param_for_amplicon gatk_snp_filters_for_amplicon '
                    'gatk_indel_filters_for_amplicon').split()
    for k in keys_to_remove:
        del config_out[k]

    ##get stuff from args

    #handle input files: fastq, bam, analysis_ready_bam
    #one bam for one sample
    if args.fastq_dirs:
        sm_n_fastq_stuff = get_fastq_input_config(args.fastq_dirs, args.fastq_sample_sheets,
                    args.include_samples, args.exclude_samples, args.use_random_read_group_info)
        config_out.update(sm_n_fastq_stuff)
    elif args.inbam_dirs or args.inbams:
        sm_n_inbam_stuff = get_inbam_input_config(args.inbams, args.inbam_dirs, args.include_samples, 
                args.exclude_samples)
        config_out.update(sm_n_inbam_stuff)
    elif (args.analysis_ready_bam_dirs or args.analysis_ready_bams) and \
            not (args.gvcf_dirs or args.gvcfs):
        sm_n_inbam_stuff = get_inbam_input_config(args.analysis_ready_bams,
                args.analysis_ready_bam_dirs, args.include_samples, args.exclude_samples)
        sm_n_inbam_stuff['sample_to_analysis_ready_bam'] = sm_n_inbam_stuff['sample_to_bam']
        del sm_n_inbam_stuff['sample_to_bam']
        del sm_n_inbam_stuff['read_groups']
        config_out.update(sm_n_inbam_stuff)
        sample_to_analysis_ready_bam_index = find_bam_index(
                sm_n_inbam_stuff['sample_to_analysis_ready_bam'])
        config_out.update({'sample_to_analysis_ready_bam_index': sample_to_analysis_ready_bam_index})

    #handle adapter_fasta 
    r1_adapter, r2_adapter = parse_adapter_fasta_file(args.adapter_fasta)
    config_out.update({'r1_adapter': r1_adapter}) #r1_adapter and r2_adapter are None, if no adapter_fasta is provided
    config_out.update({'r2_adapter': r2_adapter})

    #handle passed gvcfs
    if args.gvcf_dirs or args.gvcfs:
        all_gvcf_in_dirs = []
        for dir_ in args.gvcf_dirs:
            all_gvcf_in_dirs.extend([os.path.join(dir_, f) for f in os.listdir(dir_)
                if f.endswith('.g.vcf.gz')])
        passed_gvcfs = args.gvcfs + all_gvcf_in_dirs
        config_out.update({'passed_gvcfs': passed_gvcfs})

        samples_in_vcfs = []
        for avcf in passed_gvcfs:
            samples_in_vcfs.extend(get_sample_list_in_vcf(avcf))
        config_out['sample_names'] = samples_in_vcfs

        # sm_n_inbam_stuff = get_inbam_input_config(args.analysis_ready_bams,
                # args.analysis_ready_bam_dirs, samples_in_vcfs, [])
        # bams_for_readBackedPhasing = sm_n_inbam_stuff['sample_to_bam'].values()
        # config_out.update({'bams_for_readBackedPhasing': bams_for_readBackedPhasing})

    config_out['sample_names'] = sorted(config_out['sample_names'])
    sample_names = config_out['sample_names']

    #handle other stuff from args
    adict = {}
    # max_threads >= threads_on_aggregate_samples >= threads_per_sample
    max_threads, threads_per_sample, threads_on_aggregate_samples = args.max_threads, \
            args.threads_per_sample, args.threads_on_aggregate_samples
    _max_threads_single_job_sacgf = 46 #SACGF has only 48 cpu cores
    adict['t_per_sm'] = threads_per_sample if threads_per_sample < _max_threads_single_job_sacgf \
                            else _max_threads_single_job_sacgf
    adict['t_max'] = threads_on_aggregate_samples if threads_on_aggregate_samples < _max_threads_single_job_sacgf \
                            else _max_threads_single_job_sacgf
    args_dict = vars(args)
    needed_keys = ('base_name data_type workdir callers stop_after '
            'skips email '
            'not_annotate not_gene_annotate '
            'add_chris_gene_list add_VAF_column spreadsheet_format trio force_GATK_hard_filter force_extra_joint_gt '
            'keep_temp_files '
            'just_get_bam_stats '
            'freebayes_call_param freebayes_filter_param gatk_hc_param gatk_gtGVCF_param '
            'gatk_snp_filters gatk_indel_filters').split() #bed files handled separately; keys needed from the command line options for adict
    for k in needed_keys:
        if k in args_dict:
            adict[k] = args_dict[k]

    #handle bed files:
    call_region_bed = args.call_region_bed
    depthofcoverage_bed = args.depthofcoverage_bed
    qualimap_bed = args.qualimap_bed
    amp_inserts_bed = args.amp_inserts_bed
    amp_regions_bed = args.amp_regions_bed
    if amp_inserts_bed:
        if not depthofcoverage_bed:
            depthofcoverage_bed = amp_inserts_bed
        if not call_region_bed:
            call_region_bed = amp_inserts_bed
    if amp_regions_bed and not qualimap_bed:
        qualimap_bed = amp_regions_bed

    adict['call_region_bed'] = call_region_bed
    adict['depthofcoverage_bed'] = depthofcoverage_bed
    adict['qualimap_bed'] = qualimap_bed
    adict['amp_inserts_bed'] = amp_inserts_bed
    adict['amp_regions_bed'] = amp_regions_bed
    config_out.update(adict)

    config_out['random_str'] = random_file_name()

    # #skip GATK_BQSR if data_type is Panel or Amplicon
    # if config_out["data_type"] in [Panel, Amplicon] and \
            # (config_out['sample_to_units'] or config_out['sample_to_bam']) and \
            # config_out['stop_after'] != 'gatk_realn':
        # skips = config_out["skips"]
        # config_out["skips"] = skips + ['gatk_bqsr'] if 'gatk_bqsr' not in skips else skips
        # #choices=['marking_duplicates', 'gatk_realn', 'gatk_bqsr', 'bam_stats', 'qualimap', 'variant_stats']

    #config_out checker
    if not config_out.get('sample_to_units', None) and \
            not config_out.get('sample_to_bam', None) and \
            not config_out.get('sample_to_analysis_ready_bam', None) and \
            not config_out.get('passed_gvcfs', None):
        exit('\nERROR: No input (fastq, bam, or g.vcf) files are found.')

    #write config_out

    # parser.add_argument('--ccb_queue', action='store_true',
            # help=('Use this option to force to submit jobs onto the CCB queue'))
    # parser.add_argument('--non_scratch_dir_to_write_log_files',
            # help=('Write config file config.${base_name}.yaml to this dir, '
                # 'if --ccb_queue is set and the --out_dir is at /scratch'))

    if args.ccb_queue and args.workdir.startswith('/scratch/'): #will write config file to non_scratch_dir_to_write_log_files
        if not args.non_scratch_dir_to_write_log_files:
            exit('--non_scratch_dir_to_write_log_files is needed when --ccb_queue is set and --out_dir is at /scratch/')
        log_dir = args.non_scratch_dir_to_write_log_files
    else: #just creats logs/ and writes config file there
        log_dir = join(args.workdir, 'logs/')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    snake_config = join(log_dir, 'config.%s.yaml' % args.base_name)
    outhandle = open(snake_config, 'w')
    outhandle.writelines(yaml.dump(config_out, default_flow_style=False))
    outhandle.close()

    #run commands
    hostname = gethostname()
    # is_on_tango = True if hostname.startswith('tango') else False
    max_job_number = '300'

    snake_file = os.path.dirname(os.path.abspath(__file__)) + '/workflows/bwa_gatk_freebayes_annotate.snake'
    log_dir = '~/../../..' + log_dir if log_dir.startswith('/scratch/') else log_dir
    snake_cmd = (
        '{snakemake} all -s {snake_file} --configfile {snake_config} --latency-wait 60 '
        '--max-jobs-per-second 5 --restart-times 3 ' #use --restart-times 3 to avoid the sbatch's socket timed out error
        '--resources io={io_jobs} thread={max_threads} '
        '--printshellcmds -j {max_jobs} --jobname "{{params.job_name}}.{{jobid}}" '
        '--keep-going --rerun-incomplete --timestamp --nolock ').format(
                snakemake=default_config.snakemake,
                snake_config=snake_config,
                snake_file=snake_file, io_jobs=args.max_io_heavy_jobs,
                max_threads=max_threads, max_jobs=max_job_number)
    # if is_on_tango:
    if args.ccb_queue:
        #alway use the ccb queue
        cluster_cmd_part = (
            '''--cluster "sbatch --ntasks-per-node {{params.t}} '''
            '''-p ccb '''
            '''-e {log_dir}/{{params.job_name}}.err -o {log_dir}/{{params.job_name}}.out '''
            '''--mem {{params.mem}} --time {{params.walltime}} --workdir {log_dir}" ''').format(log_dir=log_dir)
    else: #alway use the tango queue and the -C "dc:ep" for forcing to use the VSAN scratch (Tango's default scratch)
        cluster_cmd_part = (
            '''--cluster "sbatch --ntasks-per-node {{params.t}} '''
            '''-p tango -C 'dc:ep' '''
            '''-e {log_dir}/{{params.job_name}}.err -o {log_dir}/{{params.job_name}}.out '''
            '''--mem {{params.mem}} --time {{params.walltime}} --workdir {log_dir}" ''').format(log_dir=log_dir)
    # else:#TODO: make if work on SACGF
        # cluster_cmd_part = (
            # '--cluster "qsub -o {log_dir} -j oe '
            # '-l nodes=1:ppn={{params.t}} -l mem={{params.mem}}mb,vmem={{params.vmem}}mb '
            # '-l walltime={{params.walltime}}"').format(log_dir=log_dir)
    cmd = snake_cmd + cluster_cmd_part

    #Submit the cmd to the resource management system:
    send_email = '--mail-user %s --mail-type END,FAIL ' % args.email if args.email else ''
    # if is_on_tango:
    if args.ccb_queue:
        submitting_cmd = ("sbatch -p ccb --job-name {base_name}_console "
                "--mem 5000 --time 240:00:00 --workdir {log_dir} "
                "{send_email} "
                "--wrap '{cmd}' ").format(base_name=args.base_name, log_dir=log_dir,
                        send_email=send_email, cmd=cmd)
    else:#default tango queue
        submitting_cmd = ("sbatch -p tango -C 'dc:ep'  --job-name {base_name}_console "
                "--mem 5000 --time 240:00:00 --workdir {log_dir} "
                "{send_email} "
                "--wrap '{cmd}' ").format(base_name=args.base_name, log_dir=log_dir,
                        send_email=send_email, cmd=cmd)
    # else:#TODO: make if work on SACGF
        # exit('not work on SACGF??')
    if not args.dryrun:
        print
        stdout, stderr, returncode = run_and_handle_returncode(submitting_cmd, print_cmd=False)
        print 'Job management console was successfully submitted with Job ID ', stdout
    else:
        print cmd + ' --dryrun'

