"""
Some functions here rely on $SACGF_DIR
"""
import os
import gzip
from itertools import izip
import vcf
from cyvcf2 import VCF
import pandas as pd
from sacgf.util.file_utils import remove_a_file, line_count_wc, random_file_name, \
        absolute_file_path, random_file_at_same_dir
from sacgf.util.subprocess_utils import run_subprocess_cmd, run_and_handle_error, run_and_handle_returncode
from scripts.pipelines.variant_annotation_and_filtration.variant_and_gene_annotation_pipeline_parameters import BGZIP_CMD, \
    TABIX_CMD, GATK_JAR, B37_REF, BCFTOOLS, VAR_ANNOT_FILE_DIR, IGVTOOLS_CMD


# TODO KEY_PATTERN is not needed anymore
# TODO LOAD_JAVA7_COMMAND won't be needed in the near future
KEY_PATTERN = '%s:%s:%s'  # for 'CHROM:POS:REF'
LOAD_JAVA7_COMMAND = 'module load java/java-1.7.09' # the command in SACGF server
ALL_CHRS = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']

def append_new_info(orig_line, formated_info_to_add):
    """
    append the formated_info_to_add to the INFO section in the orig_line
    and return the new plain text line with the formated_info_to_add
    """
    orig_line_elments = orig_line.rstrip().split('\t')
    try:
        orig_info = orig_line_elments[7]
        if orig_info == '.':
            orig_line_elments[7] = formated_info_to_add
        else:
            orig_line_elments[7] = orig_info + ';' + formated_info_to_add
    except IndexError:
        import sys; sys.exit('Python IndexError (malformed VCF?) in infile at line:\n%s\n', orig_line)
    return '\t'.join(orig_line_elments) + '\n'

def fix_vcf_header_for_virmid(virmid_vcf):
    """
    Remove the illegal space in the header line in the VCF file produced by VIRMID.
    The original virmid_vcf will be changed in place.
    """
    virmid_vcf = os.path.expanduser(virmid_vcf)
    cmd = "sed -i 's/ID=gt,Number=A,Type=String, Description=/ID=gt,Number=A,Type=String,Description=/' %s" % virmid_vcf
    run_subprocess_cmd(cmd)

def get_outhandle_with_new_headers(invcf, outvcf, *new_headers):
    """
    Based on the invcf, appending new header lines and putting them on top of the last header line (a.k.a. #CHROM POS...)
    Returns the file handle of outvcf (not gzipped format), which can continue writing or just be closed.
    '' for new_headers means returning the out handle with the original header lines of invcf
    """
    new_header_lines = '\n'.join([line.rstrip() for line in new_headers])
    infilehandle = gzip.open(invcf) if invcf.endswith('.vcf.gz') else open(invcf)
    outhandle = open(os.path.expanduser(outvcf), 'w')
    header_number = count_header_lines(invcf)
    for i in range(header_number):
        if i == header_number - 1: # the one before the last header
            if new_header_lines:
                outhandle.write(new_header_lines + '\n')
            else:
                pass
        outhandle.write(infilehandle.next())
    return outhandle

    
def save_to_vcf(template_filepath, vcf_outpath, records):
    """
    Save the records of pyvcf to a vcf file with meta-data and header 
    from the template.
    The template_filepath is where the records are from, normally. 
    records is a list of record objects of pyvcf. 
    """

    vcf_outpath = os.path.expanduser(vcf_outpath)
    template_filepath = os.path.expanduser(template_filepath)

    outfile_handle = vcf.Writer(open(vcf_outpath, 'w'), 
            template=vcf.Reader(filename=template_filepath), lineterminator='\n')

    # pyvcf 0.6.3 only allows to write record one by one
    for record in records:
        outfile_handle.write_record(record)

    outfile_handle.close()

def get_vcf_writer_handle(template_filepath, vcf_outpath):
    """ return a pyvcf writer handle for later use """

    # create outfile_handle with template for saving a VCF file
    template_handle = open(os.path.expanduser(template_filepath), 'rU')
    outfile_handle = vcf.Writer(open(os.path.expanduser(vcf_outpath), 'w'), 
            template=vcf.Reader(template_handle), lineterminator='\n')
    template_handle.close()
    return outfile_handle

def get_vcf_iterator_in_plain_text(fpath_vcf):
    """ returns the iterator/in-handle for a vcf file in the plain text mode, with skipping the header lines"""
    fpath_vcf = os.path.expanduser(fpath_vcf)
    inhandle = gzip.open(fpath_vcf) if fpath_vcf.endswith('.vcf.gz') else open(fpath_vcf)
    for i in range(count_header_lines(fpath_vcf)):
        inhandle.next()
    return inhandle

def get_sample_list(fpath_vcf):
    """returns a list of the samples in fpath_vcf, in the same order as listed in the header"""
    return vcf.Reader(filename=os.path.expanduser(fpath_vcf)).samples

def get_the_first_variant(fpath_vcf):
    """
    return the first variant in pyvcf object
    """
    return vcf.Reader(filename=fpath_vcf).next()

def get_info_tags(fpath_vcf):
    """returns a list of all info tags in the header lines"""
    return vcf.Reader(filename=fpath_vcf).infos.keys()


def count_multiple_ALT(fpath_vcf):
    """count and return how many records have multiple ALT alleles"""
    count_multi_ALT = 0
    count_total = 0
    with open(os.path.expanduser(fpath_vcf)) as infile:
        reader = vcf.Reader(infile)
        for record in reader:
            count_total += 1
            if len(record.ALT) > 1:
                count_multi_ALT += 1
    proportion_of_mult_ALT = float(count_multi_ALT) / count_total
    return count_multi_ALT, proportion_of_mult_ALT


def load_vcf_into_dict(fpath_vcf):
    """
    Load all variations of case into a dict,
    returns {'CHROM:POS:REF': [(ALT, record), (ALT_2, record_2)], }
    Be careful if the VCF file is too huge to load into memory.
    """
    all_variants = {}
    with open(os.path.expanduser(fpath_vcf)) as source:
        for record in vcf.Reader(source):
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT
            key = KEY_PATTERN % (chrom, pos, ref)
            try:
                all_variants[key].append((alt, record))
            except KeyError:
                all_variants[key] = [(alt, record)]
    return all_variants

def count_header_lines(fpath_vcf):
    """
    Counts and returns the number of lines starting with '#' (i.e. the meta-information 
    and header lines). The first line in the VCF file must start with '#', or zero will 
    always be returned.
    """
    count = 0
    inhandle = gzip.open(fpath_vcf) if fpath_vcf.endswith('.vcf.gz') else open(fpath_vcf)
    for line in inhandle:
        if line.startswith('#'):
            count += 1
        else:
            break
    return count

def format_existing(fpath_vcf):
    """returns True if the FORMAT section is existing, or returns False"""
    last_header = extract_header_lines(fpath_vcf)[-1]
    if len(last_header.rstrip().split('\t')) >= 9:
        return True
    else:
        return False

def having_contig_in_header(invcf):
    if vcf.Reader(filename=invcf).contigs:
        return True
    else:
        return False

def are_contigs_in_karyotypic_order(invcf):
    """
    Check if the contigs are listed in the karyotypic order in the header section.
    The contigs listed in the header are assumed to be listed in the same order as in the body
    Note it returns True if there is no contigs defined in the header section.
    """
    contig_to_value = {}
    chr_values = [i for i in range(1,26)]
    for chr_, value in zip(ALL_CHRS, chr_values):
        contig_to_value[chr_] = value

    contigs = vcf.Reader(filename=invcf).contigs.keys() 
    values_in_invcf = []
    for contig in contigs:
        try:
            value = contig_to_value[contig]
        except KeyError:
            value = 9999
        values_in_invcf.append(value)
    return sorted(values_in_invcf) == values_in_invcf

def basic_filtration_in_a_row(arow, impact_snpeff=['HIGH', 'MODERATE', 'LOW'],
        gerp_cutoff=2, cadd_cutoff=10, maf_cutoff=0.01):
    """
    (Nov/2014 version)
    The Basic filtration described in the WiKi page.
    Filter for a rare and functional variant in a pandas series with correct dtype.
    Used this after loading a tsv/excel file with annotated variants (Nov/2014 version).
    Return 'PASS' if it satisfies the criteria, else '' is returned
    """

    if pd.isnull(arow['1KG_ALT_Freq']) or min(arow['1KG_ALT_Freq'], (1-arow['1KG_ALT_Freq'])) < maf_cutoff:
        if (arow.IMPACT_SNPEFF in impact_snpeff and arow.EFFECT_SNPEFF != 'synonymous_variant') or \
                (arow.TYPE == 'SNP' and (arow.GERP >=2 or arow.CADD >= 10)):
            return 'PASS'
    return ''

def stringent_filter(arow, impact_snpeff=['HIGH', 'MODERATE'], maf_cutoff=0.005):
    """
    Return 'PASS' if:
        "PASS, HIGH or MODERATE, and < 0.5% in both 1KG and ESP"
    else returns ''
    """

    if arow.FILTER == 'PASS':
        esp_maf = arow['ESP_MAFs(%)']
        if pd.isnull(esp_maf):
            esp_maf = None
        else:
            esp_maf = float(esp_maf.split(',')[-1]) * 100

        if (pd.isnull(arow['1KG_ALT_Freq']) or \
                min(arow['1KG_ALT_Freq'], (1-arow['1KG_ALT_Freq'])) < maf_cutoff) and (not esp_maf or \
                esp_maf < maf_cutoff):
            if (arow.IMPACT_SNPEFF in impact_snpeff):
                return 'PASS'
    return ''

def extract_header_lines(fpath_vcf, with_one_record_line=False):
    """
    returns a list in which each header is one element. 
    header line starts with '#'
    at least one record lines are needed to serve as a template for the vcf.Writer()
    """
    fpath_vcf = os.path.expanduser(fpath_vcf)
    source = gzip.open(fpath_vcf) if fpath_vcf.endswith('.vcf.gz') else open(fpath_vcf)
    #with open(os.path.expanduser(fpath_vcf)) as source:
    returned_list = []
    for line in source:
        if line.startswith('#'):
            returned_list.append(line)
        else:
            if with_one_record_line:
                returned_list.append(line)
                break
            else:
                break
    return returned_list

def count_variant_number(fpath_vcf):
    """count the number of variant in a VCF file. Works with both .vcf and .vcf.gz"""
    return line_count_wc(fpath_vcf) - count_header_lines(fpath_vcf)

def all_are_hom_major_allele(record, af_1kg, threshold=0.5):
    """
    Return True if all successfully called samples are homozygous for the major (>0.5, default) allele.
    record is the variant record by pyvcf; 
    af_1kg is from the 1000 Genomes Project, and is None when not available.
    """
    # No af_1kg available
    if af_1kg is None:
        return False
    # returns False when one or more of the called samples are heterozygous
    elif record.num_het > 0:
        return False
    # all called samples are homozygous for ALT, and the ALT frequency > threshold
    elif record.num_hom_ref == 0 and af_1kg > threshold:
        return True
    # all called samples are homozygous for REF, and the REF frequency > threshold 
    # we usually don't keep non-variants, so this case won't happen
    elif record.num_hom_alt == 0 and (1- af_1kg) > threshold:
        return True
    else: # return False when requried information is missing
        return False

def is_gzip_or_bcf(invcf):
    if invcf.endswith('.vcf.gz') or invcf.endswith('.bcf'):
        return True
    else:
        return False

def remove_variants_without_sample_genotypes(invcf, outvcf):
    """remove the variants where all samples are "./" for genotype"""
    
    outvcf_handle = get_vcf_writer_handle(invcf, outvcf)
    with open(os.path.expanduser(invcf)) as source:
        reader = vcf.Reader(source)
        for record in reader:
            if record.call_rate == 0:
                continue
            else:
                outvcf_handle.write_record(record)
    outvcf_handle.close()

def remove_info_fields(invcf, outvcf, info_fields):
    """
    Deprecated in favour to bcftools.
    Infor fields in info_fields will be removed from a VCF file (i.e. removes annotations).
    The same feature from SnpSift does not remove the lines in the meta-information section, plus 
    SnpSift only adds and does not add/modify lines in a VCF when doing annotation, so sometimes 
    cause some unexpected issues.
    Note there must be a dot when no info field according to the VCF spec. -- 
    "The VCF specification requires a valid (non-zero length) info field", according to igvtools
    """
    meta_line_pattern = '##INFO=<ID=%s,'
    info_element_pattern = '%s=' # note there is Type=Flag as well
    num_meta_info = count_header_lines(invcf)
    num_total_file_line = line_count_wc(invcf)
    if num_meta_info == num_total_file_line:
        print 'No record lines are in the VCF file %s, so no OUTVCF is generaged.' % invcf
        return False
    outhandle = open(os.path.expanduser(outvcf), 'w')
    with open(os.path.expanduser(invcf)) as source:
        # handling the meta-info section
        for i in range(num_meta_info):
            meta_line = source.next()
            hit_meta = False
            for info_field in info_fields:
                if meta_line.startswith(meta_line_pattern % info_field):
                    hit_meta = True
                    break
            if hit_meta:
                continue
            else:
                outhandle.write(meta_line)
        # now handling the record lines
        new_lines = ''
        new_line_num = 0
        for line in source:
            # rstrip() is needed, as INFO can be the last column
            line_elements = line.rstrip().split('\t')
            info = line_elements[7]
            info_elements = info.split(';')
            new_info_elements = []
            for element in info_elements:
                hit_info = False
                for info_field in info_fields:
                    if element.startswith(info_element_pattern % info_field) or element == info_field:
                        hit_info = True
                        break
                if hit_info:
                    continue
                else:
                    new_info_elements.append(element)
            new_info = ';'.join(new_info_elements)
            if not new_info:
                new_info = '.' # add '.' if there is no info element left.
            # new line_elements
            line_elements[7] = new_info
            new_line = '\t'.join(line_elements) + '\n'
            new_lines = ''.join([new_lines, new_line])
            #new_lines += new_line
            new_line_num += 1
            if new_line_num == 10000:
                outhandle.write(new_lines)
                new_lines = ''
                new_line_num = 0
        if new_line:
            outhandle.write(new_lines)
        outhandle.close()
    return True

def remove_vcf_file(*vcf_file_path):
    """remove .vcf, .vcf.gz or bcf files as well as their index or tab-index files"""
    for fpath in vcf_file_path:
        # remove the VCF file
        remove_a_file(os.path.expanduser(fpath))
        # remove the corresponding index file
        remove_a_file(os.path.expanduser(fpath) + '.idx')
        remove_a_file(os.path.expanduser(fpath) + '.tbi')
        remove_a_file(os.path.expanduser(fpath) + '.csi')



def create_gatk_vcf_list_statements(vcf_file_list):
    """
    return for GATK usage something like "--variant /sacgf/variant1.vcf --variant /temp/adfa.vcf"
    """
    statements = ''
    for vcf in vcf_file_list:
        statements += '--variant %s ' % vcf
    return statements

def combine_vcf_by_gatk(vcf_file_list, output_vcf, loadJava7=False):
    """
    Only for hg19 now.
    Run GATK to combine VCF files in vcf_file_list. The samples are REQURIED to be different in all
    files in the vcf_file_list
    """

    if loadJava7:
        load_java7 = LOAD_JAVA7_COMMAND
    else:
        load_java7 = ''

    gatk_command = '''
    %(load_java7)s
    java -Xmx4g -jar %(GATKCOMMAND)s \
        -R %(REF)s \
        -T CombineVariants \
        %(vcf_file_statement)s \
        -o %(output_vcf)s 
            ''' % {
                    'GATKCOMMAND': GATK_JAR,
                    'REF': B37_REF,
                    'output_vcf': output_vcf,
                    'load_java7': load_java7,
                    'vcf_file_statement': create_gatk_vcf_list_statements(vcf_file_list)
                    }
    run_subprocess_cmd(gatk_command)
    return 'Done'

def select_variants_by_gatk(invcf, outvcf, is_freebayes_somatic_vcf=False, sample_list=None,
        intervals=None, interval_padding=False, exclude_sample=False, hugeVCF=False, num_threads=None):
    """
    Run GATK -T SelectVariants, specifing the sample_list and/or intervals.
    When is_freebayes_somatic_vcf is False (the default), it uses --removeUnusedAlternates and --excludeNonVariants.
    intervals is a file path (e.g. A BED file), a single interval string ('1:100-200'),
    a single point position ('14:45645954') or a list of string like ['1:100-200', '14:45645954']
    outvcf must be a vcf.gz file.
    Only for hg19 now.
    """
    #TODO remove the tmp_vcf_out when GATK fix the * non-variant issue
    tmp_vcf_out = random_file_at_same_dir(outvcf, prefix='before_removing_star', extension='vcf.gz')

    sample_statement = '' 
    if sample_list:
        if exclude_sample:
            for sample in sample_list:
                sample_statement += '--exclude_sample_name %s ' % sample
        else:
            for sample in sample_list:
                sample_statement += '--sample_name %s ' % sample

    interval_statements = ''
    if intervals:
        if type(intervals) is list:
            for interval in intervals:
                interval_statements += ' --intervals %s ' % interval
        else: # it is a file path or a single interval string
            interval_statements += ' --intervals %s ' % intervals

    if interval_padding:
        interval_statements += ' --interval_padding %s ' % interval_padding

    if hugeVCF:
        java_Xmx = '-Xmx14g'
    else:
        java_Xmx = '-Xmx4g'

    if num_threads:
        nt_statement = '--num_threads %s' % num_threads
    else:
        nt_statement = ''

    if is_freebayes_somatic_vcf:
        remove_unused_alternates_option = ''
    else:
        remove_unused_alternates_option = '--removeUnusedAlternates --excludeNonVariants'

    # GATK SelectVariants does normalizing selected variants
    gatk_command = '''
    java %(java_Xmx)s -jar %(GATKCOMMAND)s \
        -R %(REF)s \
        --variant %(invcf)s \
        -o %(output_vcf)s \
        %(sample_statement)s  %(interval_statements)s \
        %(nt_statement)s \
        -T SelectVariants \
        %(remove_unused_alternates_option)s
    ''' % {
            'GATKCOMMAND': GATK_JAR,
            'java_Xmx': java_Xmx,
            'REF': B37_REF,
            'invcf': invcf,
            'output_vcf': tmp_vcf_out,
            'sample_statement': sample_statement,
            'nt_statement': nt_statement,
            'interval_statements': interval_statements,
            'remove_unused_alternates_option': remove_unused_alternates_option
            }
    run_subprocess_cmd(gatk_command)

    #TODO remove the tmp_vcf_out when GATK fix the * non-variant issue
    remove_star_cmd = '''
    {bcftools} view -O z --exclude "ALT='*'" {invcf} > {outvcf}
    {bcftools} index -t {outvcf}
    '''.format(bcftools=BCFTOOLS, invcf=tmp_vcf_out, outvcf=outvcf)
    run_subprocess_cmd(remove_star_cmd)
    remove_vcf_file(tmp_vcf_out)

    return 'Done'

def select_samples_by_bcftools(samples, invcf, outvcf, exclude=False):
    """
    select samples followed by lef-alignment at indels
    invcf and outvcf must endswith .vcf.gz
    """
    #FIXME variants called on non-decoy ref probally have issue on Y PAR regions.
    invcf = absolute_file_path(invcf)
    outvcf = absolute_file_path(outvcf)
    if not (is_gzip_or_bcf(invcf) and is_gzip_or_bcf(outvcf)):
        import sys; sys.exit('Both invcf and outvcf must be .vcf.gz')

    sample_str = ','.join(samples)
    if exclude:
        sample_str = '^' + sample_str

    cmd ="""
    {bcftools} view --exclude-uncalled --min-ac 1 --trim-alt-alleles \\
        --samples {sample_str} \\
        {invcf} | {bcftools} norm -f {ref} -O z -o {outvcf}
    """.format(bcftools=BCFTOOLS, sample_str=sample_str, ref=B37_REF,
            invcf=invcf, outvcf=outvcf)
    run_and_handle_returncode(cmd)
    tabix_on_bgzipped_vcf(outvcf)


def create_vcf_index(vcf_file):
    """
    Uses igvtools to create an index file for vcf_file. Relies on $SACGF environmental viable.
    igv.log will be removed """
    cmds = '''
    {IGVTOOLS_CMD} index {VCF_FILE}
    rm -f igv.log
    '''.format(IGVTOOLS_CMD=IGVTOOLS_CMD, VCF_FILE=vcf_file)
    run_subprocess_cmd(cmds)

def create_bgzip_tabix_vcf(vcf_file, outfile=None):
    """
    run bgzip and tabix to the vcf_file.
    By defualt, the vcf_file will be replaced by compressed file, with '.gz' appended to the file name.
    """
    if outfile and not outfile.endswith('.vcf.gz'):
        import sys; sys.exit('outfile does not end with ".vcf.gz": %s' % outfile)

    vcf_file = absolute_file_path(vcf_file)
    outfile = absolute_file_path(outfile)

    if outfile:
        bgzip_and_tabix_cmd = '''
    {BGZIP_CMD} -f -c {vcf_file} > {outfile}
    {TABIX_CMD} -p vcf -f {outfile}
    '''.format(BGZIP_CMD=BGZIP_CMD, TABIX_CMD=TABIX_CMD, vcf_file=vcf_file, outfile=outfile)
    else:
        bgzip_and_tabix_cmd = '''
    {BGZIP_CMD} -f {vcf_file}
    {TABIX_CMD} -p vcf -f {vcf_file}.gz
    '''.format(BGZIP_CMD=BGZIP_CMD, TABIX_CMD=TABIX_CMD, vcf_file=vcf_file)
    run_and_handle_error(bgzip_and_tabix_cmd)

def tabix_on_bgzipped_vcf(vcf_file):
    """
    run tabix on a vcf.gz file
    """
    cmd = '{TABIX_CMD} -p vcf -f {vcf_file}'.format(TABIX_CMD=TABIX_CMD, vcf_file=vcf_file)
    run_and_handle_error(cmd)

def decompress_bgzipped_file(infile, outfile):
    """decompress bgzipped file: from .vcf.gz to .vcf"""
    infile = absolute_file_path(infile)
    outfile = absolute_file_path(outfile)
    cmd = '{BGZIP_CMD} -d -c -f {infile} > {outfile}'.format(
                BGZIP_CMD=BGZIP_CMD, infile=infile, outfile=outfile)
    run_and_handle_error(cmd)

def sort_vcf_by_igvtools(invcf, outvcf):
    """Uses igvtools to sort a VCF file"""
    cmds = '''
    %s sort %s %s
    rm -f igv.log
    ''' % (IGVTOOLS_CMD, invcf, outvcf)
    run_and_handle_returncode(cmds)

def split_and_left_align(invcf, outvcf, ref=B37_REF, ignored_bed=None):
    """
    Used for invcf without GT section
    use 'bcftools norm' to split split multiallelic sites into multiple VCF record rows, and to
    left-align indels.
    both invcf and outvcf must be .vcf.gz (with .tbi).
    Sometimes PARs in Y need to be skipped if B37_REF is used.
    """
    if not (is_gzip_or_bcf(invcf) and is_gzip_or_bcf(outvcf)):
        import sys; sys.exit('Both invcf and outvcf must be .vcf.gz file')

    invcf = absolute_file_path(invcf)
    outvcf = absolute_file_path(outvcf)

    if ignored_bed: #bed file defining the regions to skip
        ignored_str = '-T ^%s' % ignored_bed
    else:
        ignored_str = ''

    tmp_vcf = random_file_at_same_dir(outvcf, extension='vcf')
    cmd = '''
    {bcftools} norm -m -any {invcf} | {bcftools} norm {ignored_str} -f {ref}  -o {tmp_vcf} -
    '''.format(bcftools=BCFTOOLS, ref=ref, ignored_str=ignored_str, tmp_vcf=tmp_vcf, invcf=invcf)
    run_and_handle_returncode(cmd)

    #sort tmp_vcf
    tmp_sorted_vcf = random_file_at_same_dir(outvcf, prefix='sorted', extension='vcf')
    sort_vcf_by_igvtools(tmp_vcf, tmp_sorted_vcf)

    create_bgzip_tabix_vcf(tmp_sorted_vcf, outvcf)

    # the last step: remove temp files
    remove_vcf_file(tmp_vcf, tmp_sorted_vcf)

def add_control_info(case_vcf, ctrl_vcf, case_w_ctrl_info_vcf):
    """
    Using bcftools to add CTRL_AC, CTRL_AN and CTRL_AF to case_vcf to create case_w_ctrl_info_vcf
    """
    #TODO Not sure if it works with .vcf file for case_vcf and ctrl_vcf. 
    #It is tested only with .vcf.gz and it works.
    if not case_w_ctrl_info_vcf.endswith('.vcf.gz'):
        import sys; sys.exit('outfile does not end with ".vcf.gz": %s' % case_w_ctrl_info_vcf)

    header_file = os.path.join(VAR_ANNOT_FILE_DIR, 'control_header.txt')

    #create a tsv file with necessary info from the ctrl_vcf
    #first split multiple ALT alleles, then left-align
    tmp_tsv = random_file_at_same_dir(case_w_ctrl_info_vcf, prefix='ctrl_info_only', extension='tsv')
    query_f_str = r'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n'
    awk_par_str = r'{print $0 "\t" $5/$6}'
    cmd = '''
    {bcftools} norm --multiallelics -any {ctrl_vcf} | \\
            {bcftools} norm -f {ref} - | \\
            {bcftools} view - --exclude-uncalled | \\
            {bcftools} query -f '{query_f_str}' - | \\
            awk '{awk_par_str}' > {tmp_tsv}
    '''.format(bcftools=BCFTOOLS, ctrl_vcf=ctrl_vcf, ref=B37_REF, 
            query_f_str=query_f_str, awk_par_str=awk_par_str, tmp_tsv=tmp_tsv)

    run_and_handle_returncode(cmd)
    create_bgzip_tabix_vcf(tmp_tsv)

    #run bcftools annotate
    tmp_tsv_zip=tmp_tsv+'.gz'
    cmd = '''
    {bcftools} annotate -a {tmp_tsv_zip} --columns CHROM,POS,REF,ALT,CTRL_AC,CTRL_AN,CTRL_AF \\
            --header-lines {header_file} -O z {case_vcf} > {case_w_ctrl_info_vcf}
    '''.format(bcftools=BCFTOOLS, tmp_tsv_zip=tmp_tsv_zip, case_vcf=case_vcf, 
            header_file=header_file, case_w_ctrl_info_vcf=case_w_ctrl_info_vcf)
    run_and_handle_returncode(cmd)

    tabix_on_bgzipped_vcf(case_w_ctrl_info_vcf)

    remove_a_file(tmp_tsv_zip)
    remove_a_file(tmp_tsv_zip+'.tbi')

def standard_filtration_with_ccc(case_samples,vcf_file_list,outvcf,standard_filter='basic',
        tolerated_allele_num=0, ignored_samples=[], max_cpu_num=1, keep_interm_files=False, loadJava7=False):
    """
    The case samples are defined by case_samples; The control samples are non-case with the exclusion of ignored_samples.

    First, combine all vcf_file_list if more than one VCF files are there; 
    Second, perform standard_filter if specified;
    Third, extract case and control samples and then performs ccc. 
    """

    out_dir = os.path.dirname(outvcf)

    # combine all vcf_file
    if len(vcf_file_list) > 1:
        combined_vcf = os.path.join(out_dir, 'combined.%s.vcf' % random_file_name())
        combine_vcf_by_gatk(vcf_file_list,combined_vcf,loadJava7=loadJava7)
        combined_vcf_is_tmp = True
    else:
        combined_vcf = vcf_file_list[0]
        combined_vcf_is_tmp = False

    # perform standard_filter
    if standard_filter:
        filtered_vcf = os.path.join(out_dir, 'filtered.%s.vcf' % random_file_name())
        import standard_filters
        standard_filters.main(combined_vcf,filtered_vcf,var_filter=standard_filter)
        filtered_vcf_is_tmp = True
    else:
        filtered_vcf = combined_vcf
        filtered_vcf_is_tmp = False

    # extract case and control samples
    case_vcf = os.path.join(out_dir, 'case.%s.vcf' % random_file_name())
    control_vcf = os.path.join(out_dir, 'control.%s.vcf' % random_file_name())
    select_variants_by_gatk(filtered_vcf,case_vcf,sample_list=case_samples,loadJava7=loadJava7)
    select_variants_by_gatk(filtered_vcf,control_vcf, sample_list=case_samples + ignored_samples, 
            exclude_sample=True,loadJava7=loadJava7)

    # perform ccc
    import filter_case_against_control_vcf
    filter_case_against_control_vcf.main(case_vcf,control_vcf,outvcf,
            tolerated_allele_num=tolerated_allele_num, max_cpu_num=max_cpu_num)

    # remove intermediate files
    if not keep_interm_files:
        if combined_vcf_is_tmp:
            remove_vcf_file(combined_vcf)
        if filtered_vcf_is_tmp:
            remove_vcf_file(filtered_vcf)
        remove_vcf_file(case_vcf, control_vcf)

    return True

def validate_vcf_file(infile):
    """
    validate the .vcf or .vcf.gz file, for annotation pipeline
    Uses this before adding pipeline annotation.
    Returns '' if passing all checkers, or a string of ERROR MESSAGE will be returned
    Only the first record line will be checked.
    """
    infile = absolute_file_path(infile)
    if not os.path.isfile(infile):
        return '%s is not a file' % infile
    else:
        inhandle = gzip.open(infile) if infile.endswith('.vcf.gz') else open(infile)
        if not inhandle.next().startswith('##fileformat=VCFv4'):
            return 'The first line does not start with "##fileformat=VCFv4" in %s' % infile
        else:
            last_header = None
            for line in inhandle:
                if line.startswith('#CHROM'):
                    last_header = line
                    break
            if not last_header:
                return 'There is no #CHROM header line in %s' % infile
            else:
                try:
                    first_record_line = inhandle.next()
                except StopIteration:
                    return 'There is no record line in %s' % infile #TODO Not to change this message 
                if len(first_record_line.rstrip()) == 0:
                    return 'There is no record line in %s' % infile
                else:
                    last_header_columns = last_header.rstrip().split('\t')
                    first_record_columns = first_record_line.rstrip().split('\t')
                    if len(last_header_columns) < 8: #8 fixed, mandatory columns
                        return ('The number of columns in the last header line is less than 8, '
                                    'or it is not tab-delimited, in %s' % infile)
                    elif len(first_record_columns ) < 8:
                        return ('The number of columns in the first record line is less than 8, '
                                    'or it is not tab-delimited, in %s' % infile)
                    elif  len(last_header_columns) != len(first_record_columns ):
                        return ('The numbers of columns are not the same between the last header '
                                    'line and the first record line in %s' % infile)
                    elif len(last_header_columns) > 8:
                        if last_header_columns[8] != 'FORMAT':
                            return '"FORMAT" is not the 9th column in the header line in %s' % infile
                        elif not first_record_columns[8].startswith('GT'):
                            return '"GT" is not the first tag under the FORMAT column in %s' % infile
    return ''

def validate_vcf_spec_by_gatk(infile):
    """To perform only VCF format tests, by GATK"""
    gatk_command = '''
    java -Xmx4g -jar %(GATKCOMMAND)s \
        -R %(REF)s \
        -T ValidateVariants \
        --validationTypeToExclude ALL \
        --variant %(infile)s 
            ''' % {
                    'GATKCOMMAND': GATK_JAR,
                    'REF': B37_REF,
                    'infile': infile
                    }
    run_subprocess_cmd(gatk_command) # GATK writes stdout to stderr

def having_same_info_items(vcf_1, vcf_2, exclude_tags = ['HGNC_SNPEFF']):
    """
    Check if two vcf files sharing the same tag-value pairs in the INFO section.
    Assumed the tow vcf files have the same variant at each record line.
    """
    info_tags_vcf1 = get_info_tags(vcf_1)
    info_tags_vcf2 = get_info_tags(vcf_2)
    if sorted(info_tags_vcf1) != sorted(info_tags_vcf2):
        print 'The two VCF files do NOT share the same info tags in the header section'
        return False
    tags_to_compare = [t for t in info_tags_vcf1 if t not in exclude_tags]
    for v1, v2 in izip(VCF(vcf_1), VCF(vcf_2)):
        chr_1, pos_1 = v1.CHROM, v1.POS
        chr_2, pos_2 = v2.CHROM, v2.POS
        if chr_1 != chr_2 or pos_1 != pos_2:
            print 'The two VCF do NOT have the same chr or pos:'
            print 'vcf_1: chr %s, pos %s' % chr_1, pos_1
            print 'vcf_2: chr %s, pos %s' % chr_2, pos_2
            return False
        # if chr_1 == '1' and pos_1 == 150969992:
            # import ipdb; ipdb.set_trace()
        for tag in tags_to_compare:
            if v1.INFO.get(tag) != v2.INFO.get(tag):
                print 'For tag %s, the two VCF files do not share the same value at CHR %s POS %s' % (tag, chr_1, pos_1)
                print 'vcf_1 value %s, vcf_2 value %s' % (str(v1.INFO.get(tag)), str(v1.INFO.get(tag)))
                return False
    return True
