from sacgf.util.file_utils import bioinformatics_dirname, mk_path
import ConfigParser
import inspect
import os

def get_config_default(cfg, section, option, default=None):
    if cfg.has_option(section, option):
        return cfg.get(section, option)
    else:
        return default


class Config(object):
    def __init__(self, cfg_file_name=None, **kwargs):
        ''' @param cfg_file_name: configuration file
            @param defaults: overwrite_params: dict of dicts - {'section' : {'option' : value}}
            
            Loads:
            1. ~/sacgf.cfg (system wide config)
            2. cfg_file - optional project specific confg file
            3. ../../config/reference.cfg (in source control)
        '''
        defaults = {"output.base.dir" : None}
        cfg = ConfigParser.SafeConfigParser(defaults)
        self.cfg = cfg

        # ~/sacgf.cfg - base directories
        cfg.read(self.get_sacgf_cfg())
        sacgf_mandatory = {
            "data" : "data.root.dir",
            "reference" : "reference.base.dir",
        }
        self.check_mandatory(sacgf_mandatory, "sacgf.cfg in your home directory (system wide config)")

        # project specific config 
        if cfg_file_name:
            cfg.readfp(open(cfg_file_name))

        overwrite_params = kwargs.get("overwrite_params", {})
        for (section, options_dict) in overwrite_params.iteritems():
            for (option, value) in options_dict.iteritems():
                cfg.set(section, option, value)

        self.check_mandatory({"reference" : "organism"}, "a project-specific configuration file or in the defaults parameter")
        self.check_mandatory({"reference" : "annotation_group"}, "UCSC or Ensembl")
        self.check_mandatory({"reference" : "build"}, "a project-specific configuration file or in the defaults parameter")

        self.reference_build = cfg.get("reference", "build")

        # reference
        reference_cfg_filename = os.path.join(bioinformatics_dirname(), "config/reference.cfg")
        self.cfg.readfp(open(reference_cfg_filename))

        self.reference_genes = cfg.get("reference", "genes")
        self.reference_genome = cfg.get("reference", "genome")
        self.reference_trna = cfg.get("reference", "trna")
        self.ncbi_queries_dir = cfg.get("reference", "ncbi.queries.dir")
        self.brca_cohorts_dir = cfg.get("reference", "brca.cohorts.dir")
        self.reference_mature_mirna_fasta = cfg.get("reference", "mature.mir.fasta")
        self.gtf_to_genes_index = cfg.get("reference", "gtf_to_genes_index")
        self.reference_base_dir = cfg.get("reference", "reference.base.dir")
        self.ensembl_cpickle_data_dir = cfg.get("reference", "ensembl.cpickle.data.dir")
        self.conservation_numpy_array_dir = cfg.get("reference", "conservation.numpy.array.dir")
        self.mappability_bigwig_dir = cfg.get("reference", "mappability.bigwig.dir")

        output_base_dir = cfg.get("data", "output.base.dir")
        if output_base_dir:
            output_base_dir = os.path.expanduser(output_base_dir)
        self.output_base_dir = output_base_dir

        self.gene_summary_tsv = os.path.join(self.ncbi_queries_dir, "gene_summary", "genes_w_GeneSummary.tsv")
        self.publication_counts_dir = os.path.join(self.ncbi_queries_dir, "publication_counts")

    def check_mandatory(self, mandatory, expected_location):
        for (section, option) in mandatory.iteritems():
            if not self.cfg.has_option(section, option):
                raise ValueError("Missing required value for %s.%s. Generally you should set this in %s." % (section, option, expected_location))

    # In case we want to overwrite later
    def get_sacgf_cfg(self):
        return os.path.expanduser("~/sacgf.cfg")

    def get_output_dir(self, **kwargs):
        '''
            Return "project specific" directory to write files.
            
            If sacgf.cfg has an entry for "data / output_base_dir"
            
            AND the package of calling function is "sacgf.project_specific.XXX" then it returns:

            ${output_base_dir}/XXX
            
            Otherwise, return current working directory
        '''
        depth = kwargs.get("depth", 1)
        frm = inspect.stack()[depth]
        modWhole = inspect.getmodule(frm[0])
        f = modWhole.__file__
        parts = f.split(os.path.sep)
        output_dir = None
        if self.output_base_dir:
            try:
                pi = parts.index("project_specific")
                name = parts[pi + 1]
                output_dir = os.path.join(self.output_base_dir, name)
            except ValueError:
                pass

        return output_dir or os.getcwd()

    def out(self, filename):
        ''' Shortcut to return a filename for get_output_dir() / filename
            Creates directory if path doesn't exist
        '''
        full_path = os.path.join(self.get_output_dir(depth=2), filename)
        mk_path(filename)
        return full_path

class DefaultConfig(Config):
    def __init__(self, **kwargs):
        '''
        Key word arguments:
            organism (default=Homo_sapiens)
            annotation_group (default=UCSC)
            build (default=hg19)
        '''
        reference_params = {"organism" : kwargs.get("organism", "Homo_sapiens"),
                            "annotation_group" : kwargs.get("annotation_group", "UCSC"),
                            "build" : kwargs.get("build", "hg19")}
        super(DefaultConfig, self).__init__(None, overwrite_params={"reference" : reference_params})

class EnsemblConfig(Config):
    def __init__(self, **kwargs):
        '''
        Key word arguments:
            organism (default=Homo_sapiens)
            annotation_group (default=Ensembl)
            build (default=GRCh37)
        '''
        reference_params = {"organism" : kwargs.get("organism", "Homo_sapiens"),
                            "annotation_group" : "Ensembl",
                            "build" : kwargs.get("build", "GRCh37")}
        super(EnsemblConfig, self).__init__(None, overwrite_params={"reference" : reference_params})

