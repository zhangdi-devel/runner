import os, sys
import variant_tools.project as project
import variant_tools.importer as importer
import variant_tools.exporter as exporter
import variant_tools.update as update
import variant_tools.variant as variant
import variant_tools.compare as compare
import variant_tools.association as association
import variant_tools.pipeline as pipeline
import variant_tools.phenotype as phenotype

class Args:
    def __init__(self):
        self.verbosity = 1
    def Manifest(self):
        for k, v in vars(self).items():
            sys.stderr.write('{}: {}\n'.format(k, v))
    
class Init(Args):
    def __init__(self, Project='test', Build='hg19', Force=False):
        Args.__init__(self)
        self.project = Project
        self.build = Build
        self.force = Force
        self.parent = False
        self.children = False
        
    def Run(self):
        self.Manifest()
        project.init(self)

class Import(Args):
    def __init__(self, Input_files='', Format='', Build='hg19', Sample_name=[], Force=False, Jobs=1):
        Args.__init__(self)
        self.input_files = Input_files
        self.build = Build
        self.format = Format
        self.sample_name = Sample_name
        self.force = Force
        self.jobs = Jobs

    def Run(self):
        self.Manifest()
        importer.importVariants(self)

class Update(Args):
    def __init__(self, Table='variant', From_file='', Format='', Sample_name=[], Jobs=1, Set=[], From_stat=[], Samples=[], Genotypes=[]):
        Args.__init__(self)
        self.table = Table
        self.from_file = From_file
        self.format = Format
        self.build = Build
        self.jobs = Jobs
        self.sample_name =Sample_name
        self.set = Set
        self.from_stat = From_stat
        self.samples = Samples
        self.genotypes = Genotypes

    def Run(self):
        if self.from_file == '' and self.from_stat == []:
            sys.stderr('Must choose a way to update: From_file or From_stat\n')
            return
        self.Manifest()
        update.update(self)

class Phenotype(Args):
    def __init__(self, From_file='', Set=[], From_stat=[], Output=[], Jobs=1, Genotypes=[], Samples=[], Header=[], Delimiter='\t', Na='NA', Limit=-1):
        Args.__init__(self)
        self.from_file = From_file
        self.output = Output
        self.jobs = Jobs
        self.set = Set
        self.from_stat = From_stat
        self.samples = Samples
        self.genotypes = Genotypes
        self.header = Header
        self.delimiter = Delimiter
        self.na = Na
        self.limit = Limit

    def Run(self):
        if self.from_file == '' and self.from_stat == []:
            sys.stderr('Must choose a way to update: From_file or From_stat\n')
            return
        self.Manifest()
        phenotype.phenotype(self)

class Use(Args):
    def __init__(self, Source='', As='', Linked_by=[], Anno_type='', Linked_fields=[], Files=[], Rebuild=False, Jobs=1):
        Args.__init__(self)
        self.source = Source
        setattr(self, 'as', As) #hack ...
        self.linked_by = Linked_by
        self.anno_type = Anno_type
        self.linked_fields = Linked_fields
        self.files = Files
        self.rebuild = Rebuild
        self.jobs = Jobs

    def Run(self):
        self.Manifest()
        annotation.use(self)
        
class Compare(Args):
    def __init__(self, Tables=[], Union=[], Intersection=[], Difference=[], Expression=[], Mode='', Samples=[], Count=False, A_diff_B=[], B_diff_A=[], A_and_B=[], A_or_B=[]):
        self.tables = Tables
        self.union = Union
        self.intersection = Intersection
        self.difference = Difference
        self.expression = Expression
        self.mode = Mode
        self.samples = Samples
        self.count = Count
        self.A_diff_B = A_diff_B
        self.B_diff_A = B_diff_A
        self.A_and_B = A_and_B
        self.A_or_B = A_or_B

    def Run(self):
        self.Manifest()
        compare.compare(self)

class Associate(Args):
    def __init__(self, Variants='', Phenotypes=[], Covariates=[], Var_info=[], Geno_info=[], Geno_name='GT', Methods=[], Group_by=[], Samples=[], Genotypes=[], Discard_samples=[], Discard_variants=[], To_db='', Delimiter=None, Force=False, Jobs=1, Unknown_args=[]):
        Args.__init__(self)
        self.variants = Variants
        self.phenotypes = Phenotypes
        self.covariates = Covariates
        self.var_info = Var_info
        self.geno_info = Geno_info
        self.geno_name = Geno_name
        self.methods =Methods
        self.group_by = Group_by
        self.samples = Samples
        self.genotypes = Genotypes
        self.discard_samples = Discard_samples
        self.discard_variants = Discard_variants
        self.to_db = To_db
        self.delimiter = Delimiter
        self.force = Force
        self.jobs = Jobs
        self.unknown_args = Unknown_args

    def Run(self):
        self.Manifest()
        association.associate(self)

class GeneralOutput(Args):
    def __init__(self, Header=[], Delimiter=None, Na='.', Limit=None, Build='hg19', Group_by=[], All=False, Order_by=[]):
        Args.__init__(self)
        self.header = Header
        self.delimiter = Delimiter
        self.na = Na
        self.limit = Limit
        self.build = Build
        self.group_by = Group_by
        self.all = All
        self.order_by = Order_by
        
class Select(GeneralOutput):
    def __init__(self, From_table='', Condition=[], Samples=[], To_table='', Count=False, Output=[], *args, **kwargs):
        GeneralOutput.__init__(self, *args, **kwargs)
        self.from_table = From_table
        self.condition = Condition
        self.samples = Samples
        self.to_table = To_table
        self.count = Count
        self.output = Output

    def Run(self):
        self.Manifest()
        variant.select(self)

class Exclude(Select):
    def Run(self):
        self.Manifest()
        variant.select(self, reverse=True)

class Output(GeneralOutput):
    def __init__(self, Table='', Fields=[], *args, **kwargs):
        GeneralOutput.__init__(self, *args, **kwargs)
        self.table = Table
        self.fields = Fields

    def Run(self):
        self.Manifest()
        variant.output(self)

class Export(Args):
    def __init__(self, Table='', Filename=None, Output=None, Samples=[], Format='', Build='hg19', Header=[], Jobs=1):
        Args.__init__(self)
        self.table = Table
        self.filename =Filename
        self.output = Output
        self.samples = Samples
        self.format = Format
        self.build = Build
        self.header = Header
        self.jobs = Jobs

    def Run(self):
        self.Manifest()
        exporter.export(self)

class Remove(Args):
    def __init__(self, Type='', Items=[]):
        Args.__init__(self)
        self.type = Type
        self.items = Items

    def Run(self):
        self.Manifest()
        project.remove(self)

class Execute(Args):
    def __init__(self, Specfile='', Pipelines=[], Input=[], Output=[], Jobs=1, Delimiter='\t'):
        Args.__init__(self)
        self.specfile = Specfile
        self.Pipelines = pipelines
        self.input = Input
        self.output = Output
        self.jobs = Jobs
        self.delimiter = Delimiter

    def Run(self):
        self.Manifest()
        pipeline.execute(self)

class Admin(Args):
    def __init__(self, Update_resource='', Mirror_repository='', Merge_samples=False, Rename_samples=[], Rename_table=[], Describe_table=[], Validate_build=False, Validate_sex=False, Save_snapshot=[], Extra_files=[], Load_snapshot='', Set_runtime_option=[], Reset_runtime_option=''):
        Args.__init__(self)
        self.update_resouce = Update_resource
        self.mirror_repository = Mirror_repository
        self.merge_samples = Merge_samples
        self.rename_samples = Rename_samples
        self.rename_table = Rename_table
        self.describe_table = Describe_table
        self.validate_table = Validate_table
        self.validate_sex = Validate_sex
        self.save_snapshot = Save_snapshot
        self.extra_files = Extra_files
        self.load_snapshot = Load_snapshot
        self.set_runtime_option = Set_runtime_option
        self.reset_runtime_option = Reset_runtime_option

    def Run(self):
        self.Manifest()
        project.admin(self)
