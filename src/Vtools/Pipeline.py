#!/usr/bin/python
import sys, os, sqlite3
from . import Command as rvc
from ..Utils import cd

###QC in one command
def quickrun(runtime_parameters):
    init_project(runtime_parameters)
    genotype_level(runtime_parameters)
    variant_level(runtime_parameters)
    sample_level(runtime_parameters)
    
##init project, import vcf and pheno
def init_project(runtime_parameters):
    workding_dir = '{}/1genotype_level'.format(runtime_parameters.get('general', 'maindir'))
    os.system('mkdir -p {}/tmp'.format(workding_dir))
    with cd(workding_dir):
        #init
        rvc.Init(Project='{}'.format(runtime_parameters.get('general', 'project'))).Run()
        #set sqlite parameters
        rvc.Admin(Set_runtime_option=['sqlite_pragma=synchronous=OFF,journal_mode=MEMORY',
                                      'temp_dir={}/tmp'.format(workding_dir)]).Run()
        #import genotypes
        rvc.Import(Input_files=[runtime_parameters.get('vtools', 'vcf')],
                   Format=runtime_parameters.get('vtools', 'format'),
                   Jobs=runtime_parameters.getint('vtools', 'jobs')).Run()
        #import phenotypes
        rvc.Phenotype(From_file=runtime_parameters.get('vtools', 'pheno')).Run()
        
##genotype level
def genotype_level(runtime_parameters):
    working_dir = '{}/1genotype_level'.format(runtime_parameters.get('general', 'maindir'))
    with cd(working_dir):
        #join dbSNP
        rvc.Use(Source='dbSNP').Run()
        #annotate with dbSNP
        rvc.Select(From_table='variant',
                   To_table=['variant_in', 'variants in dbSNP version 138'],
                   Condition=['dbSNP.chr IS NOT NULL']).Run()
        rvc.Compare(Tables=['variant', 'variant_in'],
                    Difference=['variant_out', 'variants not in dbSNP version 138']).Run()
        
#variant level
def variant_level(runtime_parameters):
    working_dir = '{}/2variant_level'.format(runtime_parameters.get('general', 'maindir'))
    mother_dir = '{}/1genotype_level'.format(runtime_parameters.get('general', 'maindir'))
    os.system('mkdir -p {}/tmp'.format(workding_dir))
    os.system('rsync -a {0}/{1}_genotype.DB {0}/{1}.proj {2}/'.format(mother_dir,
                                                                      runtime_parameters.get('general', 'project'),
                                                                      working_dir))
    gd = runtime_parameters.get('genotype', 'gd').split(',')
    gq = runtime_parameters.get('genotype', 'gq')
    jobs = runtime_parameters.get('vtools', 'jobs')
    with cd(working_dir):
        #0. set sqlite parameters
        rvc.Admin(Set_runtime_option=['sqlite_pragma=synchronous=OFF',
                                      'temp_dir={}/tmp'.format(workding_dir),
                                      'journal_mode=MEMORY']).Run() 
        #1. remove low quality genotypes
        rvc.Remove(Type='genotypes', Items=['GD<{} OR GD>={} OR GQ<{}'.format(gd[0], gd[1], gq)])
        #2. generate genotype stats
        rvc.Update(From_stat=['maf=maf()',
                              'GT=#(GT)',
                              'alt=#(alt)',
                              'hom=#(hom)',
                              'het=#(het)',
                              'other=#(other)',
                              'wtGT=#(wtGT)',
                              'mutGT=#(mutGT)',
                              'mis=#(missing)'],
                   Jobs=jobs).Run()
        #missing rate
        rvc.Update(Set=['mr=mis/((mis + GT) * 1.0)']).Run()
        #here is a hack, to loop over the batch, we must know the batches first
        conn = sqlite3.connect('{}.proj'.format(runtime_parameters.get('general', 'project')))
        c = conn.cursor()
        batch_column_name = runtime_parameters.get('pheno', 'batch')
        c.execute('SELECT DISTINCT {} FROM sample'.format(batch_column_name))
        batches = map(lambda x: x[0], c.fetchall())
        conn.close()
        for batch in batches:
            rvc.Update(From_stat=['maf_b{}=maf()'.format(batch),
                                  'GT_b{}=#(GT)'.format(batch),
                                  'alt_b{}=#(alt)'.format(batch),
                                  'hom_b{}=#(hom)'.format(batch),
                                  'het_b{}=#(het)'.format(batch),
                                  'other_b{}=#(other)'.format(batch),
                                  'wtGT_b{}=#(wtGT)'.format(batch),
                                  'mutGT_b{}=#(mutGT)'.format(batch),
                                  'mis_b{}=#(missing)'.format(batch)],
                       Samples=['{}={}'.format(batch_column_name, batch)],
                       Jobs=jobs).Run()
            #missing rate
            rvc.Update(Set=['mr_b{0}=mis_b{0}/((mis_b{0} + GT_b{0}) * 1.0)'.format(batch)], Jobs=jobs).Run()
        #hwe p value
        rvc.Update(Set=['hwe_pvalue=HWE_exact(GT, het, hom, other)'], Jobs=jobs).Run()
        #3. remove zombie variants, which are not variants any more because of the genotype level QC
        rvc.Select(From_table='variant',
                   To_table=['_snv1', 'mother: variant, after genotype level QC, some variants are gone with the removed genotypes'],
                   Condition=['GT!=wtGT', 'GT!=mutGT']).Run()
        #4. remove variants on sex and mitochondri chromosomes
        rvc.Select(From_table='_snv1',
                   To_table=['_snv2', 'mother: _snv1, remove variants on sex and mitochondrial chromosomes'],
                   Condition=['chr!="X"', 'chr!="Y"', 'chr!="M"']).Run()
        #5. remove variants with high missing rates
        mr = runtime_parameters.getfloat('variant', 'missing')
        mr_b = runtime_parameters.getfloat('variant', 'batch_missing')
        rvc.Select(From_table='_snv2',
                   To_Table=['_snv3',
                             'mother: _snv2, removed variants with high missing rate (>={}%)'.format(mr * 100)],
                   Condition=['mr<{}'.format(mr)]).Run()
        #6. remove high missing rate variants in any batch
        rvc.Select(From_table='_snv3',
                   To_table=['_snv4',
                             'mother: _snv3, removed variants with high missing rate (>={}%) in any batch'.format(mr_b * 100)],
                   Condition=map(lambda x: 'mr_b{}<{}'.format(x, mr_b), batches)).Run()
        #7. remove batch specific variants, which are novel, only in one batch, and not a singletom
        for batch in batches:
            rvc.Select(From_table='_snv4',
                       To_table=['_spec_b{}'.format(batch),
                                 'mother: _snv3, variants only in batch {}, not in dbSNP version 138, and not a singleton'.format(batch)],
                       Condition=['alt_b{}=alt'.format(batch), 'alt>1', 'dbSNP.chr IS NULL']).Run()
        #all batch specific variants
        rvc.Compare(Tables=map(lambda x: '_spec_b{}'.format(x), batches),
                    Union=['_spec_all',
                           'mother: _snv3, union of _spec_b[batch], all batch specific variants together']).Run()
        rvc.Compare(Tables=['_snv4', '_spec_all'],
                    Difference=['_snv5',
                                'mother: _snv4, removed batch specific variants in _spec_all. _snv5 = _snv4 - _spec_all']).Run()
        #8. remove variants hwe_pvalue < 0.001, for sample_level QC only
        hwe = runtime_parameters.get('variant', 'hwe_mds')
        rvc.Select(From_table='_snv5',
                   To_table=['_snv6',
                             'mother: _snv5, removed variants hwe_pvalue < 0.001, for sample level QC only'.format(hwe)],
                   Condition=['hwe_pvalue>={}'.format(hwe)]).Run()
                
#sample level
def sample_level(runtime_parameters):
    working_dir = '{}/3sample_level'.format(runtime_parameters.get('general', 'maindir'))
    mother_dir = '{}/2variant_level'.format(runtime_parameters.get('general', 'maindir'))
    os.system('mkdir -p {}/tmp'.format(workding_dir))
    os.system('rsync -a {0}/{1}_genotype.DB {0}/{1}.proj {2}/'.format(mother_dir,
                                                                      runtime_parameters.get('general', 'project'),
                                                                      working_dir))
    with cd(working_dir):
        #0. set sqlite parameters
        rvc.Admin(Set_runtime_option=['sqlite_pragma=synchronous=OFF',
                                      'temp_dir={}/tmp'.format(workding_dir),
                                      'journal_mode=MEMORY']).Run()
        #1. remove unrelevant variants, so we can compute sample level counts based on the variants we want
        rvc.Compare(Tables=['variant', '_snv6'],
                    Difference=['_snv6_c', 'mother: variant, the complement of _snv6']).Run()
        rvc.Remove(Type='variants', Items=['_snv6_c'])
        #total gt call for each individual
        jobs = runtime_parameters.get('vtools', 'jobs')
        rvc.Phenotype(From_stat=['_GT=#(GT)'], Jobs=jobs).Run()
        mds()

#to run the mds pipeline
def mds(From_table='variant', Maf=0.01):
    rvc.Select(From_table=From_table,
               To_table=['_mds',
                         'mother: variant, removed variants with maf =< {}%'.format(Maf * 100)],
               Condition=['maf>{}'.format(Maf)]).Run()
    rvc.Execute(Specfile='KING', Extra_args='--var_table _mds').Run()
        
        
###obsolete
class Variant:
    def __init__(self, name='variant', mother='variant', cond='1'):
        self.name = name
        self.mother = mother
        self.cond = cond
        self.create = 'vtools select {} "{}" -t {}'.format(mother, cond, name)
        
    def born(self, name, cond):
        return Variant(name, self.name, cond)

    def select(self, target, gcond=1, vcond=1, scond=1,):
        cmds = '''vtools select {} "{}" --genotypes "{}" --samples "{}" {}'''.format(self.name, vcond, gcond, scond, target)
        return cmds
    
    def update(self, cmd):
        return "vtools update {} {}".format(self.name, cmd)

    def stat(self, suffix='', gcond=1, scond=1):
        cmds = '''vtools update {3} --from_stat 'maf{0}=maf()' 'GT{0}=#(GT)' 'alt{0}=#(alt)' 'hom{0}=#(hom)' 'het{0}=#(het)' 'other{0}=#(other)' 'wtGT{0}=#(wtGT)' 'mutGT{0}=#(mutGT)' 'mis{0}=#(missing)' --genotypes "{1}" --samples "{2}"'''.format(suffix, gcond, scond, self.name)
        cmds = '\n'.join([cmds, '''vtools update {1} --set 'mr{0}=mis{0}/((mis{0} + GT{0}) * 1.0)' '''.format(suffix, self.name)])
        return cmds

class MDS:
    def __init__(self, var, maf, jobname='King', scond=1):
        self.tbl = jobname
        self.jobname = jobname
        self.cmds = '{}\n'.format(var.born(self.tbl, "{}>0.01 and chr<>'X'".format(maf)).create) + \
                    'vtools execute KING --jobname {} --var_table {} --samples "{}"'.format(jobname, self.tbl, scond)
        
    def plot_with(self, pheno, scond='1'):
        return 'vtools_report plot_pheno_fields {0}_MDS1 {0}_MDS2 {1} --samples "{2}" --dot {0}_with_{1}.pdf --discrete_color Accent'.format(self.jobname, pheno, scond)

class TiTv:
    def __init__(self, srcTbl='snv', suffix='Raw', gcond=1, vcond=1, scond=1):
        self.srcTable = srcTbl
        self.tiCond = "({}) AND ((ref='A' and alt='G') OR (ref='G' and alt='A') OR (ref='C' and alt='T') OR (ref='T' and alt='C'))".format(vcond)
        self.tvCond = "({}) AND ((ref='A' and alt='C') OR (ref='A' and alt='T') OR (ref='G' and alt='C') OR (ref='G' and alt='T') OR (ref='C' and alt='A') OR (ref='C' and alt='G') OR (ref='T' and alt='A') OR (ref='T' and alt='G'))".format(vcond)
        self.tiTable = "ti" + suffix
        self.tvTable = "tv" + suffix
        self.tiCol = 'ti' + suffix
        self.tvCol = 'tv' + suffix
        self.gcond = gcond
        self.titvCol = 'titv' + suffix
        self.scond = scond
        self.suffix = suffix
        self.template = 'echo ------------------------------------------\n' + \
                        'echo Now $(date) start titv for {10} from {0}\n' + \
                        'echo mean std min max overall\n' + \
                        'vtools select {0} "{1}" -t {3}\n' + \
                        'vtools select {0} "{2}" -t {4}\n' + \
                        'vtools phenotype --from_stat "{5}=#(alt)" --genotypes "(variant_id in (select variant_id from {3})) AND ({7})"\n' + \
                        'vtools phenotype --from_stat "{6}=#(alt)" --genotypes "(variant_id in (select variant_id from {4})) AND ({7})"\n' + \
                        'vtools phenotype --set "{8}={5}/({6} * 1.0)"\n' + \
                        'vtools phenotype -s "{9}" --output "avg({8})" "stdev({8})" "min({8})" "max({8})" "sum({5})/(sum({6}) * 1.0)"\n' + \
                        'echo Now $(date) finish titv for {10} from {0}\n' + \
                        'echo ------------------------------------------\n'
        self.cmds = ''

    def __make(self):
        self.cmds = self.template.format(self.srcTable, self.tiCond, self.tvCond, self.tiTable, self.tvTable, self.tiCol, self.tvCol,
                                         self.gcond, self.titvCol, self.scond, self.suffix)        

    def getcmds(self):
        if self.cmds:
            return self.cmds
        else:
            self.__make()
            return self.cmds

    def plot(self, scond=1):
        cmds = 'vtools_report plot_pheno_fields {} '
