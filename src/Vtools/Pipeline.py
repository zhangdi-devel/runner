#!/usr/bin/python
import sys, os
from . import Command as rvc
from ..Utils import cd
###project
def quickrun(runtime_parameters):
    init_project(runtime_parameters)
    
##init project, import vcf and pheno
def init_project(runtime_parameters):
    workding_dir = '{}/1genotype_level'.format(runtime_parameters.get('general', 'maindir'))
    os.system('mkdir -p {}/tmp'.format(workding_dir))
    with cd(workding_dir):
        #init
        rvc.Init(Project='{}'.format(runtime_parameters.get('general', 'project'))).Run()
        #set sqlite parameters
        rvc.Admin(Set_runtime_option=['sqlite_pragma=synchronous=OFF',
                                      'temp_dir={}/tmp'.format(workding_dir),
                                      'journal_mode=MEMORY']).Run()
        #import genotypes
        rvc.Import(Input_files=runtime_parameters.get('vtools', 'vcf'),
                   Format=runtime_parameters.get('vtools', 'format'),
                   Jobs=runtime_parameters.get('vtools', 'jobs')).Run()
        #import phenotypes
        rvc.Phenotype(From_file=runtime_parameters.get('general', 'pheno')).Run()
        
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
                              'mis=#(missing)'], Jobs=jobs).Run()
        #here is a hack, to loop over the batch, we must know the max(batch) first
        
        for batch in range(1,5):
            rvc.Update(From_stat=['mafQC_b{}=maf()'.format(batch),
                                  'GTQC_b{}=#(GT)'.format(batch),
                                  'altQC_b{}=#(alt)'.format(batch),
                                  'homQC_b{}=#(hom)'.format(batch),
                                  'hetQC_b{}=#(het)'.format(batch),
                                  'otherQC_b{}=#(other)'.format(batch),
                                  'wtGTQC_b{}=#(wtGT)'.format(batch),
                                  'mutGTQC_b{}=#(mutGT)'.format(batch),
                                  'misQC_b{}=#(missing)'.format(batch)],
                       Samples=['batch={}'.format(batch)]).Run()
        rvc.Update(Set=['mrQC_b{0}=misQC_b{0}/((misQC_b{0} + GTQC_b{0}) * 1.0)'.format(batch)]).Run()
        #total gt call for each individual
        rvc.Phenotype(From_stat=['num=count(*)'], Jobs=jobs).Run()
        #missing rate
        rvc.Update(Set=['mr=mis/((mis + GT) * 1.0)']).Run()
        #hwe p value
        rvc.Update(Table='variant', Set=['hwe_pvalue=HWE_exact(GT, het, hom, other)']).Run()
        #remove zombie variants, which are not variants any more because of the genotype level QC
        rvc.Select(From_table='variant',
                   To_table=['_snv1', 'variants after genotype level QC, some variants are gone with the removed genotypes'],
                   Condition=['GT!=wtGT', 'GT!=mutGT']).Run()
        
#sample level



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
