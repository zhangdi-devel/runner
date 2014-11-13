#!/usr/bin/python
from Utils import *
from random import choice
#from Conf import *
        
class GemReads:
    def __init__(self, dp, sample_id):
        self.ref = env.GenePath
        self.dp = dp
        self.sample_id = sample_id
        self.num = 1 + (dp * env.GeneLength ) / (env.ReadLength * (2 if env.PairedEnd else 1)) 
        self.len = env.ReadLength
        self.qual = env.BaseQualOffset
        self.mod = env.GemModel
        self.fsize = env.FragmentLength
        self.fsizeSD = env.FragmentLengthSD
        self.pair = '-p' if env.PairedEnd else ''
        self.ReadsDir = os.path.join(env.ReadsDir, 'DP'+str(dp))
        self.outbase = os.path.join(self.ReadsDir, 'sample'+str(sample_id))
        runCommand('mkdir -p {}'.format(self.ReadsDir))
        self.SNPDir = os.path.join(env.SNPDir, 'DP'+str(dp))
        self.SNP = os.path.join(self.SNPDir, 'sample'+str(sample_id)+'.txt')
        runCommand('mkdir -p {}'.format(self.SNPDir))
        self.Read1 = self.outbase + '_fir.fastq'
        self.Read2 = self.outbase + '_sec.fastq'
        self.exon_region = reduce(lambda x,y: x+y, map(lambda e: range(e.start, e.end + 1), env.gene.exons))

    def generateSNP(self):
        with cd(self.SNPDir):
            with open(env.GenePath, 'r') as gf, open(basename(self.SNP), 'w') as sf:
                sf.write('0.5\n')
                pos = choice(self.exon_region) - env.offset
                name, seq = '', ''
                for line in gf:
                    if line.startswith('>'):
                        name = line.strip()[1:]
                    else:
                        seq += line.strip()
                Ref = seq[pos]
                Alt = choice(env.mut[Ref])
                ##GemReads uses 1-based position
                sf.write('0.5\t{}\t{}\t{}\n'.format(name, pos + 1, Alt))

    def run(self):
        self.generateSNP()
        cmd = ("GemReads.py -r {} -n {} -g {} -l {} -m {} -q {} " + \
              "-o {}").format(env.GeneLink, self.num, self.SNP, env.ReadLength,
                             env.GemModel, env.BaseQualOffset, self.outbase) + \
              " -u {} -s {} -p".format(env.FragmentLength, env.FragmentLengthSD) \
              if env.PairedEnd else ""
        #env.log(cmd)
        runCommand(cmd)

class BWA:
    def __init__(self, reads):
        dp = reads.dp
        sample_id = reads.sample_id
        self.reads = reads
        self.Index = env.BWAindex
        self.Read1 = reads.Read1
        self.Read2 = reads.Read2
        self.SAI1 = '3bam/DP{}/sample{}_fir.sai'.format(dp, sample_id)
        self.SAI2 = '3bam/DP{}/sample{}_sec.sai'.format(dp, sample_id)
        self.SAM = '3bam/DP{}/sample{}.sam'.format(dp, sample_id)
        self.SAMwithRG = '3bam/DP{}/sample{}.RG.sam'.format(dp,sample_id)
        self.BAM = '3bam/DP{}/sample{}.bam'.format(dp, sample_id)
        self.Sorted = '3bam/DP{}/sample{}.sorted'.format(dp, sample_id)
        self.SortedBAM = '3bam/DP{}/sample{}.sorted.bam'.format(dp, sample_id)
        #self.addRG = '\'if (m/^\@PG/) {print "@RG\tID:1\tSM:1\n";} if ($_ !~ m/^\@/) {chomp; print "\tRG:Z:1"}\''
        runCommand('mkdir -p 3bam/DP{}'.format(dp))
        

    def run(self):
        cmds = []
        cmds.append('bwa aln {} {} > {} 2> /dev/null'.format(self.Index, self.Read1, self.SAI1))
        cmds.append('bwa aln {} {} > {} 2> /dev/null'.format(self.Index, self.Read2, self.SAI2))
        cmds.append('bwa sampe {} {} {} {} {} > {} 2> /dev/null'.format(self.Index, self.SAI1, self.SAI2, self.Read1, self.Read2, self.SAM))
        for cmd in cmds:
            #env.log(cmd)            
            os.system(cmd)
       
class samtools:
    def __init__(self, bwa):
        self.SAM = bwa.SAM
        self.bwa = bwa
        self.SAMwithRG = bwa.SAMwithRG
        self.BAM = bwa.BAM
        self.Sorted = bwa.Sorted
        self.SortedBAM = bwa.SortedBAM
        self.VCFdir = '5vcf/DP{}'.format(bwa.reads.dp)
        self.outbase =  os.path.join(self.VCFdir, 'sample{}'.format(bwa.reads.sample_id))
        self.rawBCF = self.outbase + '.raw.bcf'
        self.fltVCF = self.outbase + '.flt.vcf'
        runCommand('mkdir -p {}'.format(self.VCFdir))

    def addRG(self):
        with open(self.SAM, 'r') as s, open(self.SAMwithRG, 'w') as rg:
            for line in s:
                if line.startswith('@PG'):
                    rg.write('{}@RG\tID:1\tSM:1\n'.format(line))
                elif line.startswith('@'):
                    rg.write(line)
                else:
                    rg.write(line.strip())
                    rg.write('\tRG:Z:1\n')

    def run(self):
        cmds = []
        self.addRG()
        cmds.append('samtools view -Sb {} > {} 2> /dev/null'.format(self.SAMwithRG, self.BAM))
        cmds.append('samtools sort {} {} 2> /dev/null'.format(self.BAM, self.Sorted))
        cmds.append('samtools index {} 2> /dev/null'.format(self.SortedBAM))
        cmds.append('samtools mpileup -uf {} {} 2> /dev/null| bcftools view -bvcg - > {} 2> /dev/null'.format(env.RefFasta, self.SortedBAM, self.rawBCF))
        cmds.append('bcftools view {} 2> /dev/null | vcfutils.pl varFilter -D 1000 > {} 2> /dev/null'.format(self.rawBCF, self.fltVCF))
        for cmd in cmds:
            #env.log(cmd)            
            os.system(cmd)
        

class GATK:
    def __init__(self, bwa):
        self.bwa = bwa
        self.MarkDupJar = os.path.join(env.PICARDPath, 'MarkDuplicates.jar')
        self.GATKJar = os.path.join(env.GATKPath, 'GenomeAnalysisTK.jar')
        self.rawBAM = bwa.SortedBAM
        self.outbase = '4gatk/DP{}/sample{}'.format(bwa.reads.dp, bwa.reads.sample_id)
        self.MarkDupBAM = self.outbase + 'MarkDup.bam'
        self.rawVCF = self.outbase + 'both.raw.vcf'
        self.filteredVCF = self.outbase + 'filtered.vcf'
        self.SNV = self.outbase + 'snv.vcf'
        runCommand('mkdir -p 4gatk/DP{}'.format(bwa.reads.dp))
        
#        self.MemLimit = env.MemLimit

    def MarkDup(self):
        cmd = 'java -Xmx{} -XX:ParallelGCThreads={} -jar {} '.format(env.MemLimit, env.threads, self.MarkDupJar) + \
              'I={} O={} M={} VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true'.format(self.rawBAM, self.MarkDupBAM, self.outbase + 'MarkDup.metrics')
        runCommand(cmd)
  
    def samIndex(self):
        cmd = 'samtools index {}'.format(self.MarkDupBAM)
        runCommand(cmd)

    def UnifiedGenotyper(self, bam):
        cmd = 'java -Xmx{} -XX:ParallelGCThreads={} -jar {} '.format(env.MemLimit, env.threads, self.GATKJar) + \
              '-R {} '.format(env.RefFasta) + \
              '-T UnifiedGenotyper ' + \
              '-nt {} '.format(env.threads) + \
              '-I {} '.format(bam) + \
              '-stand_call_conf 50 ' + \
              '-stand_emit_conf 50 ' + \
              '-mbq 20 ' + \
              '-A QualByDepth ' + \
              '-A HaplotypeScore ' + \
              '-A FisherStrand ' + \
              '-A MappingQualityZero ' + \
              '-A RMSMappingQuality ' + \
              '-A MappingQualityRankSumTest ' +\
              '-A ReadPosRankSumTest ' + \
              '-A BaseQualityRankSumTest ' + \
              '-o {} ' \
              '-glm BOTH'
        runCommand(cmd, show_stderr=True)

    def VariantFiltration(self):
        cmd = 'java -Xmx{} -XX:ParallelGCThreads={} -jar {} '.format(env.MemLimit, env.threads, self.GATKJar) + \
              '-T VariantFiltration ' + \
              '-R {} '.format(env.RefFasta) + \
              '-o {} '.format(self.filteredVCF) + \
              '--variant {} '.format(self.rawVCF) + \
              '--filterExpression "AB < 0.2 || MQ0 > 50" ' + \
              '--filterName "Nov09filters" ' + \
              '--mask mask.vcf ' + \
              '--maskName InDel'
        runCommand(cmd, show_stderr=True)
        
    def SelectSNV(self):
        cmd = 'java -Xmx{} -XX:ParallelGCThreads={} -jar {} '.format(env.MemLimit, env.threads, self.GATKJar) + \
              '-R {} '.format(env.RefFasta) + \
              '-T SelectVariants ' + \
              '--variant {} '.format(self.filteredVCF) + \
              '-o {} '.format(self.SNV) + \
              '-selectType SNP ' + \
              '-select "DP > 10 && QD > 10.0 && MQ > 50.0 && FS < 13.0 && HaplotypeScore < 13.0 && MQRankSum > -1.96 && ReadPosRankSum > -1.96"'
        runCommand(cmd, show_stderr=True)

    def run(self):
        self.MarkDup()
        self.samIndex()
        self.UnifiedGenotyper(self.MarkDupBAM)
        self.VariantFiltration()
        self.SelectSNV()

def runSample(dp,n1, n2):
    for i in range(n1, n2):
        if (i+1)%10 == 0:
            env.log('{} completed for depth={} ids={}..{}'.format(i-n1+1, dp, n1, n2))
        reads = GemReads(dp, i+1)
        reads.run()
        bwa = BWA(reads)
        bwa.run()
        vcf = samtools(bwa)
        vcf.run()
                

if __name__ == '__main__':
    dps = 20 #[5, 10, 15, 20, 30, 40, 50, 60, 80, 100]
    #parmap(lambda x: runSample(x, 2), dps, 1)
    n1, n2 = map(int, sys.argv[1:3])
    p = (n2 - n1)/8
    parmap(lambda x: runSample(20, n1 + x*p, n1 + (x+1)*p), range(8), 8)
