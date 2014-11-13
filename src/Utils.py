#!/usr/bin/python

import sys, os, subprocess, glob, re, shlex
from contextlib import contextmanager
from multiprocessing import Pool, Process, Queue, Lock, Value, cpu_count
from collections import OrderedDict, defaultdict, Counter
from distutils.file_util import copy_file
from shutil import move
from os.path import basename
import mysql.connector as mc

class Region:
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end

    def len(self):
        return self.end - self.start + 1

    def seq(self):
        return ''

class Gene(region):
    def __init__(self, chr, start, end, strand, name):
        region.__init__(self, chr, start, end, strand)
        self.name = name
        self.exons = []

    def add_exon(self, exon):
        self.exons.append(exon)


class Exon(region):
    def __init__(self, chr, start, end, strand, gene=''):
        region.__init__(self, chr, start, end, strand)
        self.gene = gene

###
# Conf
###
#global configuration
class Conf:
    def __init__(self):
        #for database connection
        self.dbhost = 'hp-nas2-ib'
        self.dbport = 3306        
        self.dbuser = 'csg'
        self.dbpass = ''
        self.dbname = 'hp19_seq'
        self.__width_cache = 1
        #general
        self.path = {'PATH':os.environ['PATH']}
        self.home = os.getenv('HOME', '/mnt/nfs/zhangdi')
        #bwa
        self.RefFasta = '/mnt/nfs/zhangdi/ref/human/g1k_phase1/fasta/human_g1k_v37.fasta'
        self.Aligner = 'bwa'
        self.BWAindex = '/mnt/nfs/zhangdi/ref/human/g1k_phase1/bwa/human_g1k_v37'
        #GATK
        self.GATKPath = os.path.join(self.home, 'software/gatk-2.6.5')
        self.PICARDPath = os.path.join(self.home, 'software/picard-tools-1.94')
        self.samtoolsPath = os.path.join(self.home, 'software/samtools-0.1.19')
        self.threads = 1
        self.MemLimit = '8G'
        #coverage
        self.mut = {'A': ['T', 'G', 'G', 'G', 'C', 'G', 'G', 'G'],
                    'T': ['A', 'C', 'C', 'C', 'G', 'C', 'C', 'C'],
                    'C': ['A', 'T', 'T', 'T', 'G', 'T', 'T', 'T'],
                    'G': ['C', 'A', 'A', 'A', 'T', 'A', 'A', 'A']}

    def log(self, msg = None, flush=False):
        #if self.debug or self.quiet:
        #    return
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = "{0:{width}}".format('\r', width = self.__width_cache + 10) + "\r" if flush else ''
        end = '' if flush else '\n'
        start = '\n' + start if msg.startswith('\n') else start
        end = end + '\n' if msg.endswith('\n') else end
        msg = msg.strip()
        sys.stderr.write(start + "\033[1;40;32mMESSAGE: {}\033[0m".format(msg) + end)
        self.__width_cache = len(msg)

def cov(bam):
    samFile = pysam.Samfile(bam, "rb")
    c = []
    for exon in conf().gene.exons:
        b = 0
        for pc in samFile.pileup(exon.chr, exon.start, exon.end):
            #print 'coverage at base {} = {}'.format(pc.pos, pc.n)
            if pc.pos >= exon.start and pc.pos <= exon.end and pc.n >= 8:
                b += 1
            #print 'coverage {}/{} = {}'.format(b, exon.len(), b/(exon.len() * 1.0))
        c.append(b)
    return c



###
# Utility function / classes
###

class StdoutCapturer(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

@contextmanager
def stdoutRedirect(to=os.devnull):
    '''
    import os

    with stdoutRedirect(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    '''
    fd = sys.stdout.fileno()

    ##### assert that Python and C stdio write using the same file descriptor
    ####assert libc.fileno(ctypes.c_void_p.in_dll(libc, "stdout")) == fd == 1

    def _redirect_stdout(to):
        sys.stdout.close() # + implicit flush()
        os.dup2(to.fileno(), fd) # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, 'w') # Python writes to fd

    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        with open(to, 'a') as file:
            _redirect_stdout(to=file)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout) # restore stdout.
                                            # buffering and flags such as
                                            # CLOEXEC may be different

#http://stackoverflow.com/a/13197763, by Brian M. Hunt
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def runCommand(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,instream = None, msg = '', upon_succ = None, show_stderr = False, return_zero = True, shell=False):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    popen_env = os.environ.copy()
    popen_env.update(env.path)
    try:
        tc = subprocess.Popen(cmd, stdin = subprocess.PIPE,
                              stdout = stdout, stderr = stderr,
                              env=popen_env, shell = shell)
        if instream:
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if return_zero:
            if tc.returncode < 0:
                raise ValueError ("Command '{0}' was terminated by signal {1}".format(cmd, -tc.returncode))
            elif tc.returncode > 0:
                raise ValueError ("{0}".format(error))
        if error.strip() and show_stderr:
            env.log(error)
    except OSError as e:
        raise OSError ("Execution of command '{0}' failed: {1}".format(cmd, e))
    # everything is OK
    if upon_succ:
        # call the function (upon_succ) using others as parameters.
        upon_succ[0](*(upon_succ[1:]))
    return out, error

class CMDWorker(Process):
    def __init__(self, queue):
        Process.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            try:
                cmd = self.queue.get()
                if cmd is None:
                    break
                else:
                    runCommand(cmd)
            except KeyboardInterrupt:
                break
            
def runCommands(cmds, ncpu):
    try:
        jobs = []
        queue = Queue()
        for i in cmds:
            queue.put(i)
        for i in range(ncpu):
            p = CMDWorker(queue)
            p.start()
            jobs.append(p)
            queue.put(None)
        for j in jobs:
            j.join()
    except KeyboardInterrupt:
        raise ValueError('Commands terminated!')

#utils that allow lambda function in mutilprocessing map
#http://stackoverflow.com/a/16071616 by klaus-se
def parmap(f, X, nprocs = cpu_count()):
    def spawn(f):
        def fun(q_in,q_out):
            while True:
                i,x = q_in.get()
                if i is None:
                    break
                q_out.put((i,f(x)))
        return fun
    #
    q_in   = Queue(1)
    q_out  = Queue()
    proc = [Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]
    [p.join() for p in proc]
    return [x for i,x in sorted(res)]
