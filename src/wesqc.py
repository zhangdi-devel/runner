#!/usr/bin/env python
import os, sys, ConfigParser
from argparse import ArgumentParser
import Runner.Vtools.Command as rvc
import Runner.Vtools.Pipeline as rvp

class runtime_parameters():
    def __init__(self, config_file):
        self.config_file = config_file
        self.config = ConfigParser.SafeConfigParser()
        self.config.read(self.config_file)
        
    def check(self):
        ##the main dir must exists
        if os.path.isdir(self.config.get('general', 'maindir')):
            sys.stdout.write('The maindir is: "{}"'.format(self.config.get('general', 'maindir')))
        else:
            sys.stderr.write('maindir "{}" does not exists\n'.format(self.config.get('general', 'maindir')))
            return False
        ##vtools input files
        for k, v in map(lambda x: (x, self.config.get('vtools', x)), ['vcf', 'format', 'pheno']):
            if os.path.isfile(v):
                sys.stdout.write('vtools input file {}: {}\n'.format(k, v))
            else:
                sys.stdout.write('vtools input file {}: "{}" does not exists\n'.format(k, v))
                return False
        return True

def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--conf', default='wesqc.conf', metavar='File',
                        help='set config file, default to wesqc.conf')
    args = parser.parse_args()
    if os.path.isfile(args.conf):
        rp = runtime_parameters(args.conf)
        if rp.check():
            rvp.quickrun(rp)
        else:
            sys.stderr.write('runtime parameter(s) error\n')
            return False
    else:
        os.stderr.write('{} not exists, please prepare config file before run.'.format(args.conf))
        return False
    return True    
    
if __name__ == '__main__':
    main()
