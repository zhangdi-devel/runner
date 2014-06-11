#!/usr/bin/python

import ConfigParser
#from argparse import ArgumentParser
from Runner.Core import main
#from Runner.Utils import env

class Conf:
    def __init__(self):
        #self.parser = ArgumentParser()
        #group = self.parser.add_mutually_exclusive_group()
        #group.add_argument('-c', '--conf', metavar='STR'
        #                   help='''Configuration file''')
        #group.add_argument('conf', metavar='STR'
        #                   help='''Configuration file''')
        self.conf = ConfigParser.ConfigParser().read('Runner.conf')
        
    #def get(self):
    #    args = self.parser.parse_args()
    #    self.conf.read(args.conf)
    #    return self.conf

if __name__ == '__main__':
    try:
        conf = Conf()
        main(conf)
    except Exception as e:
        raise
        env.eorror('{}'.format(e))
