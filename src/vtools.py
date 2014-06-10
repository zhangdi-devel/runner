#!/usr/bin/python
import sys, os


#class Genotype:

#class Sample:

class Variant:
    def __init__(self, name, mother='variant', cond='1'):
        self.name = name
        self.mother = mother
        self.cond = cond
        self.create = 'vtools select {} "{}" -t {}'.format(mother, cond, name)
        
    def born(self, name, cond):
        return Variant(name, self.name, cond)
    
    def update(self, cmd, cond):
        return "vtools update {} {}".format(self.name, cond)

#class MDS:
    


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
                        'vtools phenotype --from_stat "{5}=#(alt)" -g "(variant_id in (select variant_id from {3})) AND ({7})"\n' + \
                        'vtools phenotype --from_stat "{6}=#(alt)" -g "(variant_id in (select variant_id from {4})) AND ({7})"\n' + \
                        'vtools phenotype --set "{8}={5}/({6} * 1.0)"\n' + \
                        'vtools phenotype -s "{9}" --output "avg({8})" "stdev({8})" "min({8})" "max({8})" "sum({5})/sum({6})"\n' + \
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
