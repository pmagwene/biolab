#!/usr/bin/env python
"""
nuc2prot.py -- a Python module for translating nucleotide sequences to 
"""

import sys
import argparse
import collections, itertools

stdcodons_dna_table = { 
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 
'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 
'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 
'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 
'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 
'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 
'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 
'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 
'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 
'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*' ,'---': '-'}

AA_to_stdcodons_table = collections.defaultdict(list)
for codon,AA in stdcodons_dna_table.iteritems():
    AA_to_stdcodons_table[AA].append(codon)


ambiguous_dna_table = {
 'A': 'A',
 'B': 'CGT',
 'C': 'C',
 'D': 'AGT',
 'G': 'G',
 'H': 'ACT',
 'K': 'GT',
 'M': 'AC',
 'N': 'GATC',
 'R': 'AG',
 'S': 'CG',
 'T': 'T',
 'V': 'ACG',
 'W': 'AT',
 'X': 'GATC',
 'Y': 'CT'
}


def to_codons(nucseq, removegaps=False):
    """Gives the codons of a sequence, truncating any incomplete codons """
    if removegaps:
        nucseq = nucseq.replace('-','')
    n = len(nucseq)
    return [nucseq[i:i+3] for i in range(0,n-n%3,3)]
    
def translate_codons(codons, tbl=stdcodons_dna_table):
    AAs = []
    for codon in codons:
        if codon in tbl:
            AAs.append(tbl[codon])
        else:
            AAs.append('?')
    return AAs
    

def translate_dna(dnaseq, tbl=stdcodons_dna_table, removegaps=False):
    codons = to_codons(dnaseq, removegaps)
    AAs = translate_codons(codons, tbl)
    return ''.join(AAs)     



def thread_nucs(aaseq,nucseq,tbl=AA_to_stdcodons_table):
    codons = to_codons(nucseq)
    threaded = []
    offset = 0
    for i, AA in enumerate(aaseq):
        if AA != '-':
            codon = codons[i - offset]
            if codon in tbl[AA]:
                threaded.append(codon)
            else:
                raise Exception("Nucleotide threading error. AA/codon mismatch")
        else:
            threaded.append('---')
            offset += 1
    return ''.join(threaded)



# ------------------------------------------------------------------------
# The FASTA file code below has the following license:
#
# Copyright (C) 2003, 2004, 2006 by  Thomas Mailund <mailund@birc.au.dk>
# With some additions and corrections by P. Magwene
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
# USA.


class MalformedInput:
    "Exception raised when the input file does not look like a fasta file."
    pass

class FastaRecord:
    "Wrapper around a fasta record."
    def __init__(self, header, sequence):
        "Create a record with the given header and sequence."
        self.header = header
        self.sequence = sequence
        self.identifier = header.split()[0]        

    def __str__(self):
        result = ['>'+self.header]
        for i in xrange(0,len(self.sequence),60):
            result.append(self.sequence[i:i+60])
        return '\n'.join(result)
        

def _fasta_itr_from_file(file):
    "Provide an iteration through the fasta records in file."

    h = file.readline().strip()
    if h[0] != '>':
        raise MalformedInput()
    h = h[1:]

    seq = []
    for line in file:
        line = line.strip() # remove newline

        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))

            h = line[1:]
            seq = []
            continue

        seq += [line]

    yield FastaRecord(h,''.join(seq))


def _fasta_itr_from_name(fname):
    "Provide an iteration through the fasta records in the file named fname. "
    f = open(fname, 'rU')
    for rec in _fasta_itr_from_file(f):
        yield rec
    f.close()


def _fasta_itr(src):
    """Provide an iteration through the fasta records in file `src'.
    
    Here `src' can be either a file object or the name of a file.
    """
    if type(src) == str:
        return _fasta_itr_from_name(src)
    elif type(src) == file:
        return _fasta_itr_from_file(src)
    else:
        raise TypeError

def fasta_get_by_name(itr,name):
    "Return the record in itr with the given name."
    x = name.strip()
    for rec in itr:
        if rec.header.strip() == x:
            return rec
    return None

class fasta_itr:
    "An iterator through a sequence of fasta records."
    def __init__(self,src):
        "Create an iterator through the records in src."
        self.__itr = _fasta_itr(src)

    def __iter__(self):
        return self
    def next(self):
        return self.__itr.next()

    def __getitem__(self,name):
        return fasta_get_by_name(iter(self),name)

class fasta_slice:
    """Provide an iteration through the fasta records in file `src', from
    index `start' to index `stop'.

    Here `src' can be either a file object or the name of a file.
    """
    def __init__(self, src, start, stop):
        """Provide an iteration through the fasta records in file `src', from
        index `start' to index `stop'.

        Here `src' can be either a file object or the name of a file.
        """
        self.__itr = _fasta_itr(src)
        self.__current = 0
        self.__start = start
        self.__stop = stop

    def __iter__(self):
        return self

    def next(self):
        while self.__current < self.__start:
            # skip past first records until we get to `start'
            self.__itr.next()
            self.__current += 1

        if self.__current >= self.__stop:
            # stop after `stop'
            raise StopIteration

        self.__current += 1
        return self.__itr.next()

    def __getitem__(self,name):
        return fasta_get_by_name(iter(self),name)

def get_sequence(src,name):
    "Return the record in src with the given name."
    return fasta_itr(src)[name]


##-------------------------------------------------------
## Some additional FASTA related function by P. Magwene

class fasta_dict(dict):
    """Create a dictionary from a fasta file.
    
    dictionary provides a mapping from identifier --> FastaRecord
    """
    def __init__(self, src):
        self.__itr = _fasta_itr(src)
        for each in self.__itr:
            self.__setitem__(each.identifier, each)
            
def _fasta_itr_from_str(s):
    lines = s.splitlines()
    h = lines[0]
    if h[0] != '>':
        raise MalformedInput()
    h = h[1:]
    seq = []
    for line in lines[1:]:
        line = line.strip()
        if not len(line):
            continue
        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            h = line[1:]
            seq = []
            continue
        seq += [line]
    yield FastaRecord(h,''.join(seq))            

class fasta_dict_from_str(dict):
    def __init__(self, s):
        self.__itr = _fasta_itr_from_str(s)
        for each in self.__itr:
            self.__setitem__(each.identifier, each)


def fasta_records_as_string(recs, llen=70):
    s = ''
    for r in recs:
        s += ">%s\n" % r.header
        rstr = str(r)
        for i in range(0,len(r.sequence),llen):
            s += r.sequence[i:i+llen] + "\n"
    return s            

                    
def write_fasta_records(recs, outfile, llen=70):
    s = fasta_records_as_string(recs, llen)
    f = open(outfile, 'w')
    f.write(s)
    f.close()


def translate_fasta_file(fname, removegaps=True):
    nucrecs = [i for i in fasta_itr(fname)]
    protrecs = []
    for rec in nucrecs:
        protseq = translate_dna(rec.sequence, removegaps=removegaps)
        newrec = FastaRecord(rec.header, protseq)
        protrecs.append(newrec)
    return protrecs

        
def thread_fasta_nucs_on_prots(nucfile, protfile):
    nucrecs = [i for i in fasta_itr(nucfile)]
    protrecs = [i for i in fasta_itr(protfile)]
    if len(nucrecs) != len(protrecs):
        raise Exception("Nucleotide and protein file must contain same number of sequences.")
    threadedrecs = []
    for i in range(len(protrecs)):
        protseq = protrecs[i].sequence
        nucseq = nucrecs[i].sequence
        tseq = thread_nucs(protseq,nucseq)
        trec = FastaRecord(nucrecs[i].header, tseq)
        threadedrecs.append(trec)
    return threadedrecs

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Translate FASTA nucleotide seqs to protein sequences.")   
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)    
    parser.add_argument('-o', '--outfile',nargs='?', type=argparse.FileType('w'), default=sys.stdout)

    parser.add_argument('--action', type=str, default='translate',
                        choices=['translate','thread'],
                        help='Type of action to perform (default=translate).')
    parser.add_argument('--protfile', type=argparse.FileType('r'), required=False,
                        help="Protein FASTA file to be used for threading routine.")
                        
    parser.add_argument('--removegaps', action="store_true", default=False)

    args = parser.parse_args()
    
    if args.action == 'translate':
        print args.removegaps
        protrecs = translate_fasta_file(args.infile, removegaps=args.removegaps)
        recstr = fasta_records_as_string(protrecs)
        args.outfile.write(recstr)
        
    elif args.action == 'thread':
        if args.protfile is None:
             raise parser.error("Must specify a protein file for threading using --protfile argument.")
        threadedrecs = thread_fasta_nucs_on_prots(args.infile, args.protfile)
        recstr = fasta_records_as_string(threadedrecs)
        args.outfile.write(recstr)        

    