#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    read = ""
    for line in gzip.open(amplicon_file, "rt"):
        if line.startswith(">"):
            if len(read) > minseqlen:
                yield read
                del(read)
                read = ""
            else:
                del(read)
                read = ""
            continue
        read += line.rstrip("\n")
    if len(read) > minseqlen:
        yield read


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    seq_dict = {}
    for seq in read_fasta(amplicon_file, minseqlen):
        seq_dict.setdefault(seq, 0)
        seq_dict[seq] += 1

    for seq, count in sorted(seq_dict.items()):
        if count > mincount:
            yield [seq, count]


def abundance_greedy_clustering(amplicon_file,
                                minseqlen,
                                mincount,
                                chunk_size,
                                kmer_size):
    matrix = os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))
    otu_list = []
    for seq, count in dereplication_fulllength(amplicon_file,
                                               minseqlen,
                                               mincount):
        # should be done with next()
        if len(otu_list) == 0:
            otu_list.append([seq,count])
        else :
            otu_status = True
            for (otu, occ_in_list) in otu_list:
                if get_identity(nw.global_align(seq,
                                                otu,
                                                gap_open=-1,
                                                gap_extend=-1,
                                                matrix=matrix)) >= 97.0:
                    otu_status = False
                    break
            if otu_status :
                otu_list.append([seq, count])
    return otu_list


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for kmer in cut_kmer(sequence, kmer_size):
        kmer_dict.setdefault(kmer, [])
        kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates (kmer_dict, sequence, kmer_size):
    common = []     
    
    common.append(kmer_dict[kmer] for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict.keys())
    flatten_list = [seq for sublist in common for item in sublist for seq in item]

    return [seq[0] for seq in Counter(flatten_list).most_common(2)]


def detect_chimera(perc_identity_matrix):
    chimera = False
    status =0

    if statistics.mean([statistics.pstdev(segment) for segment in perc_identity_matrix]) > 5:
        for segment in perc_identity_matrix:
            status += segment[0] > segment[1]
    
        if status != 0 and status != len(perc_identity_matrix):
            chimera = True
        
    return chimera


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, 'w+') as o_file:
        i = 1
        for OTU in OTU_list:
            o_file.write('>OTU_%i occurrence:%i\n%s\n' % (i,
                                                      OTU[1],
                                                      fill(OTU[0])))
            i += 1

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici


if __name__ == '__main__':
    main()
