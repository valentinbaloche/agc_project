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

__author__ = "Valentin Baloche & Alix de Thoisy"
__copyright__ = "Universite de Paris"
__credits__ = ["Valentin Baloche & Alix de Thoisy"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Valentin Baloche & Alix de Thoisy"
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
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
                        type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen',
                        type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount',
                        type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size',
                        type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size',
                        type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """Reads a fasta file and returns seq whose length < minseqlen
    Parameters:
        amplicon_file (str): path to fasta file
        minseqlen (int): minimum length
    
    Yields:
        str : sequence on fasta file
    """
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
    """
    Parameters:
        amplicon_file (str): fasta file path
        minseqlen (int): minimum length
    
    Yields:
        (str, int): unique seq and number of occurences
    """
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
    """Returns OTU list with seq having less than 97% similiraties
    with the ones in the list.
    """
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
    """Updates k-mer dictionaries
    """
    for kmer in cut_kmer(sequence, kmer_size):
        kmer_dict.setdefault(kmer, [])
        kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates (kmer_dict, sequence, kmer_size):
    """Search for sequence's parents in k-mer dictionary
    """
    common = []
    
    common.append(kmer_dict[kmer] for kmer in cut_kmer(sequence, kmer_size) \
                                                if kmer in kmer_dict.keys())
    flatten_list = [seq for sublist in common for item in sublist \
                                                    for seq in item]

    return [seq[0] for seq in Counter(flatten_list).most_common(2)]


def detect_chimera(perc_identity_matrix):
    """Whether sequence is a chimera
    """
    chimera = False
    status =0

    if statistics.mean([statistics.pstdev(segment) for segment \
                                    in perc_identity_matrix]) > 5:
        for segment in perc_identity_matrix:
            status += segment[0] > segment[1]
    
        if status != 0 and status != len(perc_identity_matrix):
            chimera = True
        
    return chimera


def get_unique(ids):
    """Not used.
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    """Not used.
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Get list of non-overlapping chunks
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers
    """
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Identity percentage between 2 sequences
    """
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Main function of the program, yields otu that are no chimera
    """
    matrix = os.path.abspath(os.path.join("MATCH"))
    otu_list = abundance_greedy_clustering(amplicon_file,
                                           minseqlen,
                                           mincount,
                                           chunk_size,
                                           kmer_size)
    kmer_dict = {}
    # first 2 are considered non-chimeric by default
    for sequence in otu_list:
        for i in range (2):
            get_unique_kmer(kmer_dict,
                            sequence[0],
                            otu_list.index(sequence),
                            kmer_size)
            
        mates = search_mates(kmer_dict, sequence[0], kmer_size)
        
        
        if len(mates) > 1:
            chunks_seq = get_chunks(sequence[0], chunk_size)
            chunks_parent1 = get_chunks(otu_list[(mates[0])][0], chunk_size)
            chunks_parent2 = get_chunks(otu_list[(mates[1])][0], chunk_size)
        
            # similarity matrix
            mates_matrix = []
            for i in range(len(chunks_seq)):
                id_1 = get_identity(nw.global_align(chunks_seq[i],
                                                    chunks_parent1[i],
                                                    gap_open=-1,
                                                    gap_extend=-1,
                                                    matrix=matrix))
                id_2 = get_identity(nw.global_align(chunks_seq[i],
                                                    chunks_parent2[i],
                                                    gap_open=-1,
                                                    gap_extend=-1,
                                                    matrix=matrix))
                mates_matrix.append([id_1, id_2])
            if detect_chimera(mates_matrix):
                otu_list.remove(sequence)
            else:
                get_unique_kmer(kmer_dict,
                                sequence[0],
                                otu_list.index(sequence),
                                kmer_size)
    for otu in otu_list:
        yield otu


def fill(text, width=80):
    """Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """Write sequences in a fasta file
    """

    with open(output_file, 'w+') as o_file:
        i = 1
        for otu in OTU_list:
            o_file.write('>OTU_%i occurrence:%i\n%s\n' % (i,
                                                      otu[1],
                                                      fill(otu[0])))
            i += 1

#==============================================================
# Main program
#==============================================================
def main():
    """main fonction
    """
    # Get arguments
    args = get_arguments()
    
    otu_list = []
    for otu in chimera_removal(args.amplicon_file,
                               args.minseqlen,
                               args.mincount,
                               args.chunk_size,
                               args.kmer_size):
        otu_list.append(otu)

    write_OTU(otu_list, args.output_file)

if __name__ == '__main__':
    main()
