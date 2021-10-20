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
    """
    prend deux arguments correspondant au fichier fasta.gz et à la longueur
    minimale des séquences et retourne un générateur de séquences de longueur l
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
    Prend trois arguments correspondant au fichier fasta.gz,  la longueur minimale
    des séquences et leur comptage minimum. Elle fait appel au générateur fourni 
    par read_fasta et retourne un générateur des séquences uniques ayant une 
    occurrence >=mincount
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
    """
    utilise la fonction précédente pour générer les séquences > minseqlen et dont
    l'occurence > mincount et les ajoutent à une liste (otu_list) si elles ont moins
    de 97% de similarité avec les séquences déjà présentent dans la liste.

    renvoie une liste d'otu

    les arguments chunk_size et kmer_size sont inutiles mais doivent rester
    pour pouvoir passer le test
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
    """
    met à jour un dictionnaire de kmers à partir des kmers générés par la
    séquence, les kmers générés correspondent aux clés et la valeur associé
    correspond à une liste d'id_seq pour pouvoir identifier dans quelles séquences
    ces kmers ont été observés.

    renvoie le dictionnaire mis à jour
    """

    for kmer in cut_kmer(sequence, kmer_size):
        kmer_dict.setdefault(kmer, [])
        kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates (kmer_dict, sequence, kmer_size):
    """
    génère des kmers à partir d'une séquence, recherche les kmers identiques dans
    le dictionnaire de kmers et répertorie les id_seq correspondant pour ensuite
    cherche les deux id_seq qui ressortent le plus souvent et qui sont les 2 séquences
    parentes

    renvoie les id des 2 séquences parentes

    """
    common = []
    
    common.append(kmer_dict[kmer] for kmer in cut_kmer(sequence, kmer_size) \
                                                if kmer in kmer_dict.keys())
    flatten_list = [seq for sublist in common for item in sublist \
                                                    for seq in item]

    return [seq[0] for seq in Counter(flatten_list).most_common(2)]


def detect_chimera(perc_identity_matrix):
    """
    prend une matrice donnant par segment le taux d'identité entre la séquence
    candidate et deux séquences parentes et retourne un booléan indiquant si la
    séquence candidate est une chimère (True) ou ne l'est pas (False)
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
    """
    fonction présente de base, non utilisée
    """
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    """
    fonction présente de base, non utilisée
    """
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """
    prend une séquence et une longueur de segment chunk_size et retourne une liste
    de sous-séquences de longueur chunk_size non chevauchantes
    """
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
    """Prend une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    exploite toutes les fonctions développées avant pour générer les séquences,
    ces séquences sont ajoutées dans la liste des otu si elles présentent moins de
    97% de similarité (abundance_greedy_clustering). Puis on teste dans cette liste
    d'otu le caractère chimérique (si oui -> elles sont retirées de la liste des otu)
    (sinon -> elles sont conservés dans la liste des otu et utilisées pour compléter
    le dictionnaire des kmers)

    renvoie un générateur d'otu (serait plus simple de renvoyer une liste d'otu
    directement mais le test utilise un next())
    """
    matrix = os.path.abspath(os.path.join("MATCH"))
    otu_list = abundance_greedy_clustering(amplicon_file,
                                           minseqlen,
                                           mincount,
                                           chunk_size, 
                                           kmer_size)
    kmer_dict = {}
    
    #les 2 premiers otu sont considérés comme non chimériques:
    for sequence in otu_list:
        for i in range (2):
            get_unique_kmer(kmer_dict,
                            sequence[0],
                            otu_list.index(sequence),
                            kmer_size)
            
        #pour les autres séquences on cherche dans la liste des otu des séquences similaires(grâce aux kmers) puis on compare les chunks:
        mates = search_mates(kmer_dict, sequence[0], kmer_size)
        
        
        if len(mates) > 1:
            chunks_seq = get_chunks(sequence[0], chunk_size)
            chunks_parent1 = get_chunks(otu_list[(mates[0])][0], chunk_size)
            chunks_parent2 = get_chunks(otu_list[(mates[1])][0], chunk_size)
        
            #on crée la matrice de similarité:
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
    
            #on regarde si la séquence est une chimère:
            #si oui:
            if detect_chimera(mates_matrix):
                otu_list.remove(sequence)
            
            #sinon:
            else:
                get_unique_kmer(kmer_dict,
                                sequence[0],
                                otu_list.index(sequence),
                                kmer_size)

    # on renvoie les éléments de otu_list 1 par 1 puisque le test utilise un next()
    for otu in otu_list:
        yield otu


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """
    prend une liste d'otu et le chemin vers un fichier de sortie et affiche les otu
    au format demandé
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
    """
    Main program function

    récupère les arguments

    utilise la fonction chimera_removal qui utilise toutes les fonctions pour générer les
    otu qu'on restocke dans une liste

    on utilise cette liste d'otu pour écrire le fichier
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
