#!/usr/bin/env python

#############################################################
#
# utils.py
#
# Author : Miravet-Verde, Samuel
# Last updated : 11/15/2017
#
# [2017] - Centre de Regulació Genòmica (CRG) - All Rights Reserved
#############################################################

import os.path
import collections
import numpy as np
from Bio import SeqIO

def splitn_str(your_string, n):
    """ Given a string, returns a list with that string splitted each n characters """

    return [your_string[i:i+n] for i in range(0, len(your_string), n)]


def indexes(lista, values):
    """
    Given a list and its values return a list with the indexes of the values
    in lista
    """

    indexes = []
    for element in values:
        indexes.append(lista.index(element))
    return indexes

def list_generator(filename, index=0):
    """
    Given a file, returns a set with all the values from the column[index]
    """
    results = []
    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            results.append(line[index])
    return results

def list_NA_generator(filename, index=0):
    """
    Given a file, returns a set with all the values from the column[index]
    """
    results = []
    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            try:
                if len(line[index]) > 0:
                    results.append(line[0])
            except:
                pass
    return results

def set_generator(filename, index):
    """
    Given a file, returns a set with all the values from the column[index]
    """
    results = set()
    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            results.add(line[index])
    return results


def dic_generator(filename, key_index, value_index=None, header=False):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    results = {}
    with open(filename, 'r') as fi:
        for line in fi:
            if header:
                header=False
            else:
                line = line.strip().split()
                if not value_index:
                    results[int(line[key_index])] = float(line[key_index+1])
                else:
                    results[int(line[key_index])] = float(line[value_index])
    return results

def new_dic_generator(filename, key_index, value_index):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    results = {}
    with open(filename, 'r') as fi:
        results = {int(k):float(v) for k, v in [l.split()[0:2] for l in fi.readlines()]}
    return results

def str_dic_generator(filename, key_index, value_index=False, header=False, split_by=None):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    if header==True:
        header = 1

    results = {}
    with open(filename, 'r') as fi:
        for line in fi:
            if header != 0:
                header-=1
            else:
                if split_by:
                    line = line.strip().split(split_by)
                else:
                    line = line.strip().split()
                if value_index or value_index == 0:
                    results[line[key_index]] = line[value_index]
                else:
                    results[line[key_index]] = line[key_index+1]
    return results

def genes_coordinates(caps = True):
    """
    Given a three colums file returns a dictionary where the first column is the key and the other are the values
    in a list
    """

    results_dic = {}

    with open('/home/smiravet/crg/transpgrowth/datasets/essentials/gene_coordinates.txt', 'r') as fi:
        for line in fi:
            line = line.strip().split()
            if caps:
                results_dic[line[0].upper()] = [line[1], line[2]]
            else:
                results_dic[line[0]] = [line[1], line[2]]

    return results_dic

def dict2file(dictionary, filename):
    """
    Writes a file where the first column is the key and the second the values
    """

    directorytosave = '/home/smiravet/crg/transpgrowth/datasets/'

    fo = open(directorytosave+filename+'.txt', 'w')

    od = collections.OrderedDict(sorted(dictionary.items()))
    for k, v in od.iteritems():
        fo.write(str(k)+'\t'+str(v)+'\n')

    fo.close()

def load_multifasta(inFile):
    """ Return a dictionary wit the sequences from a multifasta file """
    your_sequences = {}
    handle = open(inFile, 'rU')
    for record in SeqIO.parse(handle, "fasta"):
        your_sequences[record.id]=str(record.seq)
    handle.close()
    return your_sequences


def load_multifasta_info(inFile):
    """ Return a dictionary wit the sequences from a multifasta file """
    your_sequences = {}
    handle = open(inFile, 'rU')
    for record in SeqIO.parse(handle, "fasta"):
        your_sequences[record.id+'//'+record.description]=str(record.seq)
    handle.close()
    return your_sequences


def load_genome(genome):
    # Determine the file type:
    if genome.endswith('gb') or genome.endswith('gbk') or genome.endswith('genbank'):
        tipo = 'genbank'
    else:
        tipo = 'fasta'
    handle = open(genome, 'rU')
    for record in SeqIO.parse(handle, tipo):
        return str(record.seq)
    handle.close()


def load_annotation(inFile):
    annotation = {}
    with open(inFile) as fi:
        for line in fi:
            line = line.strip().split()
            ide  = line[0]
            if len(line[1:]) == 2:
                st, en = sorted([int(x) for x in line[1:]])
                l      = [st, en]
            else:
                st, en = sorted([int(x) for x in line[1:] if x not in ['+', '-']])
                strand = [str(x) for x in line[1:] if x in ['+', '-']]
                l      = [st, en] + strand
            annotation[ide] = l
    return annotation


def strand_load_annotation(inFile):
    annotation = {}
    with open(inFile) as fi:
        for line in fi:
            line = line.strip().split()
            ide  = line[0]
            if len(line[1:]) == 2:
                st, en = [int(x) for x in line[1:]]
                l      = [st, en]
            else:
                st, en = [int(x) for x in line[1:] if x not in ['+', '-']]
                strand = [str(x) for x in line[1:] if x in ['+', '-']]
                l      = [st, en] + strand
            annotation[ide] = l
    return annotation


def gb2annotation(inFile):
    annotation = {}
    c = 1
    for rec in SeqIO.parse(inFile, "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    strand = '-' if feature.location.strand==-1 else '+'
                    start  = int(feature.location.start)
                    stop   = int(feature.location.end)
                    try:
                        genename = feature.qualifiers["locus_tag"][0]
                    except:
                        try:
                            genename = feature.qualifiers["gene"][0]
                        except:
                            genename = 'unknown'+str(c)
                            c += 1
                    annotation[genename] = [start, stop, strand]
    return annotation

def lists2dict(listA, listB):
    """ Given two lists of the same length, merge them in one dictionary """
    return dict(zip(listA, listB))


def remove_column(array, index):
    """ Remove the index column from a numpy array, index can be a list"""

    return np.delete(array, np.s_[index], axis=1)

##### FOR SEQUENCES

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    return ''.join([complement[k] if k in complement else 'N' for k in seq][::-1])

# [2017] - Centre de Regulació Genòmica (CRG) - All Rights Reserved
