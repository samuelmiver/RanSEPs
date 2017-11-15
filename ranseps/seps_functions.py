#!/usr/bin/env python

#############################################################
#
# seps_functions.py
#
# Author : Miravet-Verde, Samuel
# Written : 03/06/2017
# Last updated : 11/15/2017
#
#############################################################

#####################
#      IMPORTS      #
#####################

import os
import sys
import glob
import math
import random
import pickle
import inspect
import os.path
import itertools
import collections
import utils as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from propy.PyPro import GetProDes
# To compute specifici features
from propy.QuasiSequenceOrder import GetSequenceOrderCouplingNumberSW as SOSW
from propy.CTD import CalculateDistributionHydrophobicity as DHF
from propy.CTD import CalculateDistributionSecondaryStr as DSS

# Others
from scipy import interp
from operator import mul
from collections import Counter
from numpy.random import choice
from joblib import Parallel, delayed                # To parallelize
from numpy import genfromtxt, savetxt
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve, average_precision_score

# Bio handling
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

# For aligning stuff
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

# For computing RNA structures
from RNA import fold

# Graphing
# To autofit the figure sizes
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 40})

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
          'figure.figsize': (15, 15),
          'axes.labelsize': 'large',
          'axes.titlesize':'large',
          'xtick.labelsize':'large',
          'ytick.labelsize':'large'}
pylab.rcParams.update(params)

#####################
# USEFUL  VARIABLES #
#####################

codes = {'magalactiae':'NC_013948.1'   ,  'mbovis':'NC_014760.1'    ,
         'mNMbovis':'NZ_CP011348.1'    , 'mcapricolum':'NC_007633.1',
         'mgenitalium':'NC_000908.2'   , 'mhigh':'NC_017502.1'      ,
         'mhyopneumoniae':'NC_006360.1', 'mlow':'NC_004829.2'       ,
         'mmycoides':'NZ_CP001668.1'   , 'mpneumoniae':'NC_000912.1'}

#####################
#  BASIC FUNCTIONS  #
#####################

def zstandarization(value, distribution=None):
    """
    Apply a z standarization to the value.
    Value can be another distribution.
    """
    if distribution:
        return (np.array(value)-np.mean(np.array(distribution)))/np.std(np.array(distribution))
    else:
        return (np.array(value)-np.mean(value))/np.std(value)

def minmaxstandarization(value, distribution=None):
    """
    Apply a min max standarization to the value.
    Value can be another distribution
    """
    if distribution:
        return (np.array(value)-min(distribution))/(max(distribution)-min(distribution))
    else:
        return (np.array(value)-min(value))/(max(value)-min(value))

def quantile_grouping(distribution):
    """
    Apply a min max standarization to the value.
    Value can be another distribution
    """
    new = []
    a = np.percentile(distribution, 20)
    b = np.percentile(distribution, 40)
    c = np.percentile(distribution, 60)
    d = np.percentile(distribution, 80)

    for i in distribution:
        if i >= d:
            new.append(1.0)
        elif i >= c:
            new.append(0.75)
        elif i >= b:
            new.append(0.5)
        elif i >= a:
            new.append(0.25)
        else:
            new.append(0.0)

    return new

def central_limit(distribution, i, n):
    """
    Apply central limit theorem to an array.
    i is the number of times you sample and n the size
    """
    sample_means = []
    for j in range(0, i):
        sample_means.append(np.random.choice(distribution, n, replace=False))
    return np.array(sample_means)


def tryptic_peptides(seq, minimum=5, partial=False):
    """
    Compute the number of tryptic number of peptides with size greater than <minimum>
    This process is performed spliting the protein in its lysines (K) or arginines (R)).

    If partial is False peptides are considered to be fully digested. If true, partial
    digestions are returned.
    """

    seq = seq.replace('K', 'R')
    rs = seq.count('R')

    if partial:
        if rs == 0:
            total = 1
        elif rs == 1:
            total = 3
        else:
            total = (rs*(rs+1))/2
    else:
        total = len([x for x in seq.split('R') if len(x) > minimum])
    return total


def SEQ2TP(seq, minimum=5):
    """
    Compute the number of tryptic number of peptides with size greater than <minimum>
    This process is performed spliting the protein in its lysines (K) or arginines (R)).
    """
    seq = seq.replace('K', 'R')
    rs = seq.count('R')
    total = set([x for x in seq.split('R') if len(x) > minimum])
    return total


def check_nt(seq, length=True):
    alphabet = ['A','C','G','T']
    counts   = sum([0 if c in alphabet else 1 for c in seq])
    sobrante = len(seq)%3
    if not length:
        # In order to compute motives
        sobrante=0
    if counts+sobrante != 0:
        return False
    else:
        return True

def check_aa(seq):
    alphabet = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    counts   = [0 if c in alphabet else 1 for c in seq]
    if sum(counts) > 0 or any(x in seq for x in ['J','U','Z','B','O','X']):
        return False
    else:
        return True

#####################
# DATABASE HANDLING #
#####################

def match_gallisepticum():
    """ Match annotation between Luca identifiers and yours """

    with open('./DBs/gallisepticum_correspondance.txt', 'w') as fo:
        for yours, luca in [['mhigh', 'NC_017502.1']]:
            print yours
            c = 0
            aay = {k:v for k, v in u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/filtered_by_size_orfinder/'+yours+'_small_aa.fa').iteritems() if len(v)>=19}
            ann = {k:v[:-1] for k, v in u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+yours+'_annotation.txt').iteritems() if k in aay.keys()}
            aal = {k:v for k, v in u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/luca_dbs/'+luca+'.txt').iteritems()}

            laa = {}
            for k, v in aal.iteritems():
                if v in laa:
                    laa[v].append(k)
                else:
                    laa[v] = [k]

            # load annotation luca:
            annl = {}
            with open('/home/smiravet/crg/dbs/smprots_DB/luca_dbs/'+luca+'.txt', 'r') as fi:
                for line in fi:
                    if line.startswith('>'):
                        ide  = line.split(' [')[0][1:]
                        sten = line.split('[')[1].split(']')[0].split(' - ')
                        annl[ide] = [int(se) for se in sten]

            # Match
            for k, v in aay.iteritems():
                if len(laa[v]) == 1:
                    # one match
                    c+=1
                    fo.write(k+'\t'+laa[v][0]+'\n')
                else:
                    for i in laa[v]:
                        if len(set(ann[k]).intersection(set(annl[i]))) == 1:
                            c+=1
                            fo.write(k+'\t'+i+'\n')
            assert c==len(aay)==len(aal)==len(ann)==len(annl)

def match_luca_yours():
    """ Match annotation between Luca identifiers and yours """

    with open('./DBs/correspondance.txt', 'w') as fo:
        for yours, luca in codes.iteritems():
            print yours
            c = 0
            aay = {k:v for k, v in u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/filtered_by_size_orfinder/'+yours+'_small_aa.fa').iteritems()}
            ann = {k:v[:-1] for k, v in u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+yours+'_annotation.txt').iteritems() if k in aay.keys()}
            aal = {k:v for k, v in u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/luca_dbs/'+luca+'.txt').iteritems()}

            laa = {}
            for k, v in aal.iteritems():
                if v in laa:
                    laa[v].append(k)
                else:
                    laa[v] = [k]

            # load annotation luca:
            annl = {}
            with open('/home/smiravet/crg/dbs/smprots_DB/luca_dbs/'+luca+'.txt', 'r') as fi:
                for line in fi:
                    if line.startswith('>'):
                        ide  = line.split(' [')[0][1:]
                        sten = line.split('[')[1].split(']')[0].split(' - ')
                        annl[ide] = [int(se) for se in sten]

            # Match
            for k, v in aay.iteritems():
                if len(laa[v]) == 1:
                    # one match
                    c+=1
                    fo.write(k+'\t'+laa[v][0]+'\n')
                else:
                    for i in laa[v]:
                        if len(set(ann[k]).intersection(set(annl[i]))) == 1:
                            c+=1
                            fo.write(k+'\t'+i+'\n')
            assert c==len(aay)==len(aal)==len(ann)==len(annl)


def load_orthology(species_code, filename, min_conservation=0, substract=0):
    """ Return the conserved identifiers """

    xl = pd.ExcelFile(filename)
    df = xl.parse("orthology")

    orth_set = set()

    c = 0
    for ide in df[species_code]:
        ide = str(ide)
        if ide != '*':
            if min_conservation > 0:
                m = df['# Species'][c] - substract
            else:
                m = 0
            if m >= min_conservation:
                ide = ide.split(',')
                if len(ide) != 1:
                    for subide in ide:
                        orth_set.add(subide)
                else:
                    orth_set.add(ide[0])
        c += 1
    return orth_set


def load_MS():

    results = []
    xl = pd.ExcelFile('/home/smiravet/crg/antprotect/RF_final/DBs/MS_seps.xlsx')
    df = xl.parse("this")

    for index, row in df.iterrows():
        if row[0] == 1:
            for cons_pept in row[4:]:
                pre = []
                cons_pept = cons_pept.split(',')
                for sep in cons_pept:
                    if sep != '*':
                        results.append(str(sep))
    return results


def count_code(code, lista):
    return float(len([i for i in lista if i.startswith(code)]))


def load_correspondance():

    try:
        print 'loading correspondance'
        correspondance = u.str_dic_generator(filename='/home/smiravet/crg/antprotect/RF_final/DBs/correspondance.txt', key_index=1, value_index=0)
        print 'correspondance loaded'
    except:
        # Match luca dataset
        print 'generating correspondance'
        match_luca_yours()
        correspondance = u.str_dic_generator(filename='/home/smiravet/crg/antprotect/RF_final/DBs/correspondance.txt', key_index=1, value_index=0)
        print 'correspondance generated'

    return correspondance


def organism_info(organism, genome_size, size=19):

    print 'loading '+organism+' info'

    # Load reference
    if organism == 'mpneumoniae':
        NCBI1 = u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/id_genes/'+organism+'_pairs.txt', 0, 2, header=False, split_by='\t')
        your_NCBI = {ide:'NA_hypothetical' if ide not in NCBI1 else NCBI1[ide] for ide in u.list_generator('/home/smiravet/crg/dbs/smprots_DB/id_genes/'+organism+'_pairs2.txt', 0)}
    else:
        your_NCBI = u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/id_genes/'+organism+'_pairs.txt', 0, 2, header=False, split_by='\t')
    print 'reference annotation loaded'

    # Load sequences
    your_nt_seqs = u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/normal_orfinder/'+organism+'_small_nt.fa')
    prev_aa_seqs = u.load_multifasta('/home/smiravet/crg/dbs/smprots_DB/normal_orfinder/'+organism+'_small_aa.fa')
    your_aa_seqs, your_lengths = {}, {}
    for k, v in prev_aa_seqs.iteritems():
        length = len(v)
        if length >= size or k in your_NCBI:
            your_aa_seqs[k] = v
            your_lengths[k] = length
    assert len(set(your_NCBI).difference(set(your_aa_seqs)))==0
    print 'sequences loaded'

    # Load annotation & correspondance
    annotation =  u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+organism+'_annotation.txt')
    your_annotation = {}
    for k in your_aa_seqs.keys():
        your_annotation[k] = annotation[k]
    print 'full annotation loaded'

    # Frames and length
    your_frames, your_contra = add_frames(your_annotation, your_NCBI, genome_size)
    print 'frames loaded'

    return your_nt_seqs, your_aa_seqs, your_NCBI, your_annotation, your_lengths, your_frames, your_contra


def load_DBs(organisms, size=19, random_sample=False):
    """ Given a lists of organisms generate a pickle object with all the variables required """
    nt_seqs, aa_seqs, NCBI_all, anno_all, lengths, frames, contra_frames, rd_seqs, GCs, genome_lengths = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    conserved = []

    correspondance = load_correspondance()

    not_general = ['bsubtilis', 'ecoli']

    # Iterate organism loadingthe infor
    for organism in organisms:
        genome = u.load_genome_DB(organism)
        GCs[organism] = GC(genome)
        genome_lengths[organism] = len(genome)

        # Load conservation
        if organism not in not_general:
            pre_conserved  = load_orthology(codes[organism], '/home/smiravet/crg/antprotect/RF_final/DBs/orthology.xlsx')
            conserved     += [correspondance[ide] for ide in pre_conserved]
            print 'Conservation loaded'


        # Update all the information
        your_nt_seqs, your_aa_seqs, your_NCBI, your_annotation, your_lengths, your_frames, your_contra = organism_info(organism, genome_lengths[organism], size)

        for k in your_aa_seqs.keys():
            nt_seqs[k]       = your_nt_seqs[k]
            aa_seqs[k]       = your_aa_seqs[k]
            anno_all[k]      = your_annotation[k]
            lengths[k]       = your_lengths[k]
            frames[k]        = your_frames[k]
            if k not in your_NCBI:
                contra_frames[k] = your_contra[k]

        for k, v in your_NCBI.iteritems():
            NCBI_all[k] = v

        print 'DBs updated'

        if random_sample:
            your_rd_seqs.update(randomizame(min_len=57, max_len=1500, n=200, genome=genome, annotation=annotation, add=organism[:3], inter=False))
            print 'rd seqs loaded'

    # Saving the objects:
    with open('/home/smiravet/crg/antprotect/RF_final/DBs/objs.pickle', 'w') as f:
        pickle.dump([nt_seqs, aa_seqs, NCBI_all, anno_all, lengths, frames, contra_frames, rd_seqs, GCs, genome_lengths, conserved, correspondance], f)

    print 'object saved'


############ RANDOMIZING FUNCTIONS

def randomizame(min_len, max_len, n, genome, annotation, add, inter):
    """
    Given a minimum len, a max len for peptides, returns <n> random sequences
    50% sampled from intergenic regions of <genome> and 50% preserving the GC(<genome>)

    Genome has to be loaded with utils functions!

    Add is used for the identifiers and inter is to select or not intergenic regions
    """

    # Define strictly intergenic regions:
    GC_content = GC(genome)
    genome_positions = set(range(1, len(genome)+1))

    annotations_pos = set()         # One for each strand
    annotations_neg = set()

    for k, v in annotation.iteritems():
        if v[2] == '+':
            st, en = sorted([int(x) for x in v[0:2]])
            annotations_pos.update(range(st, en+1))
        else:
            st, en = sorted([int(x) for x in v[0:2]])
            annotations_neg.update(range(st, en+1))

    positions_pos = genome_positions.difference(annotations_pos)
    positions_neg = genome_positions.difference(annotations_neg)

    # Load the random sequences
    identifiers, random_seqs = [], []
    maxi = 5
    c    = 0
    yy   = 1    # GC iterator

    if inter:
        yy = 2
        # From positive intergenic regions
        for i in range(n/4):
            lene = random.randrange(min_len, max_len, 3) # Random len
            sete = set(random.sample(positions_pos, lene))

            # Return random seq:
            seqi = ''.join([genome[i-1] for i in sete])
            random_seqs.append(seqi)
            identifiers.append('in_plus'+add+str(c+1).zfill(maxi))
            c += 1

        # From negative intergenic regions
        for i in range(n/4):
            lene = random.randrange(min_len, max_len, 3)
            sete = set(random.sample(positions_neg, lene))
            seqi = ''.join([genome[i-1] for i in sete])
            random_seqs.append(u.reverse_complement(seqi))
            identifiers.append('in_minus'+add+str(c+1).zfill(maxi))
            c += 1

    # Randomizing from the GC content
    GC_content = GC_content/100.0
    for i in range(n/yy):
        lene = random.randrange(min_len, max_len, 3)
        seqi = ''.join(choice(['A', 'C', 'G', 'T'], lene, p=[(1-GC_content)/2, GC_content/2, GC_content/2, (1-GC_content)/2]))
        random_seqs.append(seqi)
        identifiers.append('gc'+add+str(c+1).zfill(maxi))
        c += 1

    random_seqs = {identifiers[i]:random_seqs[i] for i in range(0, len(random_seqs))}
    return random_seqs

############ FRAMES FUNCTIONS

def return_frame_plus(stN, st):
    """ Return frame of st in reference to stN """

    return (st - stN)%3

def percentage_overlap(a, b):
    """ a is the NCBI st and en, b for the SEP!"""

    dec = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)/float(b[1] - b[0] + 1)
    return dec*100.0


def define_overlaps_plus(SEPs, NCBI_anno):
    """
    Inner function to find overlaps
    """

    results = {}
    for ide, features in SEPs.iteritems():
        st, en, strand        = features
        for ideN, featuresN in NCBI_anno.iteritems():
            stN, enN, strandN = featuresN
            results[ideN] = ['NCBI', 0.0, 'NA']
            if stN <= st < en <= enN or st < stN <= en <= enN or stN <= st <= enN < en:
                if strand == strandN:
                    if strand == '+':
                        frame = return_frame_plus(stN, st)
                    else:
                        frame = return_frame_plus(st, stN)
                else:
                    frame = 'contra'
                perce = percentage_overlap([stN, enN], [st, en])
                results[ide] = [frame, perce, ideN]
            else:
                pass
    # Add NO Overlapping ones
    for ide in SEPs.keys():
        if ide not in results:
            results[ide] = ['NO'  , 0.0, 'NA']
    return results


def add_frames(annotation_dictionary, NCBI_ids, genome_size):
    """
    Given a dictionary: {geneid:[st, en, strand]...} and a list of geneid that are from NCBI
    Transform the set of proteins in a dictionary with the structure
    {geneid:[st, en, strand, loc, %, gene]...} where loc is NCBI, 0, 1, 2, NO (NO overlap), PO (partial overlap) and the percentage of
    overlap (0 in NCBI and NO) depending on the fram its located and the gene

    LIMITATION: when an annotation overlaps with more than one NCBI gene
    """
    # Define the NCBI annotation
    NCBI_plus     = {k:v for k, v in annotation_dictionary.iteritems() if v[2] == '+' and k in NCBI_ids}
    NCBI_minus    = {k:sorted(v[:2])+['-'] for k, v in annotation_dictionary.iteritems() if v[2] == '-' and k in NCBI_ids}

    # Plus strand is direct
    plus_dic      = {k:v for k, v in annotation_dictionary.iteritems() if v[2] == '+' and k not in NCBI_ids}
    results_plus  = define_overlaps_plus(plus_dic, NCBI_plus)
    contras_plus  = define_overlaps_plus(plus_dic, NCBI_minus)

    # Minus strand require a small transformation
    minus_dic     = {k:sorted(v[:2])+['-'] for k, v in annotation_dictionary.iteritems() if v[2] == '-' and k not in NCBI_ids}
    results_minus = define_overlaps_plus(minus_dic, NCBI_minus)
    contras_minus = define_overlaps_plus(minus_dic, NCBI_plus)


    # Add annotation, generate and return the results
    results, results_contra  = {}, {}
    for k, v in annotation_dictionary.iteritems():
        if v[-1] == '+':
            results[k] = v + results_plus[k]
            if k not in NCBI_plus.keys()+NCBI_minus.keys():
                results_contra[k] = v + contras_plus[k]
        else:
            results[k] = v + results_minus[k]
            if k not in NCBI_plus.keys()+NCBI_minus.keys():
                results_contra[k] = v + contras_minus[k]

    return results, results_contra


def check_frames():
    # Load databases
    with open('./DBs/objs.pickle') as f:
        nt_seqs, aa_seqs, NCBI_all, anno_all, lengths, frames, contra_frames, rd_seqs, GCs, genome_lengths, conserved, correspondance = pickle.load(f)

    NCBI = [x for x in NCBI_all.keys() if x.startswith('mpneu')]

    ncbipos, ncbineg = [], []
    for i in NCBI:
        st, en = sorted(anno_all[i][:2])
        if anno_all[i][-1] == '+':
            ncbipos += range(st, en+1)
        else:
            ncbineg += range(st, en+1)

    ncbipos = set(ncbipos)
    ncbineg = set(ncbineg)

    for k, v in frames.iteritems():
        if k.startswith('mpneu'):
            if v[3] == 'NO':
                st, en = sorted(anno_all[k][:2])
                if v[2] == '-':
                    s = set(range(st, en+1)).intersection(ncbineg)
                else:
                    s = set(range(st, en+1)).intersection(ncbipos)
                if len(s) > 0:
                    print k


def motif_probability_GC(motif, GC):
    """
    Given a motif <string> including AGCTN and a GC content (float, 100 base, e.g. 40.0)

    Computes the probability of finding that motif in the genome based on the combinatorial
    chance of each base appearing independently
    """

    probabilities = {'A':((100-GC)/100.0)/2.0,
                     'T':((100-GC)/100.0)/2.0,
                     'G':(GC/100.0)/2.0,
                     'C':(GC/100.0)/2.0,
                     'N':1.0}

    return reduce(mul, [probabilities[n] for n in motif], 1)


def valid_DNA(sequence):
    """ Check all the character in sequence are A C G or T """
    valid_dna = 'ACGT'
    sequence = sequence.upper()
    return(all(i in valid_dna for i in sequence))


def stack_energy(sequence):
    """
    Given a sequence string

    Returns stacking energy at 37C (AG37)
    """
    # First of all define the dictionary of values
    energies = {'AA':-1.00, 'AC':-1.44, 'AG':-1.28, 'AT':-0.88,
                'CA':-1.45, 'CC':-1.84, 'CG':-2.17, 'CT':-1.28,
                'GA':-1.30, 'GC':-2.24, 'GG':-1.84, 'GT':-1.44,
                'TA':-0.58, 'TC':-1.30, 'TG':-1.45, 'TT':-1.00}
    # Transform the string in a list of duplets
    duplets = [sequence[i:i+2] for i in range(0, len(sequence)-1, 1)]
    # Transform the duplets to its energies:
    stack_energies = [energies[duplet] if duplet in energies else 0.0 for duplet in duplets] # For cases like ecoli where you have R and N in the sequence
    # Sum and normalize by the length of the sequence
    if len(sequence) > 0:
        return (sum(stack_energies))/len(sequence)
    else:
        return 0.0

# RBS
def hasRBS(seq):
    """ Return True if sequence has an RBS """
    RBS_seqs = ['GGA'   , 'GAG'   , 'AGG'   ,
                'AGGA'  , 'GGAG'  , 'GAGG'  ,
                'AGGAG' , 'GGAGG' ,
                'AGAAGG', 'AGCAGG', 'AGGAGG', 'AGTAGG',
                'AGGCGG', 'AGGGGG', 'AGGTGG']
    return 1.0 if any(RBS_seq in seq for RBS_seq in RBS_seqs) else 0.0


def generate_DB_stack_energy(organism, annotation=None, extension=None, write=False):
    """
    Generates a DB with the values for stacking energy,
    If annotation is provided that dictionary is used instead of the general
    annotation of the organism is not used
    """
    # Load genome and annotation if required
    genome = u.load_genome_DB(organism)
    lg     = len(genome)
    genome += genome # circular
    if not annotation:
        annotation = u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+organism+'_annotation.txt')

    # Compute and write
    if extension:
        extension += '_' # To make nicer filenames
    else:
        extension  = ''

    # If write is activated write in a file
    if write:
        fo1 = open('/home/smiravet/crg/dbs/smprots_DB/stack_E/'+organism+'_'+extension+'stackE_RBS.txt', 'w')
        fo2 = open('/home/smiravet/crg/dbs/smprots_DB/stack_E/'+organism+'_'+extension+'stackE_m10p20.txt', 'w')

    stack_energies_dic = {}
    for k, v in annotation.iteritems():
        if v[-1] == '+':
            st = v[0]
            if st > 16:
                seq_RBS    = genome[st-16:st-1]
                seq_m10p20 = genome[st-11:st+19]
            else:
                seq_RBS    = genome[st+lg-16:st+lg-1]
                seq_m10p20 = genome[st+lg-11:st+lg+19]
        else:
            st = max(v[0:2])
            if st < lg-16:
                seq_RBS    = genome[st:st+15]
                seq_m10p20 = genome[st-23:st+10]
            else:
                seq_RBS    = genome[st+lg:st+lg+15]
                seq_m10p20 = genome[st+lg-23:st+lg+10]
        # Add the result
        sE1 = stack_energy(seq_RBS)
        RBS = hasRBS(seq_RBS)
        sE2 = stack_energy(seq_m10p20)
        stack_energies_dic[k] = [sE1, sE2, RBS]

        if write:
            fo1.write(k+'\t'+str(sE1)+'\n')
            fo2.write(k+'\t'+str(sE2)+'\n')
    if write:
        fo1.close()
        fo2.close()

    return stack_energies_dic

#####################
# FEATURE FUNCTIONS #
#####################

def run_m10p20(organism, training_ides, testing_ides, annotation=False, random=False, rand_dict=None, verbosity=False):
    """
    Look for motives in the -10 +20 region of the SEP
    If random == True --> testing ides has to be the sequences!
    """
    # 1.Load the coordinates annotations
    if not annotation:
        annotation = {}
        with open('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+organism+'_annotation.txt', "rU") as fi:
            for line in fi:
                line   = line.strip().split()
                ide    = line[0]
                strand = line[-1]
                st, en = sorted([int(x) for x in line[1:3]])
                annotation[ide] = [st, en, strand]

    # 2. Load the genome
    if os.path.isfile('/home/smiravet/crg/antprotect/smprots_DB/genomes/'+organism+'.fa'):
        filehandle = '.fa'
    elif os.path.isfile('/home/smiravet/crg/antprotect/smprots_DB/genomes/'+organism+'.fasta'):
        filehandle = '.fasta'
    else:
        filehandle = '.gb'
    genome = u.load_genome("/home/smiravet/crg/antprotect/smprots_DB/genomes/"+organism+filehandle)
    GC_content = GC(genome)/100.0

    # 3. Only take those sequences with a distance >50 from the stop
    # of the previous sequence
    def dist2closest(myList, myNumber):
        myList = [x for x in myList if x != myNumber]
        return abs(myNumber - (min(myList, key=lambda x:abs(x-myNumber))))

    def subset_ides(annotation, strand, training):
        ides, starts, ends = [], [], []
        your_set = []
        for ide, coords in annotation.iteritems():
            if ide in training:
                if coords[2] == strand:
                    ides.append(ide)
                    if strand == '+':
                        starts.append(coords[0])
                        ends.append(coords[1])
                    else:
                        starts.append(coords[1])
                        ends.append(coords[0])
        for i in range(len(ides)):
            if dist2closest(ends, starts[i]) > 50:
                your_set.append(ides[i])

        # assert len(your_set) > 25
        return your_set

    your_set = []
    for strand in ['-', '+']:
        your_set += subset_ides(annotation, strand, training_ides)

    # 4. Extract the sequences from -10 to +20 for those genes in prots list. Beware the strand!
    minus10 = []
    plus20  = []
    for p in your_set:
        if annotation[p][2] == '+':
            seq = genome[annotation[p][0]-11:annotation[p][0]+22]
        else:
            seq = str(Seq(genome[annotation[p][1]-23:annotation[p][1]+10]).reverse_complement())

        m10 = seq[0:10]
        p20 = seq[13:]

        if len(m10) == 10 and check_nt(m10, False):
            if 'N' in m10:
                print m10
            minus10.append(m10)
        if len(p20) == 20 and check_nt(p20, False):
            plus20.append(p20)

    if random:
        for ide in training_ides:
            minus10.append(rand_dict[ide][:10])
            plus20.append(rand_dict[ide][10:30])

    # 8. create the pwm and the pssm
    # Done following the example in : http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec244

    background = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
    m = motifs.create(minus10)
    pwm10 = m.counts.normalize(pseudocounts=background)  # Pseudocounts based on the GC content (40%)
    pssm10 = pwm10.log_odds(background)
    if verbosity:
        print '\n-----\npwm for -10,\nbackground & pseudocounts:\nGC proportion\n\n', pssm10, '--> consensus:\n', pssm10.consensus

    m = motifs.create(plus20)
    pwm20 = m.counts.normalize(pseudocounts=background)
    pssm20 = pwm20.log_odds(background)
    if verbosity:
        print '\n-----\npwm for +20,\nbackground & pseudocounts:\nGC proportion\n\n', pssm20, '--> consensus:\n', pssm20.consensus

    # 9. The score for a sequences is computed as the product of probabilities per base per position
    # https://en.wikipedia.org/wiki/Position_weight_matrix
    # Transform the motifs to a dictionary {position:{A:weight, C:weight...}
    def motifscore(your_seq, your_pssm):
        """ Basic function to return the product of probabilities for a sequence given a motif """
        # to know the prob of having an A in the 2nd position --> pssm20['A',1]
        # using this:
        return sum([your_pssm[char, index] for char, index in zip(list(your_seq), range(0, len(your_seq)))])
        # When the PWM elements are calculated using log likelihoods, the score of a sequence can be
        # calculated by adding (rather than multiplying) the relevant values at each position in the PWM.
        # The sequence score gives an indication of how different the sequence is from a random sequence.
        # The score is 0 if the sequence has the same probability of being a functional site and of being a
        # random site. The score is greater than 0 if it is more likely to be a functional site than a random site, 
        # and less than 0 if it is more likely to be a random site than a functional site.[5] The sequence score
        # can also be interpreted in a physical framework as the binding energy for that sequence.


    # 9.2. If random sequences passed by a dictionary:
    results = {}
    if random==10:
        # FIX THIS
        for ide, seqi in testing_ides.iteritems():
            results[ide] = [motifscore(seqi[0:10], pssm10), motifscore(seqi[13:20], pssm20)]
    else:
        for ide in testing_ides:
            if annotation[ide][2]=='+':
                seq = genome[annotation[ide][0]-11:annotation[ide][0]+22]
            else:
                seq = str(Seq(genome[annotation[ide][1]-23:annotation[ide][1]+10]).reverse_complement())

            if check_nt(seq):
                results[ide] = [motifscore(seq[0:10], pssm10), motifscore(seq[13:], pssm20)]
            else:
                results[ide] = [0.0, 0.0]

        if verbosity:
            print '\n-----\nannotations evaluated: ', len(results)
    return results


def codons2freq(subseq, freqdic):
    codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                  'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                  'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
    pair = codontable[subseq[:3]]+codontable[subseq[3:]]
    try:
        return freqdic[pair]
    except:
        return 0.0


def NCterminal(feature_sequences, aa_considered=2):
    """ Given a list of amino acidic sequences computes the bias in the 2 C and N terminal aminoacids """

    alphabet  = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', '_']
    aas_freqs = {aa:0.0 for aa in [''.join(i) for i in itertools.product(alphabet, repeat=2)]}
    nt_freqs  = {aa:0.0 for aa in [''.join(i) for i in itertools.product(alphabet, repeat=2)]}
    ct_freqs  = {aa:0.0 for aa in [''.join(i) for i in itertools.product(alphabet, repeat=2)]}

    c = 0.0
    for seq in feature_sequences:

        # Bias in C and N counter
        ndiaa  = seq[1:3]
        cdiaa  = seq[-2:]
        nt_freqs[ndiaa] += 1
        ct_freqs[cdiaa] += 1

        # Split in diamino acids, update counts for frequency
        diaa_counts = Counter([seq[i:i+2] for i in xrange(0,len(seq)-1)])
        c += sum(diaa_counts.values())

        for diaa, counts in diaa_counts.iteritems():
            aas_freqs[diaa] += counts

    # Normalize by size all the dicts and apply background
    aas_freqs = {aa:counts/c for aa, counts in aas_freqs.iteritems()}
    seq_c     = sum(ct_freqs.values())
    nt_freqs  = {aa:math.log((counts/seq_c)/aas_freqs[aa], 2) if float(counts) > 0.0 else 0.0 for aa, counts in nt_freqs.iteritems()}
    ct_freqs  = {aa:math.log((counts/seq_c)/aas_freqs[aa], 2) if float(counts) > 0.0 else 0.0 for aa, counts in ct_freqs.iteritems()}

    return nt_freqs, ct_freqs


def return_frame(size, selected_frame):
    return range(selected_frame, size, 3)


def hexamer_freqs(your_sequences, background):
    # We define a big dictionary including the four possible measures:
    alphabets  = ['A', 'C', 'G', 'T']
    all_codons = set([''.join(x) for x in itertools.product(alphabets, repeat = 6)])
    lol        = [[], [], [], []]
    # Generate all the possible 6 windows and count for each frame all
    # at the same time you count
    for seq in your_sequences:
        seq  = seq[0:-3] # Do not count the stop
        size = len(seq)
        nhex = size-5 # Possible number of hexamers
        seq  = [seq[i:i+6] for i in xrange(nhex)]
        ### HEXAMERS ###
        lol[3] = lol[3]+seq
        ### n-HEXAMERS ###
        for fr in [0, 1, 2]:
            frame = return_frame(size, fr)
            lol[fr] = lol[fr]+[seq[item] for item in frame if item < len(seq)]
    # From list of lists to frequencies dictionary
    frequencies = [Counter(i) for i in lol]
    # Normalize
    def norm_dict(dic, all_codons, background):
        """
        To normalize we do two steps:
            - Extract a freq value dividing by the length
            - log(F(Obs)/F(Exp)) being F exp = 1/4096
        This function do two extra things:
            - add pseudocounts
            - add missing hexamers 
        """
        # Normalize 
        size  = sum(dic.values())
        new_d = {}
        for h in all_codons:
            if h in dic:
                new_d[h] = math.log( ((float(dic[h])/size) / background[h]) + 1.0 )
            else:
                new_d[h] = math.log( 1.0 )
        return new_d

    # Return the big last dict    
    frequencies = [norm_dict(dic, all_codons, background) for dic in frequencies]
    return frequencies


def hexamer_score(sequence, freq_dict):
    seq  = sequence[0:-3] # Do not count the stop
    size = len(seq)
    nhex = size-5 # Possible number of hexamers
    seq  = [seq[i:i+6] for i in xrange(nhex)]
    lol  = [[], [], [], []]

    ### HEXAMERS ###
    lol[3] = lol[3]+seq
    ### n-HEXAMERS ###
    for fr in [0, 1, 2]:
        frame = return_frame(size, fr)
        lol[fr] = lol[fr]+[seq[item] for item in frame if item < len(seq)]

    # Compute scores

    scores = [sum([freq_dict[i][hexamer] if check_nt(hexamer) else 0.0 for hexamer in lol[i]])/float(len(lol[i])) for i in [0, 1, 2, 3]]

    # Return the big last dict    
    return scores


# Start codon
def startcodon(seq, codons=10):
    first      = [seq[i:i+3] for i in range(0, len(seq), 3)][:codons]
    if 'ATG' in first:
        starttype = 1.0
    elif 'GTG' in first:
        starttype = 2.0
    elif 'TTG' in first:
        starttype = 3.0
    else:
        starttype = 4.0
    return starttype



##########################
# RANDOM FOREST FUNCTIONS#
##########################

def featurizer(organism, biased=None):
    """
    Given an organism code,
    generated the database with all their sequences and features to run a RF

    biased is used in order to bias the positive and feature set to detect small prots

    use_random_set = negative set will include random sequences (NOT TESTED)
    """

    # GENERATE YOUR GOLD SET

    # Some general variables:
    orgcode = organism[:5]

    # Load all the information
    try:
        # Getting back the objects:
        with open('./DBs/objs.pickle') as f:
            nt_seqs, aa_seqs, NCBI_all, anno_all, lengths, frames, contra_frames, rd_seqs, GCs, genome_lengths, conserved, correspondance = pickle.load(f)
            print 'DBs loaded'
    except:
        load_DBs(organisms)
        print 'DBs generated'
        with open('./DBs/objs.pickle') as f:
            nt_seqs, aa_seqs, NCBI_all, anno_all, lengths, frames, contra_frames, rd_seqs, GCs, genome_lengths, conserved, correspondance = pickle.load(f)
        print 'DBs loaded'

    # Randomly subset the feature set that will be used in generatin the hexamer and motif matrices:
    ncbi_set = [ide for ide in NCBI_all.keys() if ide.startswith(orgcode)]
    rand_set = [ide for ide in rd_seqs.keys() if ide.startswith('gc'+orgcode[:3])]

    # Subset 100 sequences for building the motifs and sequence stats that are relative, the rest will be our positive set
    random.shuffle(ncbi_set)
    if not biased:
        feature_set, positive_set = ncbi_set[:100], ncbi_set[100:]
        assert len(set(feature_set).intersection(set(positive_set))) == 0
    elif biased=='simple':
        big_prots = [ide for ide in ncbi_set if lengths[ide] >  100]
        sep_prots = [ide for ide in ncbi_set if lengths[ide] <= 100]
        feature_set  = sep_prots+ncbi_set[:100-len(sep_prots)]
        # positive_set = sep_prots+ncbi_set[100-len(sep_prots):]
        random.shuffle(ncbi_set)
        positive_set = sep_prots+ncbi_set[:100-len(sep_prots)]
    elif biased=='msd_biased':
        sep_prots = ['mpneu08047', 'mpneu08637', 'mpneu20382', 'mpneu06210', 'mpneu10761', 'mpneu01832', 'mpneu14551', 'mpneu07452', 'mpneu08850', 'mpneu29605', 'mpneu15883', 'mpneu23245', 'mpneu09464', 'mpneu13108', 'mpneu00732', 'mpneu24865', 'mpneu15316', 'mpneu15721', 'mpneu14366', 'mpneu20109', 'mpneu05846', 'mpneu21426', 'mpneu16655', 'mpneu09229', 'mpneu01006', 'mpneu11068', 'mpneu21675', 'mpneu09125', 'mpneu12014', 'mpneu11140', 'mpneu11154', 'mpneu21368', 'mpneu06215', 'mpneu07462', 'mpneu07169', 'mpneu15773', 'mpneu08806', 'mpneu21328', 'mpneu19587', 'mpneu19373', 'mpneu29058', 'mpneu06899', 'mpneu11482', 'mpneu12263', 'mpneu29683', 'mpneu18319', 'mpneu24711', 'mpneu16704', 'mpneu01113', 'mpneu25154', 'mpneu13213', 'mpneu08225', 'mpneu26462', 'mpneu11780', 'mpneu01851', 'mpneu11617', 'mpneu05172', 'mpneu01813', 'mpneu28053', 'mpneu01348', 'mpneu19950', 'mpneu03885', 'mpneu20162', 'mpneu07342', 'mpneu27456', 'mpneu27793', 'mpneu08784', 'mpneu14957', 'mpneu29796', 'mpneu03761', 'mpneu26219', 'mpneu24568', 'mpneu02833', 'mpneu12831', 'mpneu04516', 'mpneu28814', 'mpneu15228', 'mpneu06124', 'mpneu23596', 'mpneu16046', 'mpneu08005', 'mpneu13441', 'mpneu14706', 'mpneu11790', 'mpneu15713', 'mpneu14334', 'mpneu04086', 'mpneu20231', 'mpneu25438', 'mpneu08391', 'mpneu23078', 'mpneu05235', 'mpneu04363', 'mpneu23403', 'mpneu01628', 'mpneu12575', 'mpneu11095', 'mpneu28122', 'mpneu13480', 'mpneu26490', 'mpneu26399', 'mpneu23062', 'mpneu06729', 'mpneu25190', 'mpneu25551', 'mpneu01173', 'mpneu10189', 'mpneu09714', 'mpneu16481']
        feature_set  = sep_prots
        positive_set = sep_prots
    elif biased=='con_biased':
        sep_prots = ['mpneu12128', 'mpneu00486', 'mpneu00483', 'mpneu22706', 'mpneu22708', 'mpneu22709', 'mpneu00673', 'mpneu20109', 'mpneu00676', 'mpneu00677', 'mpneu00674', 'mpneu24099', 'mpneu19245', 'mpneu19244', 'mpneu10761', 'mpneu14135', 'mpneu18240', 'mpneu18241', 'mpneu10026', 'mpneu27767', 'mpneu29920', 'mpneu02353', 'mpneu12755', 'mpneu09553', 'mpneu21368', 'mpneu07507', 'mpneu06215', 'mpneu07502', 'mpneu09643', 'mpneu10031', 'mpneu02610', 'mpneu17858', 'mpneu09936', 'mpneu03808', 'mpneu22484', 'mpneu28703', 'mpneu13291', 'mpneu04436', 'mpneu07342', 'mpneu26386', 'mpneu00651', 'mpneu00653', 'mpneu26389', 'mpneu10552', 'mpneu22556', 'mpneu18069', 'mpneu10556', 'mpneu01221', 'mpneu24568', 'mpneu10559', 'mpneu13195', 'mpneu04121', 'mpneu13203', 'mpneu05037', 'mpneu28717', 'mpneu22248', 'mpneu25299', 'mpneu00649', 'mpneu28686', 'mpneu20231', 'mpneu28328', 'mpneu24094', 'mpneu27545', 'mpneu02787', 'mpneu00645', 'mpneu26391', 'mpneu11095', 'mpneu24346', 'mpneu26399', 'mpneu04794', 'mpneu06686', 'mpneu18376', 'mpneu00648', 'mpneu05666', 'mpneu08266', 'mpneu29642', 'mpneu21426', 'mpneu20025', 'mpneu05665', 'mpneu13062', 'mpneu17490', 'mpneu17656', 'mpneu28330', 'mpneu28332', 'mpneu24357', 'mpneu20028', 'mpneu18280', 'mpneu08404', 'mpneu27536', 'mpneu08403', 'mpneu17865', 'mpneu14339', 'mpneu10174', 'mpneu03621', 'mpneu11068', 'mpneu29651', 'mpneu09695', 'mpneu29306', 'mpneu08745', 'mpneu14607', 'mpneu08806', 'mpneu07542', 'mpneu18611', 'mpneu19236', 'mpneu30089', 'mpneu03352', 'mpneu04314', 'mpneu29325', 'mpneu29324', 'mpneu30086', 'mpneu29323', 'mpneu03632', 'mpneu03631', 'mpneu03637', 'mpneu14284', 'mpneu09037', 'mpneu09035', 'mpneu29468', 'mpneu29469', 'mpneu07877', 'mpneu08750', 'mpneu29915', 'mpneu11617', 'mpneu14611', 'mpneu07269', 'mpneu25937', 'mpneu17876', 'mpneu12706', 'mpneu05693', 'mpneu09191', 'mpneu29310', 'mpneu29312', 'mpneu06550', 'mpneu30093', 'mpneu23331', 'mpneu03644', 'mpneu03645', 'mpneu24337', 'mpneu24336', 'mpneu09209', 'mpneu13322', 'mpneu15228', 'mpneu23270', 'mpneu14334', 'mpneu05661', 'mpneu04441', 'mpneu04442', 'mpneu19233', 'mpneu05669', 'mpneu04449', 'mpneu04448', 'mpneu29315', 'mpneu18549', 'mpneu07078', 'mpneu12538', 'mpneu19259', 'mpneu19255', 'mpneu01173', 'mpneu09211', 'mpneu09210', 'mpneu10189', 'mpneu22564', 'mpneu22565', 'mpneu22562']
        feature_set  = sep_prots
        positive_set = sep_prots
    elif biased=='con_balanced_biased':
        sep_prots = ['mpneu12128', 'mpneu00486', 'mpneu00483', 'mpneu22706', 'mpneu22708', 'mpneu22709', 'mpneu00673', 'mpneu20109', 'mpneu00676', 'mpneu00677', 'mpneu00674', 'mpneu24099', 'mpneu19245', 'mpneu19244', 'mpneu10761', 'mpneu14135', 'mpneu18240', 'mpneu18241', 'mpneu10026', 'mpneu27767', 'mpneu29920', 'mpneu02353', 'mpneu12755', 'mpneu09553', 'mpneu21368', 'mpneu07507', 'mpneu06215', 'mpneu07502', 'mpneu09643', 'mpneu10031', 'mpneu02610', 'mpneu17858', 'mpneu09936', 'mpneu03808', 'mpneu22484', 'mpneu28703', 'mpneu13291', 'mpneu04436', 'mpneu07342', 'mpneu26386', 'mpneu00651', 'mpneu00653', 'mpneu26389', 'mpneu10552', 'mpneu22556', 'mpneu18069', 'mpneu10556', 'mpneu01221', 'mpneu24568', 'mpneu10559', 'mpneu13195', 'mpneu04121', 'mpneu13203', 'mpneu05037', 'mpneu28717', 'mpneu22248', 'mpneu25299', 'mpneu00649', 'mpneu28686', 'mpneu20231', 'mpneu28328', 'mpneu24094', 'mpneu27545', 'mpneu02787', 'mpneu00645', 'mpneu26391', 'mpneu11095', 'mpneu24346', 'mpneu26399', 'mpneu04794', 'mpneu06686', 'mpneu18376', 'mpneu00648', 'mpneu05666', 'mpneu08266', 'mpneu29642', 'mpneu21426', 'mpneu20025', 'mpneu05665', 'mpneu13062', 'mpneu17490', 'mpneu17656', 'mpneu28330', 'mpneu28332', 'mpneu24357', 'mpneu20028', 'mpneu18280', 'mpneu08404', 'mpneu27536', 'mpneu08403', 'mpneu17865', 'mpneu14339', 'mpneu10174', 'mpneu03621', 'mpneu11068', 'mpneu29651', 'mpneu09695', 'mpneu29306', 'mpneu08745', 'mpneu14607', 'mpneu08806', 'mpneu07542', 'mpneu18611', 'mpneu19236', 'mpneu30089', 'mpneu03352', 'mpneu04314', 'mpneu29325', 'mpneu29324', 'mpneu30086', 'mpneu29323', 'mpneu03632', 'mpneu03631', 'mpneu03637', 'mpneu14284', 'mpneu09037', 'mpneu09035', 'mpneu29468', 'mpneu29469', 'mpneu07877', 'mpneu08750', 'mpneu29915', 'mpneu11617', 'mpneu14611', 'mpneu07269', 'mpneu25937', 'mpneu17876', 'mpneu12706', 'mpneu05693', 'mpneu09191', 'mpneu29310', 'mpneu29312', 'mpneu06550', 'mpneu30093', 'mpneu23331', 'mpneu03644', 'mpneu03645', 'mpneu24337', 'mpneu24336', 'mpneu09209', 'mpneu13322', 'mpneu15228', 'mpneu23270', 'mpneu14334', 'mpneu05661', 'mpneu04441', 'mpneu04442', 'mpneu19233', 'mpneu05669', 'mpneu04449', 'mpneu04448', 'mpneu29315', 'mpneu18549', 'mpneu07078', 'mpneu12538', 'mpneu19259', 'mpneu19255', 'mpneu01173', 'mpneu09211', 'mpneu09210', 'mpneu10189', 'mpneu22564', 'mpneu22565', 'mpneu22562']
        random.shuffle(sep_prots)
        feature_set  = sep_prots[:int(len(sep_prots)/2)]+ncbi_set[:int(len(ncbi_set)/2)]
        positive_set = sep_prots[int(len(sep_prots)/2):]+ncbi_set[int(len(ncbi_set)/2):]
    elif biased=='m_biased':
        sep_prots = ['mpneu00927', 'mpneu15094', 'mpneu15493', 'mpneu06980', 'mpneu15042', 'mpneu08444', 'mpneu05637', 'mpneu11026', 'mpneu10150', 'mpneu24567', 'mpneu10953', 'mpneu12197', 'mpneu01157', 'mpneu03479', 'mpneu02755', 'mpneu02737', 'mpneu02736', 'mpneu01908', 'mpneu05383', 'mpneu17013', 'mpneu01547', 'mpneu20327', 'mpneu01759', 'mpneu25444', 'mpneu01756', 'mpneu20688', 'mpneu05613', 'mpneu01022', 'mpneu07365', 'mpneu05728', 'mpneu08209', 'mpneu05239', 'mpneu25104', 'mpneu01437', 'mpneu15466', 'mpneu28410', 'mpneu00918', 'mpneu25694', 'mpneu25691', 'mpneu15068', 'mpneu05207', 'mpneu00321', 'mpneu06303', 'mpneu15119', 'mpneu00925', 'mpneu26727', 'mpneu26720', 'mpneu26931', 'mpneu25420', 'mpneu25103', 'mpneu05606', 'mpneu08223', 'mpneu05612', 'mpneu27185', 'mpneu21086', 'mpneu11156', 'mpneu17014', 'mpneu15127', 'mpneu13388', 'mpneu25431', 'mpneu01119', 'mpneu06832', 'mpneu21555', 'mpneu07703', 'mpneu26344', 'mpneu26426', 'mpneu20674', 'mpneu06095', 'mpneu26745', 'mpneu11184', 'mpneu25406', 'mpneu15463', 'mpneu20673', 'mpneu11663', 'mpneu06089', 'mpneu06085', 'mpneu07399', 'mpneu20323', 'mpneu11353', 'mpneu25417', 'mpneu01160', 'mpneu10052', 'mpneu22177', 'mpneu21017', 'mpneu10054']
        random.shuffle(sep_prots)
        feature_set  = sep_prots[:int(len(sep_prots)/2)]+ncbi_set[:int(len(ncbi_set)/2)]
        positive_set = sep_prots[int(len(sep_prots)/2):]+ncbi_set[int(len(ncbi_set)/2):]
    elif biased=='random':
        feature_set  = rand_set[:100]
        positive_set = rand_set[100:]
    elif biased in ['A', 'C', 'E', 'D', 'G', 'F', '--', 'H', 'K', 'J', 'M', 'L', 'O', 'N', 'P', 'S', 'R', 'U', 'T', 'V', 'I']:
        mpn2luca = u.str_dic_generator('../smprots_DB/id_genes/mpneumoniae_pairs (copy).txt', 1, 0)
        COGs = {}
        setset = []
        for ide, cog in u.str_dic_generator('../../mycorepo/COGs.txt', 0, 1).iteritems():
            if ide in mpn2luca:
                if cog==biased:
                    setset.append(mpn2luca[ide])
        random.shuffle(setset)
        feature_set  = setset+ncbi_set[:50]
        positive_set = setset+ncbi_set[50:100]

    # Generate a fasta file with the sequences in the feature set
    with open('./temp/feat_'+organism+'_sequences.fa', 'w') as fo:
        for ide in feature_set:
            if biased != 'random':
                fo.write('>'+ide+'\n'+str(nt_seqs[ide])+'\n')
            else:
                if ide in nt_seqs:
                    fo.write('>'+ide+'\n'+str(nt_seqs[ide])+'\n')
                else:
                    fo.write('>'+ide+'\n'+str(rd_seqs[ide])+'\n')

    # Define the negative set and randomly select a numbe to balance the positive set
    noconserved_SEPs_set = [ide for ide, orto in load_number_times_conserved('../blaster/results_processed.out').iteritems() if int(orto)==0 and lengths[ide] <= 100 and frames[ide][3]=='NO']
    # noconserved_SEPs_set = [ide for ide in nt_seqs.keys() if ide.startswith(orgcode) and ide not in conserved and lengths[ide]<=100 and frames[ide][3] == 'NO']    # These are all non conserved and SEPs in NONOVERLAPPING regions
    if biased==100:
        # FIX THIS
        noconserved_SEPs_set += rd_seqs

    random.shuffle(noconserved_SEPs_set)
    negative_set = noconserved_SEPs_set[:len(positive_set)]

    # Generate dictionary
    gold_sets = {}
    for k in nt_seqs.keys():
        if k in positive_set:
            gold_sets[k] = '+'
        elif k in negative_set:
            gold_sets[k] = '-'
        elif k in feature_set:
            gold_sets[k] = 'f'
        else:
            gold_sets[k] = 0

    ##########################################
    # Define the values from the feature seqs
    # CAI:
    cai = CodonAdaptationIndex()
    cai.generate_index("./temp/feat_"+organism+"_sequences.fa")

    # 2 aa N and C terminal
    ccdfr = {k:float(v) for k, v in u.str_dic_generator('/home/smiravet/crg/mycorepo/ncterm/bulk2aa_mollicutes.csv',
                                                        key_index=0, value_index=1, header=False, split_by='\t').iteritems()}
    nnter = {k:(math.log(float(v)/ccdfr[k],2) if float(v) > 0 else 0) for k, v in u.str_dic_generator('/home/smiravet/crg/mycorepo/ncterm/nterm2aa_mollicutes.csv',
                                                                                                      key_index=0, value_index=1, header=False, split_by='\t').iteritems()}
    ccter = {k:(math.log(float(v)/ccdfr[k],2) if float(v) > 0 else 0) for k, v in u.str_dic_generator('/home/smiravet/crg/mycorepo/ncterm/cterm2aa_mollicutes.csv',
                                                                                                      key_index=0, value_index=1, header=False, split_by='\t').iteritems()}
    nns = {k:codons2freq(str(seq[3:9])  , nnter) for k, seq in nt_seqs.iteritems() if k.startswith(orgcode)}
    ccs = {k:codons2freq(str(seq[-9:-3]), ccter) for k, seq in nt_seqs.iteritems() if k.startswith(orgcode)}

    # -10 +20 motifs
    # random sequences
    put_seqs_ides  = [x for x in nt_seqs.keys() if x.startswith(orgcode)]
    if biased!='random':
        orfs_m10p20_sc = run_m10p20(organism, feature_set, put_seqs_ides)
    else:
        orfs_m10p20_sc = run_m10p20(organism, feature_set, put_seqs_ides, random=True, rand_dict=rd_seqs)
    minus10 = {k:float(v[0]) for k, v in orfs_m10p20_sc.iteritems()}
    plus20  = {k:float(v[1]) for k, v in orfs_m10p20_sc.iteritems()}

    # Stack energies
    stackRBS = {k:float(v) for k, v in u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/stack_E/'+organism+'_stackE_RBS.txt', 0, 1).iteritems()}
    stackm10 = {k:float(v) for k, v in u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/stack_E/'+organism+'_stackE_m10p20.txt', 0, 1).iteritems()}

    # Hexamers
    # define the background model
    alphabets         = ['A', 'C', 'G', 'T']
    GC_content        = float(np.mean([GCs[nt_seqs[ide]] for ide in set(positive_set).union(set(feature_set))]))/100.0
    proportions       = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
    background_model  = set([''.join(x) for x in itertools.product(alphabets, repeat = 6)])
    background_model  = {cd:reduce(lambda x, y: x*y, [proportions[x] for x in cd]) for cd in background_model}
    assert 0.98 < sum(background_model.values()) <1.01

    hex_freq_dict = hexamer_freqs([str(x) for x in u.load_multifasta("./temp/feat_"+organism+"_sequences.fa").values()], background_model)

    hex0, hex1, hex2, hex3 = {}, {}, {},{}
    for ide, seq in {k:nt_seqs[k] for k in put_seqs_ides}.iteritems():
        scores = hexamer_score(seq, hex_freq_dict)
        hex0[ide] = scores[0]
        hex1[ide] = scores[1]
        hex2[ide] = scores[2]
        hex3[ide] = scores[3]

    ############################
    # FEATURIZE ALL THE SEQUENCES

    # Label the samples
    # Combine them in a dictionary with the labels 1=positive , 0=negative
    columns = ['ide',
               'start_codon', 'GC'               , 'CAI'                ,
               'Nterminal'  , 'Cterminal'        ,
               'm10_score'  , 'p20_score'        , 'm10p20_stack_energy', 'RBS_stack_energy',
               'n_hexamer'  , 'dicodon_frequency', '1_hexamer'          , '2_hexamer'       ,
               'set_type'   ]
    to_append = []

    for ide, seq in {k:nt_seqs[k] for k in put_seqs_ides}.iteritems():
        st            = startcodon(seq)
        gc            = GC(seq)
        your_cai      = cai.cai_for_gene(seq)
        nnf           = nns[ide]
        ccf           = ccs[ide]
        m10           = minus10[ide]
        p20           = plus20[ide]
        m10p20_stackE = stackm10[ide]
        RBS_stackE    = stackRBS[ide]
        your_hexn     = hex3[ide]
        your_hex0     = hex0[ide] # dicodon frequency
        your_hex1     = hex1[ide]
        your_hex2     = hex2[ide]
        to_append.append([ide, st, gc, your_cai, nnf, ccf, m10, p20, m10p20_stackE, RBS_stackE, your_hexn, your_hex0, your_hex1, your_hex2, gold_sets[ide]])

    to_append.sort(key=lambda x: x[0])
    df = pd.DataFrame(to_append, columns=columns)
    df = df.set_index(df['ide']) # vertical index = organism

    return df

def load_RF_features(organism, all_new=False, extension='', biased=None):

    # Load all the information or create it if it does not exist
    if not os.path.exists('./DBs/featurized_'+extension+organism+'.pickle') or all_new==True:
        print 'computing features'
        features_df = featurizer(organism, biased=biased)
        print 'features computed'
        with open('./DBs/featurized_'+extension+organism+'.pickle', 'w') as f:
            pickle.dump(features_df, f)
        print 'features stored'
    else:
        # Getting back the objects:
        with open('./DBs/featurized_'+extension+organism+'.pickle') as f:
            features_df = pickle.load(f)
        print 'Features loaded'

    return features_df


def feature_weights(rf, feat_dict, weightsxfold, add='', add_string=''):
    """
    Set of functions to explore the weights of the different
    features in a random forest classifier
    """

    # Measure importances
    importances = rf.feature_importances_

    fo = open('./results/'+add+'weights.txt','a')
    fo.write("Feature ranking tree "+add_string+" :\n----------\n")
    inv_map = {v: k for k, v in feat_dict.iteritems()}
    to_write = []
    for i in range(0, len(importances)):
        to_write.append([inv_map[i], importances[i]])
    to_write.sort(key=lambda x:x[1], reverse=True) # Sort by second element
    c = 1
    for el in to_write:
        fo.write(str(c)+'.'+el[0]+' '+str(el[1])+'\n')
        c += 1

        if el[0] not in weightsxfold:
            weightsxfold[el[0]]=[el[1]]
        else:
            weightsxfold[el[0]].append(el[1])

    fo.write('\n')
    fo.close()

    return weightsxfold


def print_weights(weightsxfold, add=''):
    """
    Given the weights per fold, prints the varplot for it (including bar errors)
    """

    h, ax3 = plt.subplots()
    dic = {k:[np.mean(v), np.std(v)] for k, v in weightsxfold.iteritems()}

    feats = []
    means = []
    stdev = []
    for feat, stats in sorted(dic.items(), key=lambda i: i[1][0]):
        feats.append(feat)
        means.append(stats[0])
        stdev.append(stats[1])

    # Plot
    ax3.barh(range(len(feats)), means, color='g', xerr=stdev, ecolor='g', align='center', alpha=0.5)
    plt.yticks(range(len(feats)), feats, fontsize = 'large')
    ax3.set_ylim([-1, len(feats)])
    ax3.set_xlabel('Variance explained', fontsize = 'large')
    ax3.set_title('weights, folds:'+str(len(weightsxfold.values()[0])), fontsize = 'large')
    h.savefig('./results/'+add+'_weights_barplot_'+str(len(weightsxfold.values()[0]))+'.svg')


def random_forest(X, y, feat_dict, to_exclude=None, folds=25, test_size=0.2, add='', organism='',
                  n_estimators=200, oob_score=1, n_jobs=-1, random_state=50, max_features="auto", min_samples_leaf=5):

    # Remove features we don't want in the analysis
    if to_exclude != None:
        X = u.remove_column(X, to_exclude)

    # Define the classifier:
    classifier = RandomForestClassifier(n_estimators=n_estimators, oob_score=oob_score, n_jobs=n_jobs, random_state=random_state, max_features=max_features, min_samples_leaf=min_samples_leaf)

    # Run the classifier in different manners depending on the folds selected
    if folds == 0 or folds == None:
        classifier.fit(X, y)
    else:
        # Run the classifier with kfolds and ROC curves
        cv = StratifiedShuffleSplit(y, n_splits=folds, test_size=test_size)

        # To compute the roc and auc:
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)

        # To plot the weigths
        if os.path.exists('./results/'+add+'weights.txt'):
            os.system('rm ./results/'+add+'weights.txt')      # Require as writing is at the end of the file not in a new
        weightsxfold = {}

        # ROC and PR curves
        f, ax1 = plt.subplots()
        g, ax2 = plt.subplots()
        i = 1
        for train, test in cv.split(X, y):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            probs = classifier.fit(X_train, y_train).predict_proba(X_test)
            # ROC and area under the curve
            fpr, tpr, thresholds = roc_curve(y_test, probs[:,1])
            mean_tpr += interp(mean_fpr, fpr, tpr)
            roc_auc   = auc(fpr, tpr)
            # Precssion recall
            precision, recall, thresholds = precision_recall_curve(y_test, probs[:,1])

            #Study the weight of each feature:
            add_string = '_'+str(i)
            feature_weights(classifier, feat_dict, weightsxfold, add+'_'+str(folds) ,add_string)

            # Plot them
            ax1.plot(fpr, tpr, lw=3, alpha=0.5)
            ax2.plot(recall, precision, lw=3, alpha=0.5)

            i += 1

        # plot the random hypothesis for the ROC 
        ax1.plot([0,1], [0,1], '--', lw=3, color=(0.6, 0.6, 0.6), label='Random')

        # Basic stats plus plot of the mean roc
        mean_tpr /= float(folds)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        ax1.plot(mean_fpr, mean_tpr, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=3)

        ######################
        # Plot the roc of the classifier
        ax1.set_xlim([-0.05, 1.05])
        ax1.set_ylim([-0.05, 1.05])
        ax1.set_xlabel('False Positive Rate')
        ax1.set_ylabel('True Positive Rate')
        ax1.set_title('ROC, folds = '+str(folds)+', '+organism)
        ax1.legend(loc="lower right")
        nm = './results/'+add+'_roc_'+str(folds)+'.svg'
        f.savefig(nm)

        # PR curve
        ax2.set_xlim([-0.05, 1.05])
        ax2.set_ylim([-0.05, 1.05])
        ax2.set_xlabel('Recall')
        ax2.set_ylabel('Precision')
        ax2.set_title('Precision-Recall Curve, folds = '+str(folds)+', '+organism)
        ax2.legend(loc="lower right")
        nm = './results/'+add+'_pr_'+str(folds)+'.svg'
        g.savefig(nm)

        #####################
        # Plot the weights:
        print_weights(weightsxfold, add=add)

    # return class
    return classifier


def RFClassifier(organism, feats_to_exclude=[], folds=25, test_size=0.2, extension='', all_new=False, additional_directory='', biased=None):
    """
    Run a whole RF analysis over the organism given
    """
    # Some general variables:
    print 'Training a classifier for '+organism
    orgcode = organism[:5]

    #Load features:
    features_df = load_RF_features(organism, all_new=all_new, extension=extension, biased=biased)

    # The same with the classifier and if not created, MAKE IT!
    if not os.path.exists('./DBs/clf_'+extension+'_'+organism+'_'+str(folds)+'.pickle') or all_new==True:
        print 'Creating classifier'
        # Extract X and y (features and labels)
        k1 = features_df.loc[(features_df.set_type == '+') | (features_df.set_type == '-')]

        features  = np.array(k1.iloc[:,1:-1])  # Extract an array without ide and set_type info
        labels  = np.array([1 if i=='+' else 0 for i in k1.iloc[:,-1]])    # Labels
        print 'Balanced set of '+str(len(labels))+' samples (50% positive - 50% negative cases)'

        # Names of the features
        # Check if required to exclude
        feat_names   = [str(m) for m in k1.iloc[:,1:-1].columns]
        to_exclude   = u.indexes(feat_names, feats_to_exclude)    # With this you not require order
        # Update
        feat_names = [feat for feat in feat_names if feat not in feats_to_exclude]
        numbers    = range(0,len(feat_names))
        feat_dict  = u.lists2dict(feat_names, numbers)
        print 'Features used: \n\n--> '+'\n--> '.join(feat_names)+'\n'

        # Create the classifier (will generate the figures in parallel)
        print 'Training your classifier...'

        ####
        clf = random_forest(X=features, y=labels, feat_dict=feat_dict, to_exclude=to_exclude, folds=folds, test_size=test_size, add=extension, organism=organism)  # Organism just for some naming
        ####

        with open('./DBs/clf_'+extension+'_'+organism+'_'+str(folds)+'.pickle', 'w') as f:
            pickle.dump(clf, f)
        print 'Classifier stored in ./DBs/clf_'+extension+'_'+organism+'_'+str(folds)+'.pickle'
    else:
        # Getting back the objects:
        with open('./DBs/clf_'+extension+'_'+organism+'_'+str(folds)+'.pickle') as f:
            clf = pickle.load(f)
        print './DBs/clf_'+extension+'_'+organism+'_'+str(folds)+'.pickle LOADED!'

    # Classify annotation
    print "Let's classify "+str(len(features_df))+' proteins!\nSubsetting your prots...'

    # Separate between ides and features arrays
    your_ides = np.array(features_df.iloc[:,0])      # Only ides
    your_luca = u.str_dic_generator(filename='/home/smiravet/crg/antprotect/RF_final/DBs/correspondance.txt', key_index=0, value_index=1)
    luca_ides = np.array([your_luca[i] for i in your_ides])

    to_classify = np.array(features_df.iloc[:,1:-1])   # Remove assigned class and the column of ides

    # Remove not used columns
    feat_names   = [str(m) for m in features_df.iloc[:,1:-1].columns]
    to_exclude   = u.indexes(feat_names, feats_to_exclude)    # With this you not require order
    to_classify  = u.remove_column(to_classify, to_exclude)

    # Predict probs and classes
    print "Predicting..."
    class_pred = clf.predict(to_classify)
    probs_pred = clf.predict_proba(to_classify)
    print "All predicted :D!\nLet's write the results in a file..."

    clf_results    = np.column_stack((your_ides, luca_ides, probs_pred, class_pred))
    cols           = ['ide', 'alt_ide', 'prob_0', 'prob_1', 'pred_class']
    if biased=='random':
        clf_results    = np.column_stack((your_ides, luca_ides, class_pred))
        cols           = ['ide', 'alt_ide', 'pred_class']
    clf_results_df = pd.DataFrame(clf_results, columns=cols)

    # Store everything in a csv
    clf_results_df.to_csv('./results/'+additional_directory+extension+'_'+organism+'_'+str(folds)+'_classified.csv', sep='\t', index=False)

    # Return the results
    return clf_results_df


#########
# NEW RANDOM FOREST RANSEPS
def create_directory_structure(project_folder):
    """ Given a project directory creates a subtree directory """

    os.mkdir(project_folder)

    for subfolder in ['classification_stats', 'classification_tasks', 'classifiers', 'feature_importances', 'features', 'results']:
        os.mkdir(project_folder+'/'+subfolder)


def compute_m10p20_dict(organism, nt_seqs, feature_set, annotation, random=False ):
    orfs_m10p20_sc = run_m10p20(organism, feature_set, nt_seqs.keys(), annotation)
    minus10 = {k:float(v[0]) for k, v in orfs_m10p20_sc.iteritems()}
    plus20  = {k:float(v[1]) for k, v in orfs_m10p20_sc.iteritems()}

    return minus10, plus20


def return_propy_nms():

    Des = GetProDes('MKKKKKMMAKKMAMMA')  # To have the names
    alldes=Des.GetALL()
    propy_nms, propy_des = zip(*sorted([[i, j] for i, j in alldes.iteritems()]))
    return ['ide']+list(propy_nms)


def all_propy(ide, seq):
    """
    Return ALL names and features for an aa feature if well formated,
    0.0 otherwise
    """
    try:
        Des = GetProDes(seq)
        alldes=Des.GetALL()
        propy_nms, propy_des = zip(*sorted([[i, j] for i, j in alldes.iteritems()]))
    except:
        propy_des = [0.0]*1537

    return [ide]+list(propy_des)


def parallel_propy_featurizer(aa_seqs, organism, output_dir='/media/smiravet/SAMHD/DBs/propy_features/', jobs=6):
    try:
        propy_df = pd.read_pickle(output_dir+organism+'_propy.pickle')
    except:
        # Generate the list of list with names
        propy_nms = return_propy_nms()
        # Generate the list of lists with features
        print 'Running featurization on '+str(len(aa_seqs))+' sequences for '+organism+' in '+str(jobs)+' cores...'
        to_append = Parallel(n_jobs=jobs)(delayed(all_propy)(ide, seq) for ide, seq in aa_seqs.iteritems())
        # Generate the df
        to_append.sort(key=lambda x: x[0])
        df = pd.DataFrame(to_append, columns=propy_nms)
        propy_df = df.set_index(df['ide']) # vertical index = organism

        # Store it
        propy_df.to_pickle(output_dir+organism+'_propy.pickle')
    'Finished!'

    return propy_df


def propy_organism(input_file, size=15):
    """
    Allows the parallel computation of propy feature databases
    Using external scripts ;)
    """
    organism = input_file.split('/')[-1].split('_')[0]
    aa_seqs  = {k:v for k, v in u.load_multifasta(input_file).iteritems() if len(v) >= size}

    propy_featurizer(aa_seqs, organism)
    print organism+' finished'



def short3propy(aa_seqs, organism, output_dir='/media/smiravet/SAMHD/DBs/propy_features3/'):
    """ Short version to compute only the 3 most important features in propy """

    if os.path.isfile(output_dir+organism+'_propy.pickle'):
        propy_df = pd.read_pickle(output_dir+organism+'_propy.pickle')
    else:
        propy_nms = ['ide', 'tausw9', '_HydrophobicityD1001', '_SecondaryStrD1001']
        to_append = []
        for ide, seq in aa_seqs.iteritems():
            l = len(seq)
            # to_append.append([ ide, SOSW(seq)['tausw9'] , DHF(seq)['_HydrophobicityD1001'] , DSS(seq)['_SecondaryStrD1001'] ])
            # to_append.append([ ide, 0.0 , 0.0 , 0.0 ])
            to_append.append([ ide, SOSW(seq)['tausw9']/l, DHF(seq)['_HydrophobicityD1001']/l, DSS(seq)['_SecondaryStrD1001']/l])

        # Create file
        df = pd.DataFrame(to_append, columns=propy_nms)
        propy_df = df.set_index(df['ide']) # vertical index = organism
        propy_df.drop('ide', axis=1, inplace=True)

        # Store it
        propy_df.to_pickle(output_dir+organism+'_propy.pickle')

    return propy_df


def featurizer2(organism, nt_seqs, aa_seqs, annotation, positive_set, feature_set, negative_set, extension, propy_feats=None, project_folder=None):
    """
    New version of featurizer,
    This function requires the user to provide two dictionaries of nt and aa seqs with ides and a list of
    ids that will be used to generate the positive, negative and feature set.
    Organism is required to load the genome and compute some features
    Annotation refers to a dict with ide:[st, en, strand]
    Extension is used to store the feature DBs
    """
    # Define the storing working directory
    if not project_folder:
        project_folder = '/home/smiravet/crg/antprotect/RF_final/temp/'
    else:
        project_folder += 'features/'

    # Detect if additional organisms in the set
    orgcodes = set([k[:5] for k in aa_seqs.keys()])
    if len(orgcodes) == 2:
        additional_organisms = True
        nt_seqsA, aa_seqsA, annotationA, nt_seqsB, aa_seqsB, annotationB = {},{},{},{},{},{}
        feature_setA, feature_setB = [], []
        if 'ecoli' in orgcodes:
            organismB = 'ecoli'
        else:
            organismB = 'bsubtilis'
        for k, v in annotation.iteritems():
            if k[:5] == organism[:5]:
                nt_seqsA[k]    = nt_seqs[k]
                aa_seqsA[k]    = aa_seqs[k]
                annotationA[k] = v
                feature_setA.append(k)
            else:
                nt_seqsB[k]    = nt_seqs[k]
                aa_seqsB[k]    = aa_seqs[k]
                annotationB[k] = v
                feature_setB.append(k)
    else:
        additional_organisms = False


    # Generate a fasta file with the sequences in the feature set
    with open(project_folder+'feat_'+organism+'_'+extension+'_sequences.fa', 'w') as fo:
        for ide in feature_set:
            fo.write('>'+ide+'\n'+str(nt_seqs[ide])+'\n')

    # Generate dictionary
    gold_sets = {k:'' for k in nt_seqs.keys()}
    full_sets = set(positive_set).union(set(feature_set), set(negative_set))

    for k in gold_sets.keys():
        if k in full_sets:
            if k in positive_set:
                gold_sets[k] += '+'
            if k in negative_set:
                gold_sets[k] += '-'
            if k in feature_set:
                gold_sets[k] += 'f'
        else:
            gold_sets[k] = 0

    ##########################################
    # Define the values from the feature seqs
    # GC:
    GCs = {ide:GC(seq) for ide, seq in nt_seqs.iteritems()}

    # CAI:
    cai = CodonAdaptationIndex()
    cai.generate_index(project_folder+'feat_'+organism+'_'+extension+'_sequences.fa')

    # 2 aa N and C terminal
    nt_diaa_freq, ct_diaa_freq = NCterminal([aa_seqs[ide] for ide in feature_set])
    nns = {ide:nt_diaa_freq[seq[1:3]]  for ide, seq in aa_seqs.iteritems()}  # Do not count the metionin
    ccs = {ide:ct_diaa_freq[seq[-2:]] for ide, seq in aa_seqs.iteritems()}

    # -10 +20 motifs and stack energies
    if not additional_organisms:
        minus10, plus20   = compute_m10p20_dict(organism, nt_seqs, feature_set, annotation)
        stackEEs          = generate_DB_stack_energy(organism, annotation, extension)
    else:
        minus10, plus20   = compute_m10p20_dict(organism, nt_seqsA, feature_setA, annotationA)
        stackEEs          = generate_DB_stack_energy(organism, annotationA, extension)
        minus10B, plus20B = compute_m10p20_dict(organismB, nt_seqsB, feature_setB, annotationB)
        stackEEsB         = generate_DB_stack_energy(organismB, annotationB, extension)

        for k in annotationB.keys():
            minus10[k]  = minus10B[k]
            plus20[k]   = plus20B[k]
            stackEEs[k] = stackEEsB[k]

    # Hexamers
    # define the background model
    alphabets         = ['A', 'C', 'G', 'T']
    GC_content        = float(np.mean([GCs[ide] for ide in set(positive_set).union(set(feature_set))]))/100.0
    proportions       = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
    background_model  = set([''.join(x) for x in itertools.product(alphabets, repeat = 6)])
    background_model  = {cd:reduce(lambda x, y: x*y, [proportions[x] for x in cd]) for cd in background_model}
    assert 0.98 < sum(background_model.values()) <1.01

    hex_freq_dict = hexamer_freqs([str(x) for x in u.load_multifasta(project_folder+'feat_'+organism+'_'+extension+'_sequences.fa').values()], background_model)

    hex0, hex1, hex2, hex3 = {}, {}, {},{}
    for ide, seq in nt_seqs.iteritems():
        scores = hexamer_score(seq, hex_freq_dict)
        hex0[ide] = scores[0]
        hex1[ide] = scores[1]
        hex2[ide] = scores[2]
        hex3[ide] = scores[3]

    ############################
    # FEATURIZE ALL THE SEQUENCES

    # Label the samples
    # Combine them in a dictionary with the labels 1=positive , 0=negative
    to_append = []
    for ide, seq in nt_seqs.iteritems():
        st            = startcodon(seq)
        gc            = GCs[ide]
        your_cai      = cai.cai_for_gene(seq)
        nnf           = nns[ide]
        ccf           = ccs[ide]
        m10           = minus10[ide]
        p20           = plus20[ide]
        RBS_stackE    = stackEEs[ide][0]
        RBS01         = stackEEs[ide][2]
        m10p20_stackE = stackEEs[ide][1]
        your_hexn     = hex3[ide]
        your_hex0     = hex0[ide] # dicodon frequency
        your_hex1     = hex1[ide]
        your_hex2     = hex2[ide]

        to_append.append([ide, st, gc, your_cai, nnf, ccf, m10, p20, m10p20_stackE, RBS_stackE, RBS01, your_hexn, your_hex0, your_hex1, your_hex2, gold_sets[ide]])

    # Create the df
    columns = ['ide',
               'start_codon', 'GC'               , 'CAI'                ,
               'Nterminal'  , 'Cterminal'        ,
               'm10_score'  , 'p20_score'        , 'm10p20_stack_energy', 'RBS_stack_energy', 'RBS_presence',
               'n_hexamer'  , 'dicodon_frequency', '1_hexamer'          , '2_hexamer'       ,
               'set_type']

    to_append.sort(key=lambda x: x[0])
    df = pd.DataFrame(to_append, columns=columns)
    feats_df = df.set_index(df['ide']) # vertical index = organism

    # PROPY COMPUTATION (CONSIDER REWRITING)
    if propy_feats:
        # Load
        if propy_feats in ['all', 'subset1']:
            propy_db = pd.read_pickle('/media/smiravet/SAMHD/DBs/propy_features/'+organism+'_propy.pickle')
            # Subset features
            if propy_feats=='all':
                propy_db = propy_db.drop('ide', 1)
            elif propy_feats=='subset1':
                propy_db = propy_db[['tausw9', 'taugrant15', '_SolventAccessibilityD2001', '_HydrophobicityD1001', 'QSOgrant47',
                                     '_PolarityD1001', '_NormalizedVDWVD3001', '_PolarizabilityD3001', '_SecondaryStrD1001',
                                     '_ChargeD2001']]
        else:
            propy_db = short3propy(aa_seqs, organism, project_folder)
            # This will be alreatd [tausw9, hydrophibicityD1001 and SecondaryStrD1001]

        # Remove possible duplicate
        if not propy_db.index.is_unique:
            print 'Removing dupiclates...'
            propy_db.drop_duplicates(inplace=True)

        # Concate and remove possible Nan (propy has everything!)
        feats_df = pd.concat([feats_df.iloc[:,:-1], propy_db, feats_df.iloc[:,-1:]], axis=1)
        feats_df = feats_df.dropna()

    # Store the pickle object
    feats_df.to_pickle(project_folder+'feat_'+organism+'_'+extension+'.pickle')

    return feats_df


def miniRanSeps(organism           , nt_seqs                   , aa_seqs                        , annotation       ,
                positive_set       , feature_set               , negative_set                   , to_exclude=[]    ,
                sfolds=0           , test_size=0.2             , random_state_test=None         ,
                n_estimators=200   , oob_score=1               , n_jobs=-1                      , random_state=None, max_depth=None, max_features="auto", min_samples_leaf=5,
                propy_feats=None   , extension=''              , project_folder='/tmp/'):
    """
    Run a whole RF analysis over the organism given, same idea than the previous version but
    this is thought to be iterative ;)

    sfolds represents how many times you repeat the classification test OVER THE SAME DATASET, not required as the
    kfolds is done by the iteration of this function.
    """

    # Some basics + Featurization
    orgcode = organism[:5]
    features_df = featurizer2(organism, nt_seqs, aa_seqs, annotation, positive_set, feature_set, negative_set, extension, propy_feats, project_folder)

    # Remove features we don't want in the analysis
    if len(to_exclude) > 0:
        for te in to_exclude:
            features_df = features_df.drop(te, 1)   # axis number == 1 (0 for rows and 1 for columns.)


    # CREATE IT!
    # Extract X and y (features and labels)
    k1 = features_df.loc[(features_df.set_type == '+') | (features_df.set_type == '+f') | (features_df.set_type == '-')]
    features  = np.array(k1.iloc[:,1:-1])  # Extract an array without ide and set_type info
    labels  = np.array([1 if i in ['+', '+f'] else 0 for i in k1.iloc[:,-1]])    # Labels

    # Names of the features
    # Check if required to exclude
    feat_names = [str(m) for m in k1.iloc[:,1:-1].columns]
    numbers    = range(0,len(feat_names))
    feat_dict  = dict(zip(numbers, feat_names))
    # print 'Features used: \n\n--> '+'\n--> '.join(feat_names)+'\n'

    # Define the classifier
    clf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, oob_score=oob_score, n_jobs=n_jobs, random_state=random_state, max_features=max_features, min_samples_leaf=min_samples_leaf)

    # Statistics
    fprs, tprs = [], []
    recalls, precisions = [], []

    if test_size > 0:
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)
        sfolds   = 1 if sfolds == 0 else sfolds   # To save code in the iteration
        for i in range(0, sfolds):
            # Split the datasets
            X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=test_size, random_state=random_state_test)
            # Train and predict
            probs                            = clf.fit(X_train, y_train).predict_proba(X_test)
            # ROC curve components
            fpr, tpr, thresholds             = roc_curve(y_test, probs[:,1])
            mean_tpr                        += interp(mean_fpr, fpr, tpr)
            roc_auc                          = auc(fpr, tpr)                        # If you need the value of AUC for that classification
            # Precision recall
            precision, recall, thresholds    = precision_recall_curve(y_test, probs[:,1])

            fprs.append(fpr)
            tprs.append(tpr)
            recalls.append(recall)
            precisions.append(precision)

        fpr_tpr = [np.mean(fprs, axis=0), np.mean(tprs, axis=0)]
        recall_precision = [np.mean(recalls, axis=0), np.mean(precisions, axis=0)]

        mean_tpr /= float(sfolds)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
    else:
        # TO avoid problems if not test_size defined
        fpr_tpr, recall_precision, mean_tpr, mean_auc = None, None, None, None


    # Create the classifier (will generate the figures in parallel)
    clf.fit(features, labels)
    with open(project_folder+'classifiers/clf_'+organism+'_'+extension+'.pickle', 'w') as f:
        pickle.dump(clf, f)

    # FEATURE IMPORTANCES
    # Measure importances
    importances = clf.feature_importances_
    to_write = []
    for i in range(0, len(importances)):
        to_write.append([feat_dict[i], importances[i]])
    to_write.sort(key=lambda x:x[1], reverse=True) # Sort by second element
    c = 1
    with open(project_folder+'feature_importances/'+organism+'_'+extension+'_weights.txt','w') as fo:
        fo.write("Feature ranking tree "+organism+'_'+extension+" :\n----------\n")
        for el in to_write:
            fo.write(str(c)+'.'+el[0]+'\t'+str(el[1])+'\n')
            c += 1

    # CLASSIFY!
    # Separate between ides and features arrays
    your_ides = np.array(features_df.iloc[:,0])      # Only ides
    to_classify = np.array(features_df.iloc[:,1:-1])   # Remove assigned class and the column of ides

    # Predict probs and classes
    class_pred = clf.predict(to_classify)
    probs_pred = clf.predict_proba(to_classify)

    clf_results    = np.column_stack((your_ides, probs_pred, class_pred))
    cols           = ['ide', 'prob_0', 'prob_1', 'pred_class']
    clf_results_df = pd.DataFrame(clf_results, columns=cols)
    # Store everything in a csv
    clf_results_df.to_csv(project_folder+'classification_tasks/'+organism+'_'+extension+'_classified.csv', sep='\t', index=False)

    # Return the results
    return clf, clf_results_df, fpr_tpr, recall_precision


def load_RF_weights(folder):
    weights_dict = {}
    started = False
    for fil in glob.glob(folder+'*_weights.txt'):
        with open(fil, 'r') as fi:
            for line in fi:
                if line[0].isdigit():
                    line = line.strip().split()
                    if started:
                        weights_dict[line[0].split('.')[-1]].append(float(line[1]))
                    else:
                        weights_dict[line[0].split('.')[-1]] = [float(line[1])]
        started = True

    # Extract information
    weights_dict = {k:[np.mean(v), np.std(v)] for k, v in weights_dict.iteritems()}
    return weights_dict


def parametrize_estimators(organism                        , nt_seqs                           , aa_seqs          , annotation            ,
                           autoset=[None, None, None, None], set_sizes=[None, None, None, None],
                           positive_set=None               , feature_set=None                  , negative_set=None, to_exclude=[]         ,
                           folds=25                        , sfolds=0                          , test_size=0.2    , random_state_test=None,
                           n_estimators=200                , oob_score=1                       , n_jobs=-1        , random_state=None     , max_depth=None, max_features="auto", min_samples_leaf=5,
                           propy_feats =None               , extension=''                      , project_folder='/tmp/'):

    """ Iteratively run RanSeps and explore the number of estimators required """

    # Create the hierarchical tree of directories to write the results
    if project_folder[-1]!='/':
        create_directory_structure(project_folder)   # os.mkdir needs the path with no '/' at the end
        project_folder += '/'
    else:
        create_directory_structure(project_folder[:-1])

    # Prepare the basics for autodefine sets if required
    if autoset[0]!=None:
        big_prots = [ide for ide in autoset[1] if len(aa_seqs[ide]) >  100]
        sep_prots = [ide for ide in autoset[1] if len(aa_seqs[ide]) <= 100]
        if len(sep_prots) < 25:
            sys.exit('No enough SEPs for train the classifier')
        neg_set_prots = autoset[2]
        if autoset[3]:
            # Can be none
            neg_add_prots = autoset[3]

    # Define the basics if test_size defined (ROC and PR curves)
    if test_size:
        froc, rocax = plt.subplots()
        fprc, prcax = plt.subplots()
        stats2write = []
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)


    # ITERATE AS MANY TIMES AS FOLDS ARE
    error_rate=[]
    folds    = 1 if folds == 0 else folds   # To save code in the iteration, ensure at least one iteration is performed.
    for i in range(0, folds):
        # Define the sets if required
        if autoset[0]!=None:
            # Define the number of seps to extract for the POSITIVE SET
            random.shuffle(big_prots)
            random.shuffle(sep_prots)
            nseps        = int(set_sizes[0]*autoset[0])
            positive_sep = sep_prots[:nseps]
            positive_big = big_prots[:set_sizes[0]-nseps]
            positive_set = positive_sep+positive_big

            # Define the number of seps to extract for the FEATURE SET
            random.shuffle(big_prots)
            random.shuffle(sep_prots)
            nseps        = int(set_sizes[1]*autoset[0])
            feature_sep  = sep_prots[:nseps]
            feature_big  = big_prots[:set_sizes[1]-nseps]
            feature_set  = feature_sep+feature_big

            # Define the negative set
            random.shuffle(neg_set_prots)
            negative_set = neg_set_prots[:set_sizes[2]]
            if autoset[-1]:
                random.shuffle(neg_add_prots)
                negative_set += neg_add_prots[:set_sizes[3]]

        # Run the classification
        clf, clf_results_df, fpr_tpr, recall_precision = miniRanSeps(organism            , nt_seqs       , aa_seqs          , annotation,
                                                                     positive_set        , feature_set   , negative_set     , to_exclude,
                                                                     sfolds              , test_size     , random_state_test,
                                                                     n_estimators        , oob_score     , n_jobs           , random_state, max_depth, max_features, min_samples_leaf,
                                                                     propy_feats         , extension+'_'+str(i+1), project_folder)
        oob_error = 1 - clf.oob_score_
        error_rate.append(oob_error)

    return error_rate


def process_clf_results(project_folder, organism):
    """ Process all the results of multiple classifications """

    probess = {}
    for fil in glob.glob(project_folder+'classification_tasks/*.csv'):
        prob1     = u.str_dic_generator(fil, 0, -2, header=True)
        for k, v in prob1.iteritems():
            if k in probess:
                probess[k].append(float(v))
            else:
                probess[k] = [float(v)]
    probess = {k:[np.mean(v), np.std(v)] for k, v in probess.iteritems()}
    with open(project_folder+'results/'+organism+'_final_prediction.txt', 'w') as fo:
        od = collections.OrderedDict(sorted(probess.items()))
        for k, v in od.iteritems():
            fo.write(k+'\t'+'\t'.join([str(val) for val in v])+'\n')


def RanSEPs(organism                        , nt_seqs                           , aa_seqs               , annotation            ,
            autoset=[None, None, None, None], set_sizes=[None, None, None, None],
            positive_set=None               , feature_set=None                  , negative_set=None     , to_exclude=[]         ,
            folds=25                        , sfolds=0                          , test_size=0.2         , random_state_test=None,
            n_estimators=200                , oob_score=1                       , n_jobs=-1             , random_state=None     , max_depth=None, max_features="auto", min_samples_leaf=5,
            propy_feats =None               , extension=''                      , project_folder='/tmp/'):

    """
    Iteratively run miniRanSeps & extract results
    autoset allows to generate the positive, negative and feature set automatically with sizes included in set sizes.
        autoset[0] = Value of autoset defines the percentage of small proteins in the positive and feature set.
        autoset[1] = list of identifiers known as annotated
        autoset[2] = list of identifiers suitable to be NEGATIVE (ex: no conserved and no overlapping
        autoset[3] = list of additional negative set in case you want to do fancy stuff combining set sizes (small big for example)
    """

    # Create the hierarchical tree of directories to write the results
    if project_folder[-1]!='/':
        create_directory_structure(project_folder)   # os.mkdir needs the path with no '/' at the end
        project_folder += '/'
    else:
        create_directory_structure(project_folder[:-1])

    # Store all the parameters in shake of reproducibility
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    with open(project_folder+'parameters.txt', 'w') as fo:
        fo.write('function name "%s"\n' % inspect.getframeinfo(frame)[2])
        for i in args[4:]:
            if values[i]:
                if isinstance(values[i], list):
                    for j in values[i]:
                        if j:
                            fo.write("    %s = %s\n" % (i, j))
                        else:
                            fo.write("    %s = None\n" % (i))
                else:
                    if values[i]:
                        fo.write("    %s = %s\n" % (i, values[i]))
                    else:
                        fo.write("    %s = None\n" % (i))

    # Prepare the basics for autodefine sets if required
    if autoset[0]!=None:
        big_prots = [ide for ide in autoset[1] if len(aa_seqs[ide]) >  100]
        sep_prots = [ide for ide in autoset[1] if len(aa_seqs[ide]) <= 100]
        neg_set_prots = autoset[2]
        if autoset[3]:
            # Can be none
            neg_add_prots = autoset[3]

    # Define the basics if test_size defined (ROC and PR curves)
    if test_size:
        froc, rocax = plt.subplots()
        fprc, prcax = plt.subplots()
        stats2write = []
        mean_tpr = 0.0
        mean_fpr = np.linspace(0, 1, 100)


    # ITERATE AS MANY TIMES AS FOLDS ARE
    print 'Training a classifier for '+organism+'\nNumber of sequences provided:\t'+str(len(nt_seqs))
    print 'Prepare for training the classfier with:\n -->'+str(set_sizes[0])+' positive sequences\n -->'+str(set_sizes[2])+' negative sequences\n'
    folds    = 1 if folds == 0 else folds   # To save code in the iteration, ensure at least one iteration is performed.
    for i in range(0, folds):
        print 'Iteration '+str(i)+'...'
        # Define the sets if required
        if autoset[0]!=None:
            # Define the number of seps to extract for the POSITIVE SET
            random.shuffle(big_prots)
            random.shuffle(sep_prots)
            nseps        = int(set_sizes[0]*autoset[0])
            positive_sep = sep_prots[:nseps]
            positive_big = big_prots[:set_sizes[0]-nseps]
            positive_set = positive_sep+positive_big

            # Define the number of seps to extract for the FEATURE SET
            random.shuffle(big_prots)
            random.shuffle(sep_prots)
            nseps        = int(set_sizes[1]*autoset[0])
            feature_sep  = sep_prots[:nseps]
            feature_big  = big_prots[:set_sizes[1]-nseps]
            feature_set  = feature_sep+feature_big

            # Define the negative set
            random.shuffle(neg_set_prots)
            negative_set = neg_set_prots[:set_sizes[2]]
            if autoset[-1]:
                random.shuffle(neg_add_prots)
                negative_set += neg_add_prots[:set_sizes[3]]

        # Run the classification
        clf, clf_results_df, fpr_tpr, recall_precision = miniRanSeps(organism            , nt_seqs       , aa_seqs          , annotation,
                                                                     positive_set        , feature_set   , negative_set     , to_exclude,
                                                                     sfolds              , test_size     , random_state_test,
                                                                     n_estimators        , oob_score     , n_jobs           , random_state, max_depth, max_features, min_samples_leaf,
                                                                     propy_feats         , extension+'_'+str(i+1), project_folder)
        # Plot if required
        if test_size:
            mean_tpr += interp(mean_fpr, fpr_tpr[0], fpr_tpr[1])    # Updated each iteration
            rocax.plot(fpr_tpr[0]         , fpr_tpr[1]         , lw=3, alpha=0.5)
            prcax.plot(recall_precision[0], recall_precision[1], lw=3, alpha=0.5)

            stats2write.append([str(fpr_tpr[0]), str(fpr_tpr[1]), str(recall_precision[0]), str(recall_precision[0])])

    ### GENERAL
    print 'Computing basic statistics\n'
    if test_size:
        # Write the results:
        with open(project_folder+'classification_stats/'+organism+'_'+extension+'_stats_'+str(folds)+'.txt', 'w') as fo:
            fo.write('FPR\tTPR\tRECALL\tPRECISION\n')
            for your_stats in stats2write:
                fo.write('\t'.join(your_stats)+'\n')
        # PR curve
        prcax.set_xlim([-0.05, 1.05])
        prcax.set_ylim([-0.05, 1.05])
        prcax.set_xlabel('Recall')
        prcax.set_ylabel('Precision')
        prcax.set_title('Precision-Recall Curve, folds = '+str(folds)+', '+organism)
        prcax.legend(loc="lower right")
        fprc.savefig(project_folder+'results/'+organism+'_'+extension+'_PR_'+str(folds)+'.svg')

        # ROC curve
        rocax.plot([0,1], [0,1], '--', lw=3, color=(0.6, 0.6, 0.6), label='Random')   # Random hypothesis

        mean_tpr    /= float(folds)
        mean_tpr[-1] = 1.0
        mean_auc     = auc(mean_fpr, mean_tpr)
        rocax.plot(mean_fpr, mean_tpr, 'k--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=5)

        rocax.set_xlim([-0.05, 1.05])
        rocax.set_ylim([-0.05, 1.05])
        rocax.set_xlabel('False Positive Rate')
        rocax.set_ylabel('True Positive Rate')
        rocax.set_title('ROC, folds = '+str(folds)+', '+organism)
        rocax.legend(loc="lower right")
        froc.savefig(project_folder+'results/'+organism+'_'+extension+'_ROC_'+str(folds)+'.svg')


    # Plot stuff and statistics
    ##########################
    # WEIGHT IMPORTANCE FIGURE

    feats, means, stdev = [], [], []
    weights_dict        = load_RF_weights(project_folder+'feature_importances/')
    wfolds              = str(len(weights_dict.values()[0]))     # To provide annotation
    for feat, stats in sorted(weights_dict.items(), key=lambda i: i[1][0]):
        feats.append(feat)
        means.append(stats[0])
        stdev.append(stats[1])

    #Plot
    fwei, weiax = plt.subplots()
    weiax.barh(range(len(feats)), means, color='g', align='center', alpha=0.4,
               xerr=stdev, error_kw={'ecolor':'g', 'elinewidth':5})
    plt.yticks(range(len(feats)), feats, fontsize = 'medium')
    weiax.set_ylim([-1, len(feats)])
    weiax.set_xlabel('Variance explained', fontsize = 'medium')
    weiax.set_title('Weights, folds:'+wfolds, fontsize = 'medium')
    fwei.savefig(project_folder+'results/'+organism+'_'+extension+'_weights_'+str(folds)+'.svg')

    print "All predicted :D!\nLet's write the results in a file located in "+project_folder+"results/"+organism+"_final_prediction.txt..."

    # Process results (mean of prob1!)
    process_clf_results(project_folder, organism)

    print "Thanks for using RanSEPs :D!\n||Stand and be True||"


# REUSE CLASSIFIERS WITH ANY PROVIDED SEQUENCE

def only_classify(organism, nt_seqs, aa_seqs, project_folder, annotation=None, extension='it1_s_', save_all=None):
    """
    Classify a set of sequences using an already trained set of classifiers
    save_all=if something, will keep dataframes in features and update the final prediction, for example: 10,
    this will be used in the new name where you store the stuff
    Annotation is not require by default but WARNING:
    classification uses multiple annotation context dependent features (-10 for example)
    that are NOT CONSIDERED at this level if you do not specify an annotation dictionary
    """

    prediction = {k:[] for k in nt_seqs.keys()}

    # Parse the classifiers
    if project_folder[-1]!='/':
        project_folder+='/'
    classifiers = glob.glob(project_folder+'classifiers/*.pickle')

    # Propy features:
    TAUs = {ide:SOSW(seq)['tausw9'] for ide, seq in aa_seqs.iteritems()}
    HYDs = {ide:DHF(seq)['_HydrophobicityD1001'] for ide, seq in aa_seqs.iteritems()}
    SECs = {ide:DSS(seq)['_SecondaryStrD1001'] for ide, seq in aa_seqs.iteritems()}

    # GC:
    GCs = {ide:GC(seq) for ide, seq in nt_seqs.iteritems()}

    # Hexamers background model
    alphabets = ['A', 'C', 'G', 'T']

    # Original annotation
    if annotation:
        organism_annotation = u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+organism+'_annotation.txt') # Load original annotation
        organism_annotation.update(annotation)

    # ITERATE CLASSIFYING AND STORING NEW STUFF
    for clf_handle in classifiers:
        with open(clf_handle) as f:
            clf = pickle.load(f)
        iteration = clf_handle.split('/')[-1].split(extension)[-1].split('.')[0]

        # PREPARE THE FEATURIZATION:
        feature_seqs_file = project_folder+'features/feat_'+organism+'_'+extension+iteration+'_sequences.fa'
        feature_seqs = u.load_multifasta(feature_seqs_file)

        # Specific features
        # CAI
        cai = CodonAdaptationIndex()
        cai.generate_index(feature_seqs_file)

        # 2 aa N and C terminal
        nt_diaa_freq, ct_diaa_freq = NCterminal(feature_seqs.values())
        nns = {ide:nt_diaa_freq[seq[1:3]] for ide, seq in aa_seqs.iteritems()}
        ccs = {ide:ct_diaa_freq[seq[-2:]] for ide, seq in aa_seqs.iteritems()}

        # -10 and +20 feature related and Stacking energies
        if annotation:
            motif_seqs      = dict(nt_seqs.items() + feature_seqs.items())
            minus10, plus20 = compute_m10p20_dict(organism, motif_seqs, feature_seqs.keys(), organism_annotation)
            stackEEs        = generate_DB_stack_energy(organism, annotation, extension)
        else:
            minus10, plus20, stackEEs = {}, {}, {}
            minus10  = {ide:0.0 for ide in nt_seqs.keys()}
            plus20   = {ide:0.0 for ide in nt_seqs.keys()}
            stackEEs = {ide:[0.0, 0.0, 0.0] for ide in nt_seqs.keys()}

        # Hexamers
        GC_content        = float(np.mean([GC(seq) for ide, seq in feature_seqs.iteritems()]))/100.0
        proportions       = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
        background_model  = set([''.join(x) for x in itertools.product(alphabets, repeat = 6)])
        background_model  = {cd:reduce(lambda x, y: x*y, [proportions[x] for x in cd]) for cd in background_model}

        hex_freq_dict = hexamer_freqs([str(x) for x in feature_seqs.values()], background_model)
        hex0, hex1, hex2, hex3 = {}, {}, {},{}
        for ide, seq in nt_seqs.iteritems():
            scores = hexamer_score(seq, hex_freq_dict)
            hex0[ide] = scores[0]
            hex1[ide] = scores[1]
            hex2[ide] = scores[2]
            hex3[ide] = scores[3]

        # Label the samples
        # Combine them in a dictionary with the labels 1=positive , 0=negative
        ides_to_clf, to_classify, to_save = [], [], []
        for ide, seq in nt_seqs.iteritems():
            st            = startcodon(seq)
            gc            = GCs[ide]
            your_cai      = cai.cai_for_gene(seq)
            nnf           = nns[ide]
            ccf           = ccs[ide]
            m10           = minus10[ide]
            p20           = plus20[ide]
            RBS_stackE    = stackEEs[ide][0]
            RBS01         = stackEEs[ide][2]
            m10p20_stackE = stackEEs[ide][1]
            your_hexn     = hex3[ide]
            your_hex0     = hex0[ide] # dicodon frequency
            your_hex1     = hex1[ide]
            your_hex2     = hex2[ide]
            tau           = TAUs[ide]
            hyd           = HYDs[ide]
            sec           = SECs[ide]
            to_append = [st, gc, your_cai, nnf, ccf, m10, p20, m10p20_stackE, RBS_stackE, RBS01, your_hexn, your_hex0, your_hex1, your_hex2, tau, hyd, sec]
            if save_all:
                to_save.append([ide]+to_append+[0]) # Last 0 is to classify
            to_classify.append(to_append)
            ides_to_clf.append(ide)

        # Create the df
        if save_all:
            columns = ['ide'        ,
                       'start_codon', 'GC'                  , 'CAI'                ,
                       'Nterminal'  , 'Cterminal'           ,
                       'm10_score'  , 'p20_score'           , 'm10p20_stack_energy', 'RBS_stack_energy', 'RBS_presence',
                       'n_hexamer'  , 'dicodon_frequency'   , '1_hexamer'          , '2_hexamer'       ,
                       'tausw9'     , '_HydrophobicityD1001', '_SecondaryStrD1001' ,
                       'set_type']
            df = pd.DataFrame(to_save, columns=columns)
            df = df.set_index(df['ide'])
            df.to_pickle(project_folder+'features/feat_'+organism+'_'+save_all+'_'+extension+iteration+'.pickle')

        ### PREDICTION
        # Predict:
        class_pred  = clf.predict(to_classify)
        probs_pred  = clf.predict_proba(to_classify)
        clf_results = np.column_stack((ides_to_clf, probs_pred, class_pred))

        # Dictionary to return 
        for ide, prob0, prob1, pred_class in clf_results:
            prediction[ide].append(float(prob1))

        # Store it if required
        if save_all:
            cols = ['ide', 'prob_0', 'prob_1', 'pred_class']
            clf_results_df = pd.DataFrame(clf_results, columns=cols)
            clf_results_df.to_csv(project_folder+'classification_tasks/'+organism+'_'+save_all+'_'+extension+iteration+'_classified.csv', sep='\t', index=False)

    # Generate the final prediction and return results if necessary:
    if save_all:
        process_clf_results(project_folder, organism)

    # Process prediction
    return {k:[np.mean(v), np.std(v)] for k, v in prediction.iteritems()}


# TEST by FRAGMENTS
def subseqs(nt_seq, aa_seq, length):
    i = random.randint(0, len(aa_seq) - length + 1)
    return [nt_seq[i*3:(i+length)*3], aa_seq[i:i+length]]


def fragment_sequences(nt_seqs, aa_seqs, size, n):
    """ Return n pairs of sequences nt and aa with aa <size> """

    # Ensure proteins are big enough to sample
    sequences = np.random.choice(nt_seqs.keys(), n)
    resulting_seqs = []

    for ide in sequences:
        resulting_seqs.append(subseqs(nt_seqs[ide], aa_seqs[ide], size))

    return resulting_seqs


##########################
#     OTHER FUNCTIONS    #
##########################

def find_alt_start_codon(seq, starts=['ATG', 'GTG', 'TTG']):
    codons = [seq[i:i+3] for i in range(3, len(seq), 3)]   # Do not take into consideration the first codon (original start)
    alt_starts = [i*3 for i in range(0, len(codons)) if codons[i] in starts]
    return alt_starts

def shortener(nt_seqs, aa_seqs, annotation, size=19):
    """
    Given 3 dictionaries with nt and aa seqs and their annotation,
    Add the shortened version of the sequences to the dictionary
    """

    new_nt_seqs, new_aa_seqs, new_annotation = {}, {}, {}

    for ide, seq in nt_seqs.iteritems():
        # Keep track of the annotation
        new_nt_seqs[ide]    = seq
        new_aa_seqs[ide]    = aa_seqs[ide]
        new_annotation[ide] = annotation[ide]
        st, en, strand      = sorted(annotation[ide][:2])+[annotation[ide][-1]]

        # Compute new starts
        alt_starts = find_alt_start_codon(seq)

        if len(alt_starts) != 0:
            c = 1
            for alt_st in alt_starts:
                extension = 'a'+str(c)
                new_ide   = ide+extension
                new_nts   = seq[alt_st:]
                new_aas   = 'M'+aa_seqs[ide][(alt_st/3)+1:]
                if strand == '+':
                    new_ann = [st+alt_st, en, strand]
                else:
                    new_ann = [st, en-alt_st, strand]

                if len(new_aas) >= size:
                    new_nt_seqs[new_ide]    = new_nts
                    new_aa_seqs[new_ide]    = new_aas
                    new_annotation[new_ide] = new_ann
                    c += 1

    return new_nt_seqs, new_aa_seqs, new_annotation

##########################
#  STATISTICS FUNCTIONS  #
##########################

def list_NA_generator(filename, index=0):
    """ index refers to the column with the identifier """
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


def samples_iterator(m):
    """ Define the identifiers found in each sample """
    samples = range(4, m)
    for i in samples:
        yield list_NA_generator('/home/smiravet/crg/antprotect/RF_final/DBs/MS_outcome_study/MS_116.csv', i)


def MS_iterator(n, m=121):
    """
    n =  sampling size
    m = last index where a sample is found
    """

    samples = [s for s in samples_iterator(m)]

    for i in range(m):
        sampling_set = set([])
        for j in np.random.choice(samples, n, replace=False):
            sampling_set = sampling_set.union(set(j))
        yield sampling_set

##########################
# CONSERVATION FUNCTIONS #
##########################


def process_blastDB(size):
    all_lengths = {}
    all_seqs = u.load_multifasta('/media/smiravet/SAMHD/metablast_results/all.fa')
    for ide, seq in all_seqs.iteritems():
        if len(seq) >= size:
            all_lengths[ide] = len(seq)

    return all_lengths


def load_blastDB_lengths(size):
    try:
        all_lengths = {k:int(v) for k, v in u.str_dic_generator('./DBs/all_lengths.txt', 0, 1).iteritems()}
    except:
        all_lengths = process_blastDB(size)
        with open('./DBs/all_lengths.txt', 'w') as fo:
            for ide, lens in all_lengths.iteritems():
                fo.write(ide+'\t'+str(lens)+'\n')

    return all_lengths


def load_all_orthologs(blast_fil, evalth=2e-8, alignth=50.0, lenth=58.0, identh=50.0, close=[], size=19):
    """
    Counts the number of times a protein is conserved based on several threshold
    filters:
        - evalue
        - alignment length (percentage of the protein that is properly aligned)
        - % of the length difference between the two annotations
        - % identity

    Close allows to discard close species in the counting
    """

    organism = blast_fil.split('/')[-1].split('_')[0][:5]
    if organism == 'resul':
        organism = 'mpneu'
    close += ['mferi', organism]

    lengths = load_blastDB_lengths(size)
    results = {k:[] for k in lengths.keys() if k.startswith(organism)}
    with open(blast_fil) as fi:
        for line in fi:
            line = line.strip().split()
            target  = line[0]
            query = line[1][:5]
            query_c = line[1]

            if query not in close and organism[:5] not in query and target in lengths and query_c in lengths:
                evl  = float(line[-2])  # eval
                ale  = int(line[3])-int(line[5])  # alignment length - gap opens
                idp  = float(line[2])

                ale = 100.0*float(ale)/lengths[target]
                ple = 100.0*(max(lengths[target],lengths[query_c])-abs(lengths[target]-lengths[query_c]))/max(lengths[target],lengths[query_c])
                if evl <= evalth and ale >= alignth and ple >= lenth and idp >= identh:
                    results[target] += [query_c]

    return {k:set(v) for k, v in results.iteritems()}


def load_number_times_conserved(blast_fil, evalth=2e-8, alignth=50.0, lenth=58.0, identh=50.0, close=[], size=19):
    """
    Counts the number of times a protein is conserved based on several threshold
    filters:
        - evalue
        - alignment length (percentage of the protein that is properly aligned)
        - % of the length difference between the two annotations
        - % identity

    Close allows to discard close species in the counting
    """
    pre_results = load_all_orthologs(blast_fil, evalth, alignth, lenth, identh, close, size)
    results = {}
    for ide, orthologs in pre_results.iteritems():
        results[ide] = len(set([v[:5] for v in orthologs]))

    return results


def return_orthologs(species_codeA, species_codeB):
    """
    Given a species codeA,
    return a dictionary {ide protein: [orthologs]}

    species_codeB refers to the organism you compare
    """

    xl = pd.ExcelFile('./DBs/orthology.xlsx')
    df = xl.parse("orthology")

    # Subset columns and remove rows where species A has no protein:
    df = df[[species_codeA, species_codeB]]
    df = df[df[species_codeA] != '*']

    orth_dict = {}

    for ideA, ideB in zip(df[species_codeA], df[species_codeB]):
        ideA, ideB = str(ideA), str(ideB)
        ideA = ideA.split(',')
        for i in ideA:
            orth_dict[i] = ideB.split(',')

    return orth_dict

def find_closest_homologs():
    """
    Compare all pairs of species to return a dictionary 
    """


def PAL2NAL(aa_alg_seq, nt_seq):
    """
    Given a pair of matched aa and nt seqs, aa seq in aligned format (NNN---NNN-N)
    Returns the nt seq including the gaps and respecting the codons
    """

    nt_seq = u.splitn_str(nt_seq, 3)
    nt_alg_seq = ''
    c = 0
    for i in range(0, len(aa_alg_seq)):
        if aa_alg_seq[i] == '-':
            nt_alg_seq += 'XXX'
        else:
            nt_alg_seq += nt_seq[c]
            c += 1

    return nt_alg_seq


def ntseq2conservation(nt_seqA, aa_seqA, nt_seqB, aa_seqB):
    """
    Given two sequences in aa and nt format, being A: interest sequence, B: reference protein
    Return a binary vector representing the conservation of each base in the nt sequence.

    The algorithm works aligning the aa sequences (supposed to be more conserved than nt sequence),
    that aligned sequences are used as reference to align the nt sequences.
    Last step compare the PAL2NAL transformed sequences to generate a vector where each bit represents
    a base of the nt sequence:
    0 - NO CONSERVED
    1 - CONSERVED
    """

    #1. Align amino acids
    matrix   = matlist.pam30
    best_alg = pairwise2.align.globaldx(aa_seqA, aa_seqB, matrix)[0]

    #2. Use that alignment as reference to align the nt sequences
    nt_seqA_alg = PAL2NAL(best_alg[0], nt_seqA)
    nt_seqB_alg = PAL2NAL(best_alg[1], nt_seqB)

    #3. Transform the interest protein to a binary string
    binary_vector = []
    nt_seqA_alg   = u.splitn_str(nt_seqA_alg, 3)
    nt_seqB_alg   = u.splitn_str(nt_seqB_alg, 3)

    for Acodon, Bcodon in zip(nt_seqA_alg, nt_seqB_alg):
        if Acodon != 'XXX':
            for i in [0,1,2]:
                if Acodon[i] == Bcodon[i]:
                    binary_vector.append(1)
                else:
                    binary_vector.append(0)

    return binary_vector


def decompose_vector(array):
    """
    Given an array of ones and zeroes
    Return the sum of 1s in 1st, 2nd and 3rd position normalized by the length (len/3)
    """

    # extract the length after ensuring it is divisible by 3
    l      = len(array)
    assert l%3 == 0
    norm   = l/3.0

    result = []
    for m in [0,1,2]:
        result.append(sum([array[n] for n in range(m, l, 3)])/norm)

    return result



def compute_conservation_ratios(ref_start, ref_stop, binary_vector, sep_start, sep_stop, strand):
    """
    Given a pair of set of coordinates, a binary vector and a strand
    Return the noverlap_vector and the overlap_vector arrays of conservation example

    Iput ref start 1, ref stop 12, sep start 3, sep stop 8, strand +:

    A 1  2  3  4  5  6  7  8  9  10 11 12
    B       1  2  3  4  5  6
      1  1  1  1  1  1  1  1  0  1  1  0  # Binary vector

    1. Compute separatedly
    overlap_vector: 1 1 1 1 1 1
    nonoverlap    : 1 1 0 1 1 0     # As we are always in triplets, we basically merge as a 2nd position will be always the second if you remove n*3 characters in the middle

    2. Decompose vectors:
    over = [1.0, 1.0, 1.0]
    nono = [1.0, 1.0, 0.0]

    PAY ATTENTION THAT THE OVERLAPPING VECTOR HAS TO BE RE-ARRANGED
    """

    if strand == '+':
        # Translate to vector scale:
        # ref start --> 1, etc
        sep_start = sep_start - ref_start + 1
        sep_stop  = sep_stop  - ref_start + 1
        ref_stop  = ref_stop  - ref_start + 1
        ref_start = 1

        # Divide the vector between overlapping and not overlapping
        OVER_vector = decompose_vector(binary_vector[sep_start-1:sep_stop])
        NOOV_vector = decompose_vector(binary_vector[:sep_start-1]+binary_vector[sep_stop:])

        # The OVER vector is not sorted properly, depending on the frame the positions have to be moved
        frame = return_frame_plus(ref_start, sep_start)
        if frame == 1:
            OVER_vector = [OVER_vector[-1],] + OVER_vector[:2]
        elif frame == 2:
            OVER_vector = OVER_vector[1:]    + [OVER_vector[0],]
        else:
            # Should not occur
            OVER_vector = ['problem','formatting','frame']
    else:
        pass

    return [NOOV_vector, OVER_vector]



def blaster(your_file, your_code, minus=19, threshold=0.00001, directory='/home/smiravet/crg/antprotect/RF_final/metaRF/blaster/', itself=False):
    """
    (modification tobe usable in a pipeline)
    minus allows to subset the sequences with length of aa minus <= length <= plus
    """

    # Starting with a directory merge all of them into a database
    if not os.path.exists(directory+'all.fa'):
        print('concatenating SEPs files...')
        cmd = 'cat /home/smiravet/crg/dbs/rbs_orfs/*_small_aa.fa > '+directory+'all.fa'
        os.system(cmd)

    # Filter by size target file
    if not os.path.exists(directory+your_code+'_target_filtered.fa'):
        print('filtering target fa by minus='+str(minus)+'...')
        fo = open(directory+your_code+'_target_filtered.fa', 'w')
        with open(your_file, 'r') as fi:
            for line in fi:
                if line.startswith('>'):
                    ide = line
                else:
                    if minus <= len(line.strip()):
                        fo.write(ide+line)
        fo.close()

    # Reformat the db file to ensure each sequence is just a line satisfying the size
    if not os.path.exists(directory+'db_filtered.fa'):
        print('reformating db ...')
        fo = open(directory+'db_filtered.fa', 'w')
        for ide, seq in u.load_multifasta(directory+'all.fa').iteritems():
            if minus <= len(seq):
                fo.write('>'+ide+'\n'+seq+'\n')
        fo.close()

    # Generate the database required to run the blas analysis
    if not os.path.exists(directory+'db_filtered.fa.phr'):
        print('making the database...')
        cmd = 'makeblastdb -in '+directory+'db_filtered.fa -dbtype prot -parse_seqids'
        os.system(cmd)

    # Run blastp
    print('running blast...')
    cmd = 'blastp -db '+directory+'db_filtered.fa -query '+directory+your_code+'_target_filtered.fa -out '+directory+'results.out -outfmt 7 -evalue '+str(threshold)+' -num_threads 5'
    os.system(cmd)

    # Process the results file
    if not os.path.exists(directory+your_code+'_output.out'):
        print('processing results...')
        cmd = "sed '/^#/ d' < "+directory+"results.out > "+directory+your_code+"_results_processed.out"
        os.system(cmd)
        print 'threshold given = '+str(threshold)

        results = {key:[] for key in u.load_multifasta(directory+your_code+'_target_filtered.fa')}
        with open(directory+your_code+"_results_processed.out", 'r') as fi:
            print threshold
            for line in fi:
                blastresults = line.strip().split()
                target = blastresults[0]
                query  = blastresults[1][:5]
                if itself:
                    results[target] += [query]   # It does not matter the intraspecific match
                else:
                    if query != target[:5]:
                        results[target] += [query]

        # Write the file
        fo = open(directory+your_code+'_output.out', 'w')
        for ide, cons in results.iteritems():
            fo.write(ide+'\t'+str(len(set(cons)))+'\n')
        fo.close()

def conservationALL(inFiles='/home/smiravet/crg/dbs/rbs_orfs/*_small_aa.fa', directory='/home/smiravet/crg/antprotect/RF_final/metaRF/blaster/', size=19):
    for fil in glob.glob(inFiles):
        print fil
        your_code = fil.split('/')[-1].replace('_small_aa.fa', '')
        blaster(fil, your_code, minus=size, threshold=0.00001, directory=directory)


########### METAFEATURES STUDIES


def load_mean_array_organism(organism, filters={'minsize':20, 'maxsize':100, 'RF':0.0}):
    """
    Generate a dataframe with the mean values for each
    protein features
    """

    list_of_files = glob.glob('/media/smiravet/SAMHD/metaRFs/it1_s_'+organism+'/features/feat_'+organism+'_it1_s_*.pickle')

    # Start the df:
    df = pd.read_pickle(list_of_files[0]).drop('ide', 1)
    # Iterate
    for fil in list_of_files[1:]:
        df2 = pd.read_pickle(fil).drop('ide', 1)
        df = pd.concat((df, df2))
    # Do the mean
    df = df.groupby(df.index).mean()

    # Add the lengths colums and return it
    probess = {k:float(v) for k, v in u.str_dic_generator('/media/smiravet/SAMHD/metaRFs/it1_s_'+organism+'/results/'+organism+'_final_prediction.txt', 0, 1).iteritems()}
    target = {ide:int(length) for ide, length in u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/aa_lengths/'+organism+'_lengths.txt', 0, 1).iteritems() if filters['minsize'] <= int(length) <= filters['maxsize']}
    if filters['RF'] > 0.0:
        target = {ide:target[ide] for ide, prob in probess.iteritems() if ide in target and prob >= filters['RF']}
    lengths = pd.DataFrame.from_dict(target, orient='index')
    lengths.columns = ['length']

    # Concat and return
    filt_df = pd.concat([df, lengths], axis=1).dropna()
    return filt_df.mean(0).to_frame().transpose().rename(index={0:organism})


def features_values_per_organism(organisms, name, min_size_th=20, max_size_th= 100, RF_th=0.0):
    """
    Given a directory with subfolders of RF results
    return a single dataframe with the mean value for
    the feature values
    """
    filename = './DBs/'+name+'.pickle'
    if os.path.isfile(filename):
        df = pd.read_pickle(filename)
    else:
        started = False
        for organism in organisms:
            if not started:
                df  = load_mean_array_organism(organism, {'minsize':min_size_th, 'maxsize':max_size_th, 'RF':RF_th})
                started = True
            else:
                df2 = load_mean_array_organism(organism, {'minsize':min_size_th, 'maxsize':max_size_th, 'RF':RF_th})
                df  = df.append(df2)
        df.to_pickle(filename)
    return df


############ REPEATED REGIONS

def repeated(ann_list, wdw_size=50):

    st, en = sorted(ann_list[:2])
    annotation = set(range(st, en+1))
    n_wdws = float(en-st-wdw_size)           # Number of possible windows in that gene

    # Count the number of windows repeated
    repeated = []
    for i in u.list_generator('/home/smiravet/crg/mycorepo/repeated_'+str(wdw_size)+'.txt', 0):
        st, en = [int(pos) for pos in i.split('-')]
        repeated += range(st, en)
    repeated = set(repeated)

    n = len(annotation.intersection(repeated))-wdw_size+1.0

    ratio = n/n_wdws*100.0

    if ratio <= 0:
        ratio = 0.0

    return ratio


def repeated_mm(ann_list, wdw_size=50, mm=1):

    st, en = sorted(ann_list[:2])
    annotation = set(range(st, en+1))
    n_wdws = float(en-st-wdw_size)           # Number of possible windows in that gene

    # Count the number of windows repeated
    repeated = []
    for i in u.list_generator('/home/smiravet/crg/mycorepo/repeated_'+str(wdw_size)+'.txt', 0):
        st, en = [int(pos) for pos in i.split('-')]
        repeated += range(st, en)
    repeated = set(repeated)

    n = len(annotation.intersection(repeated))-wdw_size+1.0

    ratio = n/n_wdws*100.0

    if ratio <= 0:
        ratio = 0.0

    return ratio



############ RBS STUDY FUNCTIONS

def retrieve_genes_DOOR(organism_code, genetype, header=True):
    """
    Given an organism code retrieves all the genetype genes in the organism
    from the dbs/operons directory

    genetypes accepted:
        m = monocistronic
        p = policistronic
        a = all

    The function returns a dictionary with the strucure {gene:[start, stop, strand]}
    Start is the lowest number, if strand == '-' this values corresponds to the stop
    position
    """

    # Load the operon structure from the opr file
    operon_structure = {}
    annotation       = {}
    with open('/home/smiravet/crg/dbs/operons/'+organism_code+'.opr', 'r') as fi:
        for line in fi:
            if header:
                header=False
            else:
                line = line.strip().split()
                annotation[line[2]] = [int(line[3]), int(line[4]), line[5]]
                oper = line[0]
                if oper in operon_structure:
                    operon_structure[oper].append(line[2])
                else:
                    operon_structure[oper] = [line[2]]

    # susbset the type you need 
    if genetype == 'm':
        subset_genes = [v[0] for k, v in operon_structure.iteritems() if len(v) == 1]
    elif genetype == 'p':
        subset_genes = []
        for k, v in operon_structure.iteritems():
            if len(v) > 1:
                for g in v:
                    subset_genes.append(g)
    else:
        subset_genes = []
        for k, v in operon_structure.iteritems():
            for g in v:
                subset_genes.append(g)

    # Return
    return {k:v for k, v in annotation.iteritems() if k in subset_genes}


def nperKb(n, size):

    return (float(n)/size)*1000


def generate_rbs_df(organisms, RBS_seqs, wdw=15):
    """
    Given a list of organisms and a list of motifs to search
    Computes the rate between monocistronic/policistronic/all genes with RBS and without
    Append growth information if available (NA if not)
    And mean MFE for the -20 motif

    Return all information dataframe

    wdw = upstream window you look
    """

    columns   = ['organism', 'genes', 'RBS','ratio', 'mean_MFE', 'GC', 'generation_time', 'growth_type', 'genome_size', 'ORFs', 'NCBI', 'SEPs', 'ORF/Kb', 'NCBI/Kb', 'SEPs/Kb']
    to_append = []

    growth_rates = {k:float(v) for k, v in u.str_dic_generator('/home/smiravet/crg/dbs/operons/growth_rates.txt', 0, 1).iteritems()}

    for organism in organisms:

        # Prepare
        genome   = u.load_genome('/home/smiravet/crg/dbs/operons/'+organism+'.fasta')
        gsize    = len(genome)
        energies = []

        # subset genes if required and count number of genes
        subset_genes = u.gb2annotation('/home/smiravet/crg/dbs/operons/'+organism+'.gb')
        m = len(subset_genes)

        # check if RBS present and compute energy of -20
        n = 0
        for k, v in subset_genes.iteritems():
            if v[-1]=='+':
                prev = genome[v[0]-wdw:v[0]+3]
            else:
                prev = u.reverse_complement(genome[v[1]-3:v[1]+wdw])

            if prev[-3:] in ['ATG', 'GTG', 'TTG'] and len(prev[:-3]) > 1:
                prev = prev[:-3]
                if len(prev) == wdw:
                    mfe = fold(prev[:-3])[1]
                    energies.append(mfe)
                    if any(rbs in prev for rbs in RBS_seqs):
                        n += 1

        # Check if growth info
        if organism in growth_rates:
            gr = growth_rates[organism]
            gt = 'F' if gr < 2.5 else 'S'
        else:
            gr = 'NA'
            gt = 'NA'

        # ORFs infor
        aa_seqs = u.load_multifasta('/home/smiravet/crg/dbs/rbs_orfs/'+organism+'_small_aa.fa')
        orfs    = len([v for k, v in aa_seqs.iteritems() if len(v) >= 19])
        seps    = len([v for k, v in aa_seqs.iteritems() if len(v) <= 100])

        # add info
        to_append.append([organism, m, n, round(float(n)/m*100,2), np.mean(energies), round(GC(genome), 2), gr, gt, gsize, orfs, m, seps, nperKb(orfs, gsize), nperKb(m, gsize), nperKb(seps, gsize)])

    to_append.sort(key=lambda x: x[5])  # Sort by GC
    df = pd.DataFrame(to_append, columns=columns)
    df = df.set_index(df['organism']) # vertical index = organism
    return df



#####################
#  TRANSCRIPTOMICS  #
#####################

def combine_expression_files(organism):
    """Returns a dictionary of experiments shape= {experimentA{+:{pos:reads}, -:{pos:reads}}}"""
    try:
        with open('./DBs/metapiles/'+organism+'.pickle') as f:
            results = pickle.load(f)
    except:
        results = {}
        for fil in glob.glob('/home/smiravet/crg/antprotect/transcriptomics_comparative/piles/'+organism+'_*'):
            if 'intergenic' not in fil:
                sample  = fil.split('/')[-1][:-4]
                exp_pos = u.dic_generator(fil, 0, 2, header=True) 
                exp_neg = u.dic_generator(fil, 0, 1, header=True)
                results[sample] = {'+':exp_pos, '-':exp_neg}
        print results['mpneumoniae_6h_1']['+'][300]
        with open('./DBs/metapiles/'+organism+'.pickle', 'w') as f:
            pickle.dump(results, f)
    return results


def plot_piles(organism, xaxis, strand, all_piles):
    for experiment, piles in all_piles.iteritems():
        expression = [math.log(piles[strand][x]+1, 2) for x in xaxis]
        plt.plot(xaxis, expression, linewidth=3, label=experiment, alpha=0.5)


def alt_startcodon(ide, nt_seqs):
    codons = [nt_seqs[ide][i:i+3] for i in range(0, len(nt_seqs[ide]), 3)]
    return codons


def plot_profile(st, en, strand, mpn_piles):
    # Expression stuff
    xaxis = range(st, en)
    plot_piles('mpneumoniae', xaxis, strand, mpn_piles)

    plt.title('['+str(st) +' - '+ str(en)+',' +strand+']')
    plt.legend(bbox_to_anchor=(1, 0.8))

    # Set labels and ticks
    plt.xlim(min(xaxis), max(xaxis))
    plt.locator_params(nbins=5)
    plt.xlabel('Genome base position')
    plt.ylabel('Expression [log2(CPM)]')
    plt.show()



def plot_ide(ide, anno_all, mpn_piles, nt_seqs, frames, peptide_map, alt_start=False, overlap=False, peptides=False):
    # Expression stuff
    st, en, strand = anno_all[ide]
    xaxis = range(st-30, en+31)
    plot_piles('mpneumoniae', xaxis, strand, mpn_piles)

    # Plot annotation
    if strand == '+':
        plt.axvline(x=st, c='g', label='Original start: '+nt_seqs[ide][:3])
        plt.axvline(x=en, c='r', label='Original stop: ' +nt_seqs[ide][-3:])
    else:
        plt.axvline(x=en, c='g', label='Original start: '+nt_seqs[ide][:3])
        plt.axvline(x=st, c='r', label='Original stop: ' +nt_seqs[ide][-3:])
    plt.title(ide+'  ['+str(st) +' - '+ str(en)+',' +strand+']')

    # Plot alternative starts
    if alt_start:
        started = False
        codons  = alt_startcodon(ide, nt_seqs)
        if strand == '+':
            for your_codon in range(len(codons[1:])):
                if codons[1:][your_codon] in ['ATG', 'TTG', 'GTG']:
                    if not started:
                        plt.axvline(x=st+(your_codon*3), c='g', linestyle='--', label='Alternative Start', alpha=0.7)
                        started = True
                    else:
                        plt.axvline(x=st+(your_codon*3), c='g', linestyle='--', alpha=0.7)
        else:
            for your_codon in range(len(codons[1:])):
                if codons[1:][your_codon] in ['ATG', 'TTG', 'GTG']:
                    if not started:
                        plt.axvline(x=en-(your_codon*3), c='g', linestyle='--', label='Alternative Start', alpha=0.7)
                        started = True
                    else:
                        plt.axvline(x=en-(your_codon*3), c='g', linestyle='--', alpha=0.7)

    # Plot overlap with NCBI genes
    if overlap:
        if frames[ide][-1] != 'NA':
            geneB = frames[ide][-1]
            stB, enB = anno_all[geneB][:2]
            plt.plot([stB, enB],[-0.5, -0.5], c="blue",linewidth=15, alpha=0.3, label='NCBI gene')

    # Plot peptides
    if peptides:
        if ide in peptide_map:
            started=False
            if strand == '+':
                for pept in peptide_map[ide]:
                    a, b = [st+(i*3)-3 for i in pept]
                    if started:
                        plt.plot([a, b],[-0.2, -0.2], c="purple", linewidth=10, alpha=0.5)
                    else:
                        started=True
                        plt.plot([a, b],[-0.2, -0.2], c="purple", linewidth=10, alpha=0.5, label='Peptide')
            else:
                for pept in peptide_map[ide]:
                    a, b = [en-(i*3)+3 for i in pept]
                    if started:
                        plt.plot([a, b],[-0.2, -0.2], c="purple", linewidth=10, alpha=0.5)
                    else:
                        started=True
                        plt.plot([a, b],[-0.2, -0.2], c="purple", linewidth=10, alpha=0.5, label='Peptide')
    # Legend
    plt.legend(bbox_to_anchor=(1, 0.8))

    # Set labels and ticks
    plt.xlim(min(xaxis), max(xaxis))
    plt.locator_params(nbins=5)
    plt.xlabel('Genome base position')
    plt.ylabel('Expression [log2(CPM)]')
    plt.show()


def overlap_ncRNA(size=20):
    """ Return a dictionary with the overlapping between annotation and ncRNAs """

    # Load databases
    nt_seqs, aa_seqs, NCBI, annotation, lengths, frames, contra = organism_info('mpneumoniae', 816394, size)

    # Load ncRNAs
    nc_RNAs_pos = []
    nc_RNAs_neg = []
    with open('./DBs/nc_mpneumoniae.txt', 'r') as fi:
        for line in fi:
            line = line.strip().split()
            if line[0] == '+':
                nc_RNAs_pos += [sorted([int(i) for i in line[2:]])]
            else:
                nc_RNAs_neg += [sorted([int(i) for i in line[2:]])]

    # Match
    results = []
    for k in aa_seqs.keys():
        stA, enA, strand = annotation[k]
        if strand == '+':
            for stB, enB in nc_RNAs_pos:
                if stB <= stA < enA <= enB:
                    results.append(k)
        else:
            for stB, enB in nc_RNAs_neg:
                if stB <= stA < enA <= enB:
                   results.append(k)
    return set(results)


def generate_mean_expression():
    for fil in glob.glob('/home/smiravet/crg/dbs/smprots_DB/expression/*_1.txt'):
        print fil
        try:
            org   = fil.split('/')[-1].split('_1')[0]
            ides  = [str(i) for i in u.list_generator(fil, 0)]
            x     = [float(i) for i in u.list_generator(fil, 2)]
            y     = [float(i) for i in u.list_generator(fil.replace('_1.', '_2.'), 2)]
            means = {i[0]:np.mean(i[1:]) for i in zip(ides, x, y)}

            fi = open(fil, 'r')
            with open('/home/smiravet/crg/dbs/smprots_DB/expression_mean/'+org+'.txt', 'w') as fo:
                for line in fi:
                    line = line.strip().split()
                    fo.write(line[0]+'\t'+line[1]+'\t'+str(means[line[0]])+'\n')
            fi.close()

            # Mean values for intergenic
            interfil = fil.replace('.txt', '_intergenic.txt')
            m1       = [float(i) for i in u.list_generator(interfil, 0)]
            m2       = [float(i) for i in u.list_generator(interfil.replace('_1.', '_2.'), 0)]

            s1       = [float(i) for i in u.list_generator(interfil, 1)]
            s2       = [float(i) for i in u.list_generator(interfil.replace('_1.', '_2.'), 1)]

            with open('/home/smiravet/crg/dbs/smprots_DB/expression_mean/'+org+'_intergenic.txt', 'w') as fo:
                fo.write(str(np.mean([m1, m2]))+'\t'+str(np.mean([s1, s2])))
        except:
            pass


# DEGENERATED PROTEINS

def contiguous(annotation, strand):
    cont_dict = {}
    orderedDict = collections.OrderedDict(sorted(annotation.iteritems(), key=lambda (k,v):(v,k)))
    sts, ens, ids = {}, [], []
    for k, v in orderedDict.iteritems():
        if v[-1]==strand:
            st, en  = sorted(v[0:2])
            sts[k]  = st
            ens.append(en)
            ids.append(k)
    array_ends = np.array(ens)
    for ide, st in sts.iteritems():
        vals = [i for i in st-array_ends if i>0]
        if len(vals) != 0:
            close = vals.index(min(vals))
            pair  = ids[close]

            # Check they are in frame
            mini, maxi = sorted(anno_all[ide][:2])

            cont_dict[ide]=pair
        else:
            pass
    return cont_dict


def find_closest_organism(organism):
    """
    Finds the closest organism based on the conservation 
    """
    possibles = {}
    annotated = u.list_generator('../../dbs/smprots_DB/id_genes/'+organism+'_pairs.txt')
    blast_fil = '/home/smiravet/crg/antprotect/RF_final/metaRF/blaster/'+organism+'_results_processed.out'
    with open(blast_fil) as fi:
        current_target=''
        for line in fi:
            line = line.strip().split()
            target  = line[0]
            if target in annotated:
                evl  = float(line[-2])
                query = line[1]
                if query[:5]!=target[:5] and evl < 1e-8:
                    if target in possibles:
                        possibles[target].append([query[:5], evl])
                    else:
                        possibles[target] = [[query[:5], evl]]
    possibles = {k:sorted(v, key = lambda x: x[1])[0][0] for k, v in possibles.iteritems()}.values()
    return max(set(possibles), key=possibles.count)    # Return the most common element


def predict_org(code):
    fil = glob.glob('../../dbs/smprots_DB/id_genes/'+code+'*_pairs.txt')[0]
    return fil

def load_close_homology(organism):
    close_conservation={}
    closest_organism = find_closest_organism(organism)
    annotatedB = u.list_generator(predict_org(closest_organism), 0)
    blast_fil = '/home/smiravet/crg/antprotect/RF_final/metaRF/blaster/'+organism+'_results_processed.out'

    with open(blast_fil) as fi:
        for line in fi:
            line = line.strip().split()
            target  = line[0]
            evl  = float(line[-2])
            query_c = line[1]
            if query_c in annotatedB:
                if target not in close_conservation:
                    close_conservation[target]=[[query_c, evl]]
                else:
                    close_conservation[target].append([query_c, evl])
    return close_conservation


def sharing_stop(annotation, close_dict, strand):

    # Prepare dictionaries and list of genes
    cont_dict = {}
    if strand=='+':
        filt = {k:sorted(v) for k, v in annotation.iteritems() if v[-1]=='+'}
        ord_ann = collections.OrderedDict(sorted(filt.iteritems(), key=lambda (k,v):(v,k)))  # Sort by the start
    else:
        filt = {k:[max(v[:2]), min(v[:2]), v[-1]] for k, v in annotation.iteritems() if v[-1]=='-'}
        ord_ann = collections.OrderedDict(sorted(filt.iteritems(), key=lambda (k,v):(v,k), reverse=True)) 
    your_genes = ord_ann.keys()

    # Iterate and return the black list

    for i in range(0, len(your_genes)-1):
        ideA = your_genes[i]
        if ideA in close_dict:
            stA, enA = sorted(ord_ann[your_genes[i]])[:2]
            for ideB in your_genes[i+1:]:
                if strand=='+':
                    en = max(ord_ann[ideB][:2])
                    if enA-en<1000 and (stA-en+1)%3:
                        if ideA in cont_dict:
                            cont_dict[ideA].append(ideB)
                        else:
                            cont_dict[ideA]=[ideB]
                else:
                    en = min(ord_ann[ideB][:2])
                    if stA-en<1000 and (stA-en-1)%3:
                        if ideA in cont_dict:
                            cont_dict[ideA].append(ideB)
                        else:
                            cont_dict[ideA]=[ideB]  
    return cont_dict

def pair_degenerated_conservation(organism):
    """
    One output is the black list including any gene that can come from a fragmentation of a 
    bigger one in mycogenitalium
    and other is for the MPN genes if they are longer
    than expected by conservation
    """

    # Load the information
    lengths    = u.str_dic_generator('/home/smiravet/crg/dbs/smprots_DB/aa_lengths/'+organism+'_lengths.txt', 0, 1)    
    annotation = u.load_annotation('/home/smiravet/crg/dbs/smprots_DB/SEPs_annotation/'+organism+'_annotation.txt')    
    annotation = {k:v for k, v in annotation.iteritems() if k in lengths and lengths[k]>=20}
    annotated  = u.list_generator('../../dbs/smprots_DB/id_genes/'+organism+'_pairs.txt')

    # Prepare the homology dictionary
    close_dict = load_close_homology(organism)

    posstop = sharing_stop(annotation, close_dict, '+')
    negstop = sharing_stop(annotation, close_dict, '-')

    # Create the black list
    black_list = []
    for dataset in [posstop, negstop]:
        for ideA, other_ides in dataset.iteritems():
            for ideB in other_ides:
                if ideA in close_dict and ideB in close_dict and ideB not in annotated:
                    setA = set([vals[0] for vals in close_dict[ideA]])
                    setB = set([vals[0] for vals in close_dict[ideB]])
                    if len(setA.intersection(setB)) > 0:
                        black_list.append(ideB)
    return set(black_list)


#### PROBAS

def probability_per_size(size, specific_GC, mode=0):

    probability = 1.0
    start_codons  = ['ATG', 'TTG', 'GTG']
    stop_codons   = ['TAG', 'TAA']
    if mode>0:
        stop_codons += ['TGA']

    start_prob = sum(motif_probability_GC(start_seq, specific_GC) for start_seq in start_codons)
    stop_prob  = sum(motif_probability_GC(stop_seq, specific_GC) for stop_seq in stop_codons)
    any_prob   = 1-stop_prob

    p = start_prob*stop_prob*(any_prob**(size-1))

    return p


def pSEP_per_GC(specific_GC, mode=0):

    summatory = 0.0
    for i in range(1, 101):
        summatory += probability_per_size(i, specific_GC, mode)

    return summatory


#####################
#         RUN       #
#####################

if '__main__' == __name__:

#    organisms  = [x.split('/')[-1].replace('.fasta', '') for x in glob.glob('/home/smiravet/crg/dbs/operons/*.fasta')]
#
#    RBSsimple  = ['AGGA'  , 'AAGG']
#    RBSs       = ['AGGA'  , 'AAGG', 'AGAGG', 'AGGG', 'GGAGGT', 'AGGAGG', 'AGGAGGT']
#    RBSs2      = ['GGAGGT', 'AGGAGG']
#    RBS7       = ['GGA'   , 'GAG'   , 'AGG'  , 'AGGA'  , 'GGAG'  , 'GAGG' , 'AGGAG' , 'GGAGG' , 'AGAAGG', 'AGCAGG','AGGAGG', 'AGTAGG', 'AGGCGG', 'AGGGGG', 'AGGTGG']
#
#    df_srbs = generate_rbs_df(organisms, RBS7, 'all')
#
#    print df_srbs

    # RFClassifier(organism='mpneumoniae', feats_to_exclude=[], folds=0, test_size=0.2, extension='test_0f_mpneu', all_new=True)
    # conservationALL()

    conservationALL(inFiles='/home/smiravet/crg/dbs/smprots_DB/tempss/mhigh_small_aa.fa', directory='/home/smiravet/crg/antprotect/RF_final/metaRF/blaster15/', size=15)
