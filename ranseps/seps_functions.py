#!/usr/bin/env python

#############################################################
#
# seps_functions.py
#
# Author : Miravet-Verde, Samuel
# Written : 03/06/2017
# Last updated : 07/19/2018
#
# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
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
import numpy as np
import pandas as pd
import ranseps_utils as u
import matplotlib.pyplot as plt

# To compute specificic features
from propy.QuasiSequenceOrder import GetSequenceOrderCouplingNumberSW as SOSW
from propy.CTD import CalculateDistributionHydrophobicity as DHF
from propy.CTD import CalculateDistributionSecondaryStr as DSS

# Others
from scipy import interp
from operator import mul
from collections import Counter, OrderedDict
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve

# Bio handling
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex

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
#  BASIC FUNCTIONS  #
#####################

def filter_dic(dic, exclude):
    for ide in exclude:
        dic.pop(ide, None)
    return dic

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


##########################
# CONSERVATION FUNCTIONS #
##########################

def load_all_orthologs(subcode, evalth=2e-8, alignth=50.0, lenth=58.0, identh=50.0, close=[], size=19):
    """
    Counts the number of times a protein is conserved based on several threshold
    filters:
        - evalue
        - alignment length (percentage of the protein that is properly aligned)
        - % of the length difference between the two annotations
        - % identity

    Close allows to discard close species in the counting
    """
    # Tajke files, measure their lengths and define some variables
    targets_fa = u.load_multifasta(subcode+'_small_aa.fa')
    dbquery_fa = u.load_multifasta(subcode+'_db_filtered.fa')

    lens_targets = {ide:len(seq) for ide, seq in targets_fa.iteritems()}
    lens_dbquery = {ide:len(seq) for ide, seq in dbquery_fa.iteritems()}

    blast_outs = subcode+'_results_processed.out'

    close += ['mferi']

    results = {k:[] for k in targets_fa.keys()}
    with open(blast_outs) as fi:
        for line in fi:
            line = line.strip().split()
            target  = line[0]
            query = line[1][:5]
            query_c = line[1]

            if query not in close and target[:5] != query:
                evl  = float(line[-2])            # eval
                ale  = int(line[3])-int(line[5])  # alignment length - gap opens
                idp  = float(line[2])

                ale = 100.0*float(ale)/lens_targets[target]
                ple = 100.0*(max(lens_targets[target],lens_dbquery[query_c])-abs(lens_targets[target]-lens_dbquery[query_c]))/max(lens_targets[target],lens_dbquery[query_c])

                if evl <= evalth and ale >= alignth and ple >= lenth and idp >= identh:
                    results[target] += [query_c]

    return {k:set(v) for k, v in results.iteritems()}


def load_number_times_conserved(subcode, evalth, alignth, lenth, identh, close, size):
    """
    Counts the number of times a protein is conserved based on several threshold
    filters:
        - evalue
        - alignment length (percentage of the protein that is properly aligned)
        - % of the length difference between the two annotations
        - % identity
    Close allows to discard close species in the counting
    """
    pre_results = load_all_orthologs(subcode, evalth, alignth, lenth, identh, close, size)
    results = {}
    for ide, orthologs in pre_results.iteritems():
        results[ide] = len(set([v[:5] for v in orthologs]))

    return results


#####################
# FRAMES  FUNCTIONS #
#####################

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

#####################
# DATABASE HANDLING #
#####################

def organism_info(organism, genome_size, size):

    your_NCBI = u.str_dic_generator(organism+'_pairs.txt', 0, 2, header=False, split_by='\t')
    # Load sequences
    your_nt_seqs = u.load_multifasta(organism+'_small_nt.fa')
    your_aa_seqs = u.load_multifasta(organism+'_small_aa.fa')
    your_lengths = {k:len(v) for k, v in your_aa_seqs.iteritems()}
    assert len(set(your_NCBI).difference(set(your_aa_seqs)))==0
    # Load annotation & correspondance
    your_annotation =  u.load_annotation(organism+'_annotation.txt')
    # Frames and length
    your_frames, your_contra = add_frames(your_annotation, your_NCBI, genome_size)
    return your_nt_seqs, your_aa_seqs, your_NCBI, your_annotation, your_lengths, your_frames, your_contra


#####################
# FEATURE FUNCTIONS #
#####################

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



def run_m10p20(genome, training_ides, testing_ides, annotation):
    """
    Look for motives in the -10 +20 region of the SEP
    If random == True --> testing ides has to be the sequences!
    """
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

    # 8. create the pwm and the pssm
    background = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
    m = motifs.create(minus10)
    pwm10 = m.counts.normalize(pseudocounts=background)  # Pseudocounts based on the GC content (40%)
    pssm10 = pwm10.log_odds(background)

    m = motifs.create(plus20)
    pwm20 = m.counts.normalize(pseudocounts=background)
    pssm20 = pwm20.log_odds(background)

    # 9. The score for a sequences is computed as the product of probabilities per base per position
    def motifscore(your_seq, your_pssm):
        """ Basic function to return the product of probabilities for a sequence given a motif """
        # to know the prob of having an A in the 2nd position --> pssm20['A',1] // using this:
        return sum([your_pssm[char, index] for char, index in zip(list(your_seq), range(0, len(your_seq)))])

    # 9.2. If random sequences passed by a dictionary:
    results = {}
    for ide in testing_ides:
        if annotation[ide][2]=='+':
            seq = genome[annotation[ide][0]-11:annotation[ide][0]+22]
        else:
            seq = str(Seq(genome[annotation[ide][1]-23:annotation[ide][1]+10]).reverse_complement())
        if check_nt(seq):
            results[ide] = [motifscore(seq[0:10], pssm10), motifscore(seq[13:], pssm20)]
        else:
            results[ide] = [0.0, 0.0]
    return results


def compute_m10p20_dict(genome, nt_seqs, feature_set, annotation):
    orfs_m10p20_sc = run_m10p20(genome, feature_set, nt_seqs.keys(), annotation)
    minus10 = {k:float(v[0]) for k, v in orfs_m10p20_sc.iteritems()}
    plus20  = {k:float(v[1]) for k, v in orfs_m10p20_sc.iteritems()}
    return minus10, plus20


def stack_energy(sequence):
    """Given a sequence stringReturns stacking energy at 37C (AG37)"""
    energies = {'AA':-1.00, 'AC':-1.44, 'AG':-1.28, 'AT':-0.88,
                'CA':-1.45, 'CC':-1.84, 'CG':-2.17, 'CT':-1.28,
                'GA':-1.30, 'GC':-2.24, 'GG':-1.84, 'GT':-1.44,
                'TA':-0.58, 'TC':-1.30, 'TG':-1.45, 'TT':-1.00}
    duplets = [sequence[i:i+2] for i in range(0, len(sequence)-1, 1)]
    stack_energies = [energies[duplet] if duplet in energies else 0.0 for duplet in duplets] # For cases like ecoli where you have R and N in the sequence
    if len(sequence) > 0:
        return (sum(stack_energies))/len(sequence)
    else:
        return 0.0


def hasRBS(seq):
    """ Return True if sequence has an RBS """
    RBS_seqs = ['GGA'   , 'GAG'   , 'AGG'   ,
                'AGGA'  , 'GGAG'  , 'GAGG'  ,
                'AGGAG' , 'GGAGG' ,
                'AGAAGG', 'AGCAGG', 'AGGAGG', 'AGTAGG',
                'AGGCGG', 'AGGGGG', 'AGGTGG']
    return 1.0 if any(RBS_seq in seq for RBS_seq in RBS_seqs) else 0.0


def generate_DB_stack_energy(genome, annotation):
    """
    Generates a DB with the values for stacking energy,
    If annotation is provided that dictionary is used instead of the general
    annotation of the organism is not used
    """
    # Load genome and annotation if required
    lg      = len(genome)
    genomec = genome+genome # circular

    stack_energies_dic = {}
    for k, v in annotation.iteritems():
        if v[-1] == '+':
            st = v[0]
            if st > 16:
                seq_RBS    = genomec[st-16:st-1]
                seq_m10p20 = genomec[st-11:st+19]
            else:
                seq_RBS    = genomec[st+lg-16:st+lg-1]
                seq_m10p20 = genomec[st+lg-11:st+lg+19]
        else:
            st = max(v[0:2])
            if st < lg-16:
                seq_RBS    = genomec[st:st+15]
                seq_m10p20 = genomec[st-23:st+10]
            else:
                seq_RBS    = genomec[st+lg:st+lg+15]
                seq_m10p20 = genomec[st+lg-23:st+lg+10]
        # Add the result
        sE1 = stack_energy(seq_RBS)
        RBS = hasRBS(seq_RBS)
        sE2 = stack_energy(seq_m10p20)
        stack_energies_dic[k] = [sE1, sE2, RBS]

    return stack_energies_dic


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
    return [sum([freq_dict[i][hexamer] if check_nt(hexamer) else 0.0 for hexamer in lol[i]])/float(len(lol[i])) for i in [0, 1, 2, 3]]


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

def create_directory_structure(project_folder):
    """ Given a project directory creates a subtree directory """
    if not os.path.exists(project_folder):
        os.mkdir(project_folder)
        for subfolder in ['classification_stats', 'classification_tasks', 'classifiers', 'feature_importances', 'features', 'results']:
            if not os.path.exists(project_folder+'/'+subfolder):
                os.mkdir(project_folder+'/'+subfolder)


def short3propy(aa_seqs):
    """ Short version to compute only the 3 most important features in propy """

    propy_nms = ['ide', 'tausw9', '_HydrophobicityD1001', '_SecondaryStrD1001']
    to_append = []
    for ide, seq in aa_seqs.iteritems():
        l = len(seq)
        to_append.append([ ide, SOSW(seq)['tausw9']/l, DHF(seq)['_HydrophobicityD1001']/l, DSS(seq)['_SecondaryStrD1001']/l])

    # Create file
    df = pd.DataFrame(to_append, columns=propy_nms)
    propy_df = df.set_index(df['ide']) # vertical index = organism
    propy_df.drop('ide', axis=1, inplace=True)
    return propy_df


def featurizer2(genome, organism, nt_seqs, aa_seqs, annotation, positive_set, feature_set, negative_set, extension, project_folder):
    """
    This function requires the user to provide two dictionaries of nt and aa seqs with ides and a list of
    ids that will be used to generate the positive, negative and feature set.
    Organism is required to load the genome and compute some features
    Annotation refers to a dict with ide:[st, en, strand]
    Extension is used to store the feature DBs
    """

    # Generate a fasta file with the sequences in the feature set
    with open(project_folder+'features/feat_'+organism+'_'+extension+'_sequences.fa', 'w') as fo:
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
    cai.generate_index(project_folder+'features/feat_'+organism+'_'+extension+'_sequences.fa')

    # 2 aa N and C terminal
    nt_diaa_freq, ct_diaa_freq = NCterminal([aa_seqs[ide] for ide in feature_set])
    nns = {ide:nt_diaa_freq[seq[1:3]] for ide, seq in aa_seqs.iteritems()}
    ccs = {ide:ct_diaa_freq[seq[-2:]] for ide, seq in aa_seqs.iteritems()}

    # -10 +20 motifs and stack energies
    minus10, plus20 = compute_m10p20_dict(genome, nt_seqs, feature_set, annotation)
    stackEEs        = generate_DB_stack_energy(genome, annotation)

    # Hexamers
    # define the background model
    alphabets         = ['A', 'C', 'G', 'T']
    GC_content        = float(np.mean([GCs[ide] for ide in set(positive_set).union(set(feature_set))]))/100.0
    proportions       = {'A':(1-GC_content)/2, 'C':GC_content/2, 'G':GC_content/2, 'T':(1-GC_content)/2}
    background_model  = set([''.join(x) for x in itertools.product(alphabets, repeat = 6)])
    background_model  = {cd:reduce(lambda x, y: x*y, [proportions[x] for x in cd]) for cd in background_model}
    assert 0.98 < sum(background_model.values()) <1.01

    hex_freq_dict = hexamer_freqs([str(x) for x in u.load_multifasta(project_folder+'features/feat_'+organism+'_'+extension+'_sequences.fa').values()], background_model)

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
        your_hex0     = hex0[ide] # dicodon freq
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
    propy_db = short3propy(aa_seqs)
    if not propy_db.index.is_unique:
        propy_db.drop_duplicates(inplace=True)

    # Concate and remove possible Nan (propy has everything!)
    feats_df = pd.concat([feats_df.iloc[:,:-1], propy_db, feats_df.iloc[:,-1:]], axis=1)
    feats_df = feats_df.dropna()

    # Store the pickle object
    feats_df.to_pickle(project_folder+'features/feat_'+organism+'_'+extension+'.pickle')

    return feats_df


def miniRanSeps(genome         , organism        , nt_seqs          , aa_seqs        , annotation,
                positive_set   , feature_set     , negative_set     , to_exclude     ,
                sfolds         , test_size       , random_state_test, n_estimators   ,
                oob_score      , n_jobs          , random_state     , max_depth      ,
                max_features   , min_samples_leaf, extension        , project_folder):
    """
    Run a whole RF analysis over the organism given, same idea than the previous version but
    this is thought to be iterative ;)

    sfolds represents how many times you repeat the classification test OVER THE SAME DATASET, not required as the
    kfolds is done by the iteration of this function.
    """

    features_df = featurizer2(genome, organism, nt_seqs, aa_seqs, annotation, positive_set, feature_set, negative_set, extension, project_folder)

    # Remove features we don't want in the analysis
    if len(to_exclude) > 0:
        for te in to_exclude:
            features_df = features_df.drop(te, 1)

    # Extract X and y (features and labels)
    k1 = features_df.loc[(features_df.set_type == '+') | (features_df.set_type == '+f') | (features_df.set_type == '-')]
    features  = np.array(k1.iloc[:,1:-1])  # Extract an array without ide and set_type info
    labels  = np.array([1 if i in ['+', '+f'] else 0 for i in k1.iloc[:,-1]])    # Labels

    # Names of the features
    feat_names = [str(m) for m in k1.iloc[:,1:-1].columns]
    numbers    = range(0,len(feat_names))
    feat_dict  = dict(zip(numbers, feat_names))

    # Define the classifier
    clf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, oob_score=oob_score, n_jobs=n_jobs,
                                 random_state=random_state, max_features=max_features, min_samples_leaf=min_samples_leaf)
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
            fpr, tpr, thresholds             = roc_curve(y_test, probs[:,1])
            mean_tpr                        += interp(mean_fpr, fpr, tpr)

            roc_auc                          = auc(fpr, tpr)                        # If you need the value of AUC for that classification
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
        fpr_tpr, recall_precision, mean_tpr, mean_auc = None, None, None, None

    # Create the classifier (will generate the figures in parallel)
    clf.fit(features, labels)
    with open(project_folder+'classifiers/clf_'+organism+'_'+extension+'.pickle', 'w') as f:
        pickle.dump(clf, f)

    # FEATURE IMPORTANCES
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
    your_ides      = np.array(features_df.iloc[:,0])      # Only ides
    to_classify    = np.array(features_df.iloc[:,1:-1])   # Remove assigned class and the column of ides

    # Predict probs and classes
    class_pred     = clf.predict(to_classify)
    probs_pred     = clf.predict_proba(to_classify)

    clf_results    = np.column_stack((your_ides, probs_pred, class_pred))
    cols           = ['ide', 'prob_0', 'prob_1', 'pred_class']
    clf_results_df = pd.DataFrame(clf_results, columns=cols)

    # Store everything in a csv and return the results
    clf_results_df.to_csv(project_folder+'classification_tasks/'+organism+'_'+extension+'_classified.csv', sep='\t', index=False)
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
        od = OrderedDict(sorted(probess.items()))
        for k, v in od.iteritems():
            fo.write(k+'\t'+'\t'.join([str(val) for val in v])+'\n')


def RanSEPs(genome          , organism         , nt_seqs        , aa_seqs        , annotation,
            autoset         , set_sizes        , positive_set   , feature_set    ,
            negative_set    , to_exclude       , folds          , sfolds         ,
            test_size       , random_state_test, n_estimators   , oob_score      ,
            n_jobs          , random_state     , max_depth      , max_features   ,
            min_samples_leaf, extension        , project_folder):
    """
    Iteratively run miniRanSeps & extract results
    autoset allows to generate the positive, negative and feature set automatically with sizes included in set sizes.
    autoset[0] = Value of autoset defines the percentage of small proteins in the positive and feature set.
    autoset[1] = list of identifiers known as annotated
    autoset[2] = list of identifiers suitable to be NEGATIVE (ex: no conserved and no overlapping
    autoset[3] = list of additional negative set in case you want to do fancy stuff combining set sizes (small big for example
    """

    # Create the hierarchical tree of directories to write the results
    if project_folder[-1]!='/':
        create_directory_structure(project_folder)   # os.mkdir needs the path with no '/' at the end
        project_folder += '/'
    else:
        create_directory_structure(project_folder[:-1])
    project_folder = project_folder.replace('//', '/')

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
    print 'Training a classifiers'
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
        clf, clf_results_df, fpr_tpr, recall_precision = miniRanSeps(genome        , organism        , nt_seqs               , aa_seqs      , annotation,
                                                                     positive_set  , feature_set     , negative_set          , to_exclude   ,
                                                                     sfolds        , test_size       , random_state_test     , n_estimators ,
                                                                     oob_score     , n_jobs          , random_state          , max_depth    ,
                                                                     max_features  , min_samples_leaf, extension+'_'+str(i+1),project_folder)
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
        # rocax.legend(loc="lower right")
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

##########################
#HANDLING OUTPUT FUNCTION#
##########################

def handle_outputs(species_code, outDir):

    # Move figures:
    c = 'cp '+outDir+'intermediary_files/rs_results/results/*.svg '+outDir
    os.system(c)

    # Generate the last file
    header = ['identifier', 'start'   , 'stop', 'strand', 'aa_length',
              'annotation', 'function', 'RanSEPs_score' , 'RanSEPs_stdv',
              'nt_seq'    , 'aa_seq']

    intDir = outDir+'intermediary_files/'+species_code

    with open(outDir+species_code+'_predicted.csv', 'w') as fo:
        fo.write('\t'.join(header)+'\n')
        dbaaseqs = u.load_multifasta(intDir+'_small_aa.fa')
        dbntseqs = u.load_multifasta(intDir+'_small_nt.fa')
        dbannots = u.load_annotation(intDir+'_annotation.txt')
        dbpaires = u.str_dic_generator(intDir+'_pairs.txt', 0, 1)
        dbfuncts = u.str_dic_generator(intDir+'_pairs.txt', 0, 2, split_by='\t')
        dbscores = u.str_dic_generator(outDir+'intermediary_files/rs_results/results/'+species_code+'_final_prediction.txt', 0, 1)
        dbstdevs = u.str_dic_generator(outDir+'intermediary_files/rs_results/results/'+species_code+'_final_prediction.txt', 0, 2)

        ides = sorted(dbaaseqs.keys())
        for ide in ides:
            if ide in dbpaires:
                p = dbpaires[ide].replace(';', ',')
                f = dbfuncts[ide].replace(';', ',')
            else:
                p = 'Putative ORF'
                f = 'Unassigned function'
            if ide not in dbscores:
                dbscores[ide] = 'Not a valid sequence'
                dbstdevs[ide] = 'Not a valid sequence'
            to_write = [ide]+[str(x) for x in dbannots[ide]]+[str(len(dbaaseqs[ide])), p, f, str(dbscores[ide]), str(dbstdevs[ide]), dbntseqs[ide], dbaaseqs[ide]]
            fo.write('\t'.join(to_write)+'\n')
        fo.close()

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
