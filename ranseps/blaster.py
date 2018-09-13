#!/usr/bin/env python

#############################################################
#
# blaster.py
#
# Author : Miravet-Verde, Samuel
# Written : 03/06/2017
# Last updated : 08/21/2018
# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
#############################################################

#####################
#      IMPORTS      #
#####################

import os
import sys
import glob
import numpy as np
import ranseps_utils as u
from Bio.Seq import Seq
from collections import Counter

import pkg_resources

def build_all_fa(outdir, blastp_db):

    if blastp_db:
        cmd = 'cp '+blastp_db+' '+outdir+'all.fa'
        os.system(cmd)
    else:
        resource_package = __name__  # Could be any module/package name
        resource_path = '/'.join(('dbs', '*.fa'))
        template = pkg_resources.resource_filename(resource_package, resource_path)
        cmd = 'cat '+template+' > '+outdir+'all.fa'
        os.system(cmd)


def is_fragment(ide, values, reference, stdev_lens, all_lens):
    fragment = False
    lA = all_lens[ide]
    for sseqid, pident, length, qstart, qend, sstart, send, evalue in values:
        lB = qend - qstart + 1
        if lB>lA*0.85:
            if (any(ref in sseqid for ref in reference)):
                lC = all_lens[sseqid]
                if sstart<20 and abs(lC-send)>stdev_lens:
                    fragment = sseqid
                elif sstart>stdev_lens and send>lC-20:
                    fragment = sseqid
    return fragment


def function_known(ide, ides, all_functions):
    markers = ['inference', 'hypothetical', 'ensembl', 'locus',
               'putative', 'uncharacterized', 'homology']
    vector = {}
    for i in ides:
        if i in all_functions and i[:5]!=ide[:5]:
            f = all_functions[i].lower()
            if any(marker in f for marker in markers):
                if 2 in vector:
                    vector[2].append(i+'||'+all_functions[i])
                else:
                    vector[2] = [i+'||'+all_functions[i]]
            else:
                if 1 in vector:
                    vector[1].append(i+'||'+all_functions[i])
                else:
                    vector[1] = [i+'||'+all_functions[i]]
    if 1 in vector:
        return [1, vector[1][0]]
    elif 2 in vector:
        return [2, vector[2][0]]
    else:
        return [3, 0]


def run_blaster(outdir, min_size, species_code, threshold, threads, blastp_db):
    """
    minus allows to subset the sequences with length of aa minus <= length <= plus
    """

    dir_handle   = outdir+species_code
    your_aa_seqs = outdir+species_code+'_small_aa.fa'
    build_all_fa(outdir, blastp_db) # Build DB

    # Reformat the db file to ensure each sequence is just a line satisfying the size
    ides = []
    all_seqs = {}
    for ide, seq in u.load_multifasta(outdir+'all.fa').iteritems():
        if min_size <= len(seq):
            all_seqs[ide] = seq
    for ide, seq in u.load_multifasta(your_aa_seqs).iteritems():
        ides.append(ide)
        all_seqs[ide] = seq
    for ide, seq in u.load_multifasta(outdir+'ncbi_aa.fa').iteritems():
        seq = 'M'+str(Seq(seq).translate())[1:-1].replace('*', 'W')
        all_seqs[ide] = seq
    all_lens = {}
    with open(dir_handle+'_db_filtered.fa', 'w') as fo:
        for ide, seq in all_seqs.iteritems():
            all_lens[ide] = len(seq)
            fo.write('>'+ide+'\n'+seq+'\n')

    # All functional seps in the db
    resource_package = __name__
    resource_path = '/'.join(('dbs', 'all_pairs.txt'))
    template = pkg_resources.resource_filename(resource_package, resource_path)
    all_functions = {k:v.lower() for k, v in u.str_dic_generator(template, 0,2, split_by='\t').iteritems()}
    annotated_in_db = set(all_functions.keys())

    # Generate the database required to run the blas analysis
    cmd = 'makeblastdb -in '+dir_handle+'_db_filtered.fa -dbtype prot -parse_seqids'
    os.system(cmd)

    # Run blastp
    print('running blast...')
    cmd = 'blastp -db '+dir_handle+'_db_filtered.fa -query '+your_aa_seqs+' -out '+dir_handle+'_results_processed.out -outfmt 6 -evalue '+str(threshold)+' -num_threads '+str(threads)
    os.system(cmd)

    # Load blast fil and correct the pairs
    print ('generating last requirement files')
    results = {ide:[] for ide in ides}
    NCBI_pairs = u.str_dic_generator(dir_handle+'_pairs.txt', 1, 0, header=False, split_by='\t')
    NCBI_ides  = set(NCBI_pairs.values())
    print len(NCBI_ides)
    close = []
    with open(dir_handle+'_results_processed.out') as fi:
        for line in fi:
            line   = line.strip().split()
            qseqid = line[0]
            sseqid = line[1]
            pident = float(line[2])
            length = int(line[3])
            qstart, qend, sstart, send = [int(i) for i in line[6:10]]
            evalue = float(line[10])
            # To define closely-related organisms
            if sseqid.startswith('NCBIRANSEPS'):
                if abs(all_lens[qseqid]-all_lens[sseqid])<=20 and pident>=0.85 and sseqid not in NCBI_ides:
                    NCBI_pairs[sseqid]=qseqid
            else:
                if qseqid!=sseqid:
                    hit = [sseqid, pident, length, qstart, qend, sstart, send, evalue]
                    results[qseqid].append(hit)
                    if qseqid in NCBI_ides and sseqid in annotated_in_db:
                        close.append(sseqid[:5])
    NCBI_pairs = {v:k for k, v in NCBI_pairs.iteritems()}
    NCBI_ides  = set(NCBI_pairs.keys())
    print len(NCBI_ides)

    # Define close
    thr = 0.85*len(NCBI_ides)
    close = set({k:v for k, v in Counter(close).iteritems() if v>=thr}.keys())

    # Define homology types used in priorizitation
    stdev_lens = np.std(all_lens.values())

    print ('assigning homology types...')
    fo = open(dir_handle+'_homology_types.out', 'w')
    for k, vals in results.iteritems():
        if len(vals)==0:
            typ, ff = 0, 0
        else:
            lll = set([v[0][:5] for v in vals])
            if len(lll.difference(close))>0:
                typ, ff = 3, 0
            else:
                typ, ff = 0, 0
        if k in NCBI_ides:
            typ, ff = 6, NCBI_pairs[k]
        else:
            ides = [v[0] for v in vals if v[:5]!=species_code[:5]]
            typ, ff = function_known(k, ides, all_functions)
            if typ==3 or typ==0:
                codes = [v[0][:5] for v in vals]
                fragment = is_fragment(k, vals, close, stdev_lens, all_lens)
                if codes.count(species_code[:5])>=4:
                    typ = 5
                elif fragment:
                    typ = 4
                    ff = fragment
        fo.write(k+'\t'+str(typ)+'\t'+str(ff)+'\n')
    fo.close()

    print 'close organisms codes:', close, '\n------\n'
    return list(close)

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
