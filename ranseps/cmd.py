#!/usr/bin/env python

#############################################################
#
# cmd.py
#
# Author : Miravet-Verde, Samuel
# Last updated : 07/19/2018
#
# Command line parser to run RanSEPs.
# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
#############################################################

import sys
import os
import argparse

from run_ranseps import run_ranseps

def run_all():
    """Command line program processor."""

    # Ensure it is a directory
    if options.outDir[-1]!='/':
        options.outDir+='/'

    run_ranseps(options.genome         , options.cds              , options.outDir          ,options.codon_table       , options.min_size , options.species_code,
                options.eval_thr       , options.threads          ,
                options.eval_thr2      , options.align_thr        , options.length_thr      , options.iden_thr         ,
                options.seps_percentage, options.positive_set_size, options.feature_set_size, options.negative_set_size,
                options.test_size      , options.folds            , options.other_database  , options.blastp_db)


#PARSER
parser = argparse.ArgumentParser(description = "RanSEPs provides a framework for bacterial genome re-annotation and novel small proteins (SEPs) detection adjusting the search to different genomic features that govern protein-coding capabilities.")

parser.add_argument('-g', '--genome',
                    dest="genome",
                    action="store",
                    required=True,
                    type=str,
                    help="Genome of reference in fasta or genbank complete format.\nAutodection: .gb, .gbk, .genbank = genbank; any other file considered as fasta.\nCDS argument is not mandatory, own annotation will be used.")

parser.add_argument('-c', '--CDS',
                    dest="cds",
                    action="store",
                    type=str,
                    help="Coding DNA sequences of reference genome in fasta format. If genome in genbank format, annotation of CDS will be replaced with the ones included in this file.")

parser.add_argument('-o', '--output_directory',
                    dest="outDir",
                    action="store",
                    default='./',
                    type=str,
                    help="Output directory where all the results will be stored.")

parser.add_argument('-t', '--codon_table',
                    dest="codon_table",
                    action="store",
                    default=0,
                    type=int,
                    help="Codon translation table. Accepted tables = [0, 4, 11]. 0 will use as START and STOP codons the ones observed for genes in --CDS")

parser.add_argument('-s', '--min_size',
                    dest="min_size",
                    action="store",
                    default=10,
                    type=int,
                    help="Minimum protein size considered (as number of aminoacids).")

parser.add_argument('-sp', '--species_code',
                    dest="species_code",
                    action="store",
                    default=None,
                    type=str,
                    help="Identifier code used to name the putative proteins. By default RanSEPs uses the 5 first character of the genome name. Ex: mpneumoniae.fa > mpneu")

# To run blast:
parser.add_argument('-beval', '--blast_evalue_threshold',
                    dest="eval_thr",
                    action="store",
                    default=1e-8,
                    type=float,
                    help="e-value threshold stablished to consider a hit as conserved in BlastP")

parser.add_argument('-bthreads', '--blast_threads',
                    dest="threads",
                    action="store",
                    default=12,
                    type=int,
                    help="Number of threads to run BlastP")

# To run conservation analysis and define the negative set:
parser.add_argument('-ceval', '--negative_evalue_threshold',
                    dest="eval_thr2",
                    action="store",
                    default=2e-8,
                    type=float,
                    help="e-value threshold stablished to consider a hit as conserved for the negative set definition")

parser.add_argument('-calign', '--negative_alignment_threshold',
                    dest="align_thr",
                    action="store",
                    default=0.0,
                    type=float,
                    help="Alignment threshold stablished to consider a hit as conserved for the negative set definition")

parser.add_argument('-clength', '--negative_length_threshold',
                    dest="length_thr",
                    action="store",
                    default=0.0,
                    type=float,
                    help="Length threshold stablished to consider a hit as conserved for the negative set definition")

parser.add_argument('-cident', '--negative_identity_threshold',
                    dest="iden_thr",
                    action="store",
                    default=50.0,
                    type=float,
                    help="Identity threshold stablished to consider a hit as conserved for the negative set definition")

# RanSEPs arguments

parser.add_argument('-rp', '--seps_percentage',
                    dest="seps_percentage",
                    action="store",
                    default=0.25,
                    type=float,
                    help="Fraction (base 1) of SEPs included in the positive and feature training sets")

parser.add_argument('-pn', '--positive_set_size',
                    dest="positive_set_size",
                    action="store",
                    default=100,
                    type=int,
                    help="Training positive set size for each iteration")

parser.add_argument('-fn', '--feature_set_size',
                    dest="feature_set_size",
                    action="store",
                    default=100,
                    type=int,
                    help="Feature set size for each iteration")

parser.add_argument('-nn', '--negative_set_size',
                    dest="negative_set_size",
                    action="store",
                    default=150,
                    type=int,
                    help="Training negative set size for each iteration")

parser.add_argument('-tn', '--test_size',
                    dest="test_size",
                    action="store",
                    default=0.2,
                    type=float,
                    help="Percentage (base 1) of sequences from the training sets used in testing")

parser.add_argument('-f', '--folds',
                    dest="folds",
                    action="store",
                    default=50,
                    type=int,
                    help="Training positive set size for each iteration")

parser.add_argument('-bp', '--blastp_db',
                    dest='blastp_db',
                    action="store",
                    default=None,
                    type=str,
                    help=" Uses the provided file as custom amino acidic database for the BlastP step. This allows custom searches against user-defined datasets. No conserved sequences with this file will be used in the negative training set.")

parser.add_argument('-db', '--database',
                    dest='other_database',
                    action="store",
                    default=None,
                    type=str,
                    help="Uses files from a previous prediction recomputing the prediction. Just add here the location of 'intermediary_files' and the ORFs and conservation will be loaded here")

options = parser.parse_args()

if __name__=="__main__":
    run_all()

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

