#!/usr/bin/env python

import sys
import os
import argparse

from run_ranseps import run_ranseps

def run_all():
    """Command line program processor."""

    # Ensure it is a directory
    if options.outDir[-1]!='/':
        options.outDir+='/'

    run_ranseps(options.genome      , options.cds      , options.outDir      ,
                options.codon_table , options.min_size , options.species_code,
                options.eval_thr    , options.threads)


#PARSER
parser = argparse.ArgumentParser(description = "RanSEPs provides a framework for genome re-annotation and novel small proteins detection adjusting the search to different genomic features that govern protein-coding capabilities.")

parser.add_argument('-g', '--genome',
                    dest="genome",
                    action="store",
                    required=True,
                    type=str,
                    help="Genome of reference in fasta format.")

parser.add_argument('-c', '--CDS',
                    dest="cds",
                    action="store",
                    required=True,
                    type=str,
                    help="Coding DNA sequences of reference genome in fasta format.")

parser.add_argument('-o', '--output_directory',
                    dest="outDir",
                    action="store",
                    default='./',
                    type=str,
                    help="Output directory where all the results will be stored.")

parser.add_argument('-t', '--codon_table',
                    dest="codon_table",
                    action="store",
                    default=11,
                    type=int,
                    help="Codon translation table. Accepted tables = [0, 4, 11]. 0 will use as START and STOP codons the ones observed for genes in --CDS")

parser.add_argument('-s', '--min_size',
                    dest="min_size",
                    action="store",
                    default=10,
                    type=int,
                    help="Minimum protein size considered.")

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

# To run blast:


options = parser.parse_args()

if __name__=="__main__":
    run_all()

