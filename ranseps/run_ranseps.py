#!/usr/bin/env python

#############################################################
#
# run_ranseps.py
#
# Author : Miravet-Verde, Samuel
# Last updated : 11/15/2017
#
#############################################################

import sys, os
import seps_functions as sf
import utils as u

from orfinder import run_orfinder
from generate_dbs import run_gdbs
from blaster import run_blaster

def run_ranseps(genome, cds, outDir, codon_table, min_size, species_code, eval_thr, threads):

    # Create a folder for intermediary files
    intDir = outDir+'intermediary_files/'
    if not os.path.exists(intDir):
        os.makedirs(intDir)

    # Define code
    if species_code==None:
        species_code = genome.split('/')[-1].split('.')[0][:5]

    # Genome size:
    gl = len(u.load_multifasta(genome))

    # Generate DBs
    run_orfinder(genome_handle=genome, cds=cds, outdir=intDir, ct=codon_table, min_size=min_size, species_code=species_code)
    run_gdbs(genome=genome, cds=cds, outdir=intDir, min_size=min_size, species_code=species_code, genome_length=gl)

    # Run Blast
    run_blaster(outdir=intDir, min_size=min_size, species_code=species_code, threshold=eval_thr, threads=threads)




    # Run RanSEPs



