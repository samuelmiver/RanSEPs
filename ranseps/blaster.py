#!/usr/bin/env python

#############################################################
#
# blaster.py
#
# Author : Miravet-Verde, Samuel
# Written : 03/06/2017
# Last updated : 11/16/2017
#
#############################################################

#####################
#      IMPORTS      #
#####################

import os
import sys
import glob
import utils as u

import pkg_resources



def run_blaster(outdir, min_size, species_code, threshold, threads):
    """
    minus allows to subset the sequences with length of aa minus <= length <= plus
    """

    dir_handle   = outdir+species_code
    your_aa_seqs = outdir+species_code+'_small_aa.fa'

    resource_package = __name__  # Could be any module/package name
    resource_path = '/'.join(('dbs', 'all.fa'))
    template = pkg_resources.resource_filename(resource_package, resource_path)

    # Reformat the db file to ensure each sequence is just a line satisfying the size
    fo = open(outdir+'db_filtered.fa', 'w')
    for ide, seq in u.load_multifasta(template).iteritems():
        if min_size <= len(seq):
            fo.write('>'+ide+'\n'+seq+'\n')
    fo.close()

    # Generate the database required to run the blas analysis
    cmd = 'makeblastdb -in '+outdir+'db_filtered.fa -dbtype prot -parse_seqids'
    os.system(cmd)

    # Run blastp
    print('running blast...')
    cmd = 'blastp -db '+outdir+'db_filtered.fa -query '+your_aa_seqs+' -out '+dir_handle+'_results_processed.out -outfmt 6 -evalue '+str(threshold)+' -num_threads '+str(threads)
    os.system(cmd)
