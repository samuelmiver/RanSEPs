#!/usr/bin/env python

#############################################################
#
# blaster.py
#
# Author : Miravet-Verde, Samuel
# Written : 03/06/2017
# Last updated : 07/19/2018
# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
#############################################################

#####################
#      IMPORTS      #
#####################

import os
import sys
import glob
import ranseps_utils as u

import pkg_resources

def build_all_fa(outdir):

    resource_package = __name__  # Could be any module/package name
    resource_path = '/'.join(('dbs', '*.fa'))
    template = pkg_resources.resource_filename(resource_package, resource_path)
    cmd = 'cat '+template+' > '+outdir+'all.fa'
    os.system(cmd)


def run_blaster(outdir, min_size, species_code, threshold, threads):
    """
    minus allows to subset the sequences with length of aa minus <= length <= plus
    """

    dir_handle   = outdir+species_code
    your_aa_seqs = outdir+species_code+'_small_aa.fa'
    build_all_fa(outdir) # Build DB

    # Reformat the db file to ensure each sequence is just a line satisfying the size
    fo = open(dir_handle+'_db_filtered.fa', 'w')
    for ide, seq in u.load_multifasta(outdir+'all.fa').iteritems():
        if min_size <= len(seq):
            fo.write('>'+ide+'\n'+seq+'\n')
    fo.close()

    # Generate the database required to run the blas analysis
    cmd = 'makeblastdb -in '+dir_handle+'_db_filtered.fa -dbtype prot -parse_seqids'
    os.system(cmd)

    # Run blastp
    print('running blast...')
    cmd = 'blastp -db '+dir_handle+'_db_filtered.fa -query '+your_aa_seqs+' -out '+dir_handle+'_results_processed.out -outfmt 6 -evalue '+str(threshold)+' -num_threads '+str(threads)
    os.system(cmd)

    # Define close organism
    NCBI_ides = u.list_generator(dir_handle+'_pairs.txt', 0)
    cons = []
    with open(dir_handle+'_results_processed.out') as fi:
        for line in fi:
            line   = line.strip().split()
            target = line[0]
            query  = line[1]
            ident  = float(line[2])
            evalue = float(line[-2])

            if target in NCBI_ides and evalue<=2e-8 and ident>=75.0:
                cons.append(query[:5])

    print 'close organisms codes:', set(cons), '\n------\n'
    return list(set(cons))

# 2018 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
