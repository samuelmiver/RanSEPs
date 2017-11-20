#!/usr/bin/env python

#############################################################
#
# Author : Miravet-Verde, Samuel
# Written : 06/02/2016
# Last updated : 11/17/2017
# [2017] - Centre de Regulació Genòmica (CRG) - All Rights Reserved
#############################################################

import sys, os
import re
import utils as u
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein

def run_orfinder(genome, cds, outdir, ct, min_size, species_code):
    """
    Given a genome in fasta or genbank format
    Generates the nt and aa databases with proteins larger than min size

    mode 1 = requires fasta file with genes including real START and STOP codons,
    mode 4 = use codon table 4 (mycoplasma and spiroplasma)
    mode 11 =  use codon table 11
    """

    genome = Seq(genome)
    print 'running the analysis on:'+species_code+'\n------\n'

    # Define the start and stop codons based on what we see in the annotation
    if ct==0:
        annotated_seqs = u.load_multifasta(cds)
        start_codons = set()
        stop_codons  = set()
        for fasta_identifier, sequence_of_interest in annotated_seqs.iteritems():
            start_codons.add(sequence_of_interest[:3])
            stop_codons.add(sequence_of_interest[-3:])
        start_codons = list(start_codons)
        stop_codons  = list(stop_codons)
    elif ct==4:
        start_codons = ['ATG', 'TTG', 'GTG']
        stop_codons  = ['TAG', 'TAA']
    elif ct==11:
        start_codons = ['ATG', 'TTG', 'GTG']
        stop_codons  = ['TAG', 'TAA', 'TGA']
    else:
        print '\n\n--> '+str(ct)+' is not an accepted translation table\n'
        sys.exit(1)

    print 'start codons:\n'+'\t'.join(start_codons)+'\n------\n'
    print 'stop  codons:\n'+'\t'.join(stop_codons)+'\n------\n'

    # Transform them into a regexp
    regexpstarts = '(?='+'|'.join(start_codons)+')'
    regexpstops  = '(?='+'|'.join(stop_codons)+')'

    # Some internal function
    def return_frame(size, selected_frame):
        return set(range(selected_frame-1, size, 3))

    def closest(start, ends):
        return min(filter(lambda x: start < x, ends))

    def interval(start, end, midvalues):
        return set([i for i in midvalues if start < i < end])

    def return_ORFs(frame, starts, ends, results, genome, strand):
        pass_set = set([])
        genome_len  = len(genome)

        frame_set = return_frame(genome_len*2, frame)
        start_cds = sorted([x+1 for x in starts.intersection(frame_set)])
        end_cds   = sorted([x+1 for x in ends.intersection(frame_set)])

        if len(start_cds) >= 1 and len(end_cds) >= 1:
            for st in start_cds:
                if st not in pass_set:
                    if st < max(end_cds):
                        en = closest(st, end_cds)

                        # Exclude the next start to take the largest ORF
                        try:
                            next_starts = interval(st, en, start_cds)
                            pass_set.update(next_starts)
                        except:
                            pass

                        if en >= genome_len:
                            seq = (genome+genome)[st-1:en+2]
                        else:
                            seq = genome[st-1:en+2]

                        # Process seq
                        results.append([st, en-1, str(seq), strand])
        return results


    # Define the number of stops and starts
    results = []
    total_stops  = 0.0
    total_starts = 0.0

    # Define the genomes
    list_of_genomes = [('+', genome), ('-', genome.reverse_complement())]

    # Run the prediction
    for strand, sequence in list_of_genomes:
        # Define the starts and ends
        start_indexes = set([x for x in [m.start() for m in re.finditer(regexpstarts, str(sequence))]])
        end_indexes = set([m.start() for m in re.finditer(regexpstops, str(sequence+sequence))])

        print 'strand: '+strand+'\n-------\n1st frame...'
        results = return_ORFs(1, start_indexes, end_indexes, results, sequence, strand)
        print '2nd frame...'
        results = return_ORFs(2, start_indexes, end_indexes, results, sequence, strand)
        print '3rd frame...'
        results = return_ORFs(3, start_indexes, end_indexes, results, sequence, strand)

    print 'number of SEPs found:'+str(len(results))+'\n------\n'

    # Write in fasta files
    font = open(outdir+species_code+'_small_nt.fa', 'w')

    maxi = len(str(len(results)+1))
    identifiers = [species_code+str(i+1).zfill(maxi) for i in range(0, len(results))]

    records_nt = []
    for i in range(0, len(results)):
        st, en, nt_seq, strand = results[i]
        if len(nt_seq)-3 >= min_size*3.0:
            if strand == '-':
                st = len(genome)+1-st
                en = len(genome)+1-en
            records_nt.append(SeqRecord(Seq(nt_seq, generic_dna), identifiers[i], description='|'+strand+'|'+str(st)+'..'+str(en)))

    print 'number of SEPs after filter:'+str(len(records_nt))+'\n------\n'
    SeqIO.write(records_nt, font, 'fasta')

# [2017] - Centre de Regulació Genòmica (CRG) - All Rights Reserved
