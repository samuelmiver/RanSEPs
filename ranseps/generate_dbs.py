#!/usr/bin/env python

# 2017 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

import glob
import sys, os
import os.path
import utils as u
from Bio.Seq import Seq


def extract_functions(anno_dict):

    dic = {}
    for info, seq in anno_dict.iteritems():
        try:
            info       = info.split('[')[1:]
            gene, prot = ''.join(info).split(']')[0:2]
            gene       = gene.replace('gene=','')
            protein    = prot.replace('protein=','')
            dic[gene+' |'+protein] = seq
        except:
            gene = info.split(' ')[0].replace('>', '')
            dic[gene+' |'+gene] = seq
    return dic


def pair_id(pred, anno, outfile):

    # Depurate the annotated dictionary to have
    # {'gene | protein':seq}
    anno = extract_functions(anno)

    # Find matches
    merged = pred
    merged.update(anno)
    flipped = {}
    lengths = {}
    for key, value in merged.items():
        lengths[key] = (len(value)-3)/3
        if value[-45:] not in flipped:
            flipped[value[-45:]] = [key]
        else:
            flipped[value[-45:]].append(key)

    # We need to do only a portion of the sequence as the start is variable depending on
    # annotation but the stop is conserved! We set the size in 57 which is the size for 19 aa
    results = []
    c = 0
    for k, v in flipped.iteritems():
        if len(v) > 1:
            c += 1
            if '|' in v[1] and not '|' in v[0]:
                g, p = v[1].split(' | ')
                results.append([v[0], g, p , str(lengths[v[0]])])
            elif '|' in v[0] and not '|' in v[1]:
                g, p = v[0].split(' | ')
                results.append([v[1], g, p, str(lengths[v[1]])])
            else:
                pass

    print len(anno), len(results)
    with open(outfile, 'w') as fo:
        for row in sorted(results):
            fo.write('\t'.join(row)+'\n')


def annotation_file(g_size, predd, dir_handle, min_size):

    genes = []
    for info, seq in predd.iteritems():
        info = info.replace(' ', '').split('//')[1].split('|')
        sten = [int(x) for x in info[2].split('..')]
        gene, strand, start, end  = info[0], info[1], sten[0], sten[1]
        if strand == '+':
            end += 3
        else:
            end -= 3
            if end < 0:
                end += g_size-1
        genes.append([gene, start, end, strand])

    # Handle for DB table
    fo2 = open(dir_handle+'_info.txt', 'w')
    fo2.write('genename\tassigned\tstrand\tstart\tstop\tlength\taa_seq\n')

    # Dictionary of paired ides and aa seqs
    pair_id = u.str_dic_generator(dir_handle+'_pairs.txt', 0, 1)
    aa_seqs = u.load_multifasta(dir_handle+'_small_aa.fa')

    # Generate the files
    with open(dir_handle+'_annotation.txt', 'w') as fo:
        for row in sorted(genes, key = lambda x: x[1]):
            fo.write('\t'.join([str(i) for i in row])+'\n')
            genename = row[0]
            l = len(aa_seqs[genename])
            if l >= min_size:
                # Table
                if genename in pair_id:
                    pair = pair_id[genename]
                else:
                    pair = 'putative'
                new_row = [genename, pair, row[3], row[1], row[2], l, str(aa_seqs[genename])]
                fo2.write('\t'.join([str(i) for i in new_row])+'\n')
    fo2.close()


def run_gdbs(cds, outdir, min_size, species_code, genome_length):

    dir_handle = outdir+species_code
    annotated = u.load_multifasta_info(cds)
    predicted = u.load_multifasta(dir_handle+'_small_nt.fa')
    predicted_info = u.load_multifasta_info(dir_handle+'_small_nt.fa')

    # Annotation
    print 'Generating aa fasta'
    fo1 = open(dir_handle+'_small_aa.fa', 'w')
    fo2 = open(dir_handle+'_aa_lengths.txt', 'w')
    for ide, seq in predicted_info.iteritems():
        seq = Seq(seq)
        translated = 'M'+str(seq.translate())[1:-1].replace('*', 'W')
        l = len(translated)
        if l >= min_size:
            fo1.write('>'+ide.split('//')[1]+'\n'+translated+'\n')
            fo2.write(ide.split('//')[0]+'\t'+str(l)+'\n')

    fo1.close()
    fo2.close()

    # Pair ides
    print 'Pairing annotation\n----\n'
    pair_id(predicted, annotated, dir_handle+'_pairs.txt')

    # Write annotation
    print 'generating annotation\n----\n'
    annotation_file(genome_length, predicted_info, dir_handle, min_size)

# 2017 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
