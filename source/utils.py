#!/usr/bin/env python

import os.path
import subprocess
import matplotlib.pyplot as plt
import collections
import scipy
import pylab
import numpy as np
from Bio import SeqIO
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable

# Several easy scripts in order to perform simple processes

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2


def polyfit2(x, y, degree):
    results = {}
    coeffs = np.polyfit(x, y, degree)
     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()
    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                      # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot
    return results

def splitn_str(your_string, n):
    """ Given a string, returns a list with that string splitted each n characters """

    return [your_string[i:i+n] for i in range(0, len(your_string), n)]


def indexes(lista, values):
    """
    Given a list and its values return a list with the indexes of the values
    in lista
    """

    indexes = []
    for element in values:
        indexes.append(lista.index(element))
    return indexes


def list_generator(filename, index=0):
    """
    Given a file, returns a set with all the values from the column[index]
    """
    results = []
    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            results.append(line[index])
    return results

def list_NA_generator(filename, index=0):
    """
    Given a file, returns a set with all the values from the column[index]
    """
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


def ins2positions(filename):
    """
    Given a ins file extract all the positions and returns them in a set
    """

    with open(filename, 'rU') as fi:
        return set([int(line.split()[0]) for line in fi.readlines()])


def set_generator(filename, index):
    """
    Given a file, returns a set with all the values from the column[index]
    """
    results = set()
    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            results.add(line[index])
    return results


def dic_generator(filename, key_index, value_index=None, header=False):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    results = {}
    with open(filename, 'r') as fi:
        for line in fi:
            if header:
                header=False
            else:
                line = line.strip().split()
                if not value_index:
                    results[int(line[key_index])] = float(line[key_index+1])
                else:
                    results[int(line[key_index])] = float(line[value_index])
    return results

def new_dic_generator(filename, key_index, value_index):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    results = {}
    with open(filename, 'r') as fi:
        results = {int(k):float(v) for k, v in [l.split()[0:2] for l in fi.readlines()]}
    return results

def str_dic_generator(filename, key_index, value_index=False, header=False, split_by=None):
    """
    Given a file, returns a dictionary where {key_index:key_index+1}
    """
    if header==True:
        header = 1

    results = {}
    with open(filename, 'r') as fi:
        for line in fi:
            if header != 0:
                header-=1
            else:
                if split_by:
                    line = line.strip().split(split_by)
                else:
                    line = line.strip().split()
                if value_index or value_index == 0:
                    results[line[key_index]] = line[value_index]
                else:
                    results[line[key_index]] = line[key_index+1]
    return results


def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def double_set_generator(filename, index):
    """
    Given a file, returns two sets with all the values from the column[index]. The first set includes
    those positions appearing in both replicas and the thirs those appearing in only one
    """

    rep_results = set()
    nonrep_results = set()

    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            try:
                if int(line[2]) == 2:
                    rep_results.add(line[0])
            except:
                nonrep_results.add(line[0])

    return (rep_results, nonrep_results)

def genes_coordinates(caps = True):
    """
    Given a three colums file returns a dictionary where the first column is the key and the other are the values
    in a list
    """

    results_dic = {}

    with open('/home/smiravet/crg/transpgrowth/datasets/essentials/gene_coordinates.txt', 'r') as fi:
        for line in fi:
            line = line.strip().split()
            if caps:
                results_dic[line[0].upper()] = [line[1], line[2]]
            else:
                results_dic[line[0]] = [line[1], line[2]]

    return results_dic

def return_ene_set(region_type):
    """
    Given E or NE return the genes with that feature
    """
    results_set = set()
    with open('/home/smiravet/crg/transpgrowth/datasets/essentials/goldsets.csv', 'r') as fi:
        for line in fi:
            line = line.strip().split()
            if line[1] == region_type:
                results_set.add(line[0].upper())
    return results_set


def process_ene_set(ENE_set, gene_coord_dic, percentage = 10):
    """
    Generate new dictionaries with the genes appearing in the NE E list and the 10% of the ORF removed (5% per side)
    """

    ENE_dic = {}

    # For the essential set
    for gene in ENE_set:
        # restore the start and end
        start = int(gene_coord_dic[gene][0])
        end = int(gene_coord_dic[gene][1])
        length = end - start
        bases_to_remove = int(round(length*(percentage/2)*0.01))
        new_start = start + bases_to_remove
        new_end = end - bases_to_remove

        ENE_dic[gene] = [new_start, new_end]

    return ENE_dic


def return_two_list(filename):
    """
    Given a file with position reads, returns two lists with those values in order to easy plot them
    """
    positions_list = []
    reads_list = []

    with open(filename, 'r') as fi:
        for line in fi:
            line = line.strip().split()
            position = int(line[0])
            reads = float(line[1])
            positions_list.append(position)
            try:
                if line[2] == '2':
                    reads_list.append(reads)
            except:
                reads_list.append(reads)

    return positions_list, reads_list


def mapping_figure(datasetA, datasetB, spanning = False):
    """
    Function to retrieve map reads per position for all the genome of pneumoniae
    """

    if spanning:
        gene_coord = gene_coordinates_dic()
        n += 1

        # Span ENE regions
        with open('/home/smiravet/crg/transpgrowth/datasets/essentials/goldsets.csv', 'r') as fi:
            for line in fi:
                line = line.strip().split()
                if line[-1] == 'E':
                    plt.axvspan(gene_coord[line[0].upper()][0], gene_coord[line[0].upper()][1], facecolor='b', alpha=0.1)
                else:
                    plt.axvspan(gene_coord[line[0].upper()][0], gene_coord[line[0].upper()][1], facecolor='g', alpha=0.1)

    # Set figure width to 24 and height to 9
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 40
    fig_size[1] = 9
    plt.rcParams["figure.figsize"] = fig_size

    # Extract the information from the nc and the rrn files
    dataA = return_two_list('/home/smiravet/crg/transpgrowth/datasets/'+datasetA)
    dataB = return_two_list('/home/smiravet/crg/transpgrowth/datasets/'+datasetB)

    plt.subplot(2, 1, 1)
    plt.title('mapping '+datasetA+' and '+datasetB)
    plt.xlim([0,816394])
    plt.bar(dataA[0], dataA[1], alpha = 0.6)
    plt.ylabel(datasetA.replace('.ins',''))

    plt.subplot(2, 1, 2)
    plt.xlim([0,816394])
    plt.bar(dataB[0], dataB[1], alpha = 0.6)
    plt.ylabel(datasetB.replace('.ins',''))

    plt.savefig('/home/smiravet/crg/transpgrowth/results/mapping/'+datasetA.replace('/','_').replace('.ins','')+datasetB.replace('/','_').replace('.ins','')+'.pdf')


def mapping_figure_from_dictionary(your_dictionary):
    # Extract the information from the nc and the rrn files

    listA= [int(k) for k,v in your_dictionary.iteritems()]
    listB= [int(v) for k,v in your_dictionary.iteritems()]

    plt.title('mapping')
    plt.xlim([0,816394])
    plt.bar(listA, listB, alpha = 0.6)
    plt.show()


def dict2file(dictionary, filename):
    """
    Writes a file where the first column is the key and the second the values
    """

    directorytosave = '/home/smiravet/crg/transpgrowth/datasets/'

    fo = open(directorytosave+filename+'.txt', 'w')

    od = collections.OrderedDict(sorted(dictionary.items()))
    for k, v in od.iteritems():
        fo.write(str(k)+'\t'+str(v)+'\n')

    fo.close()


def histogram(dataset, numBins = None, location = None):
    """ Plot a histogram for the dataset """

    if not numBins:
        numBins = len(dataset)/20

    fig = plt.figure()
    ax  = fig.add_subplot(111)

    ax.hist(dataset, numBins, color='green', alpha = 0.25)

    if location:
        plt.savefig(location)
    else:
        plt.show()


def load_multifasta(inFile):
    """ Return a dictionary wit the sequences from a multifasta file """
    your_sequences = {}
    handle = open(inFile, 'rU')
    for record in SeqIO.parse(handle, "fasta"):
        your_sequences[record.id]=str(record.seq)
    handle.close()
    return your_sequences


def load_multifasta_info(inFile):
    """ Return a dictionary wit the sequences from a multifasta file """
    your_sequences = {}
    handle = open(inFile, 'rU')
    for record in SeqIO.parse(handle, "fasta"):
        your_sequences[record.id+'//'+record.description]=str(record.seq)
    handle.close()
    return your_sequences


def load_genome(genome):
    # Determine the file type:
    if genome.endswith('gb') or genome.endswith('gbk') or genome.endswith('genbank'):
        tipo = 'genbank'
    else:
        tipo = 'fasta'
    handle = open(genome, 'rU')
    for record in SeqIO.parse(handle, tipo):
        return str(record.seq)
    handle.close()


def load_genome_DB(organism):
    """Uses load_genome function to return the sequence of the organism selected"""
    if os.path.exists('/home/smiravet/crg/dbs/smprots_DB/genomes/'+organism+'.fasta'):
        genome = load_genome('/home/smiravet/crg/dbs/smprots_DB/genomes/'+organism+'.fasta')
    else:
        genome = load_genome('/home/smiravet/crg/dbs/smprots_DB/genomes/'+organism+'.gb')

    return genome


def load_annotation(inFile):
    annotation = {}
    with open(inFile) as fi:
        for line in fi:
            line = line.strip().split()
            ide  = line[0]
            if len(line[1:]) == 2:
                st, en = sorted([int(x) for x in line[1:]])
                l      = [st, en]
            else:
                st, en = sorted([int(x) for x in line[1:] if x not in ['+', '-']])
                strand = [str(x) for x in line[1:] if x in ['+', '-']]
                l      = [st, en] + strand
            annotation[ide] = l
    return annotation


def strand_load_annotation(inFile):
    annotation = {}
    with open(inFile) as fi:
        for line in fi:
            line = line.strip().split()
            ide  = line[0]
            if len(line[1:]) == 2:
                st, en = [int(x) for x in line[1:]]
                l      = [st, en]
            else:
                st, en = [int(x) for x in line[1:] if x not in ['+', '-']]
                strand = [str(x) for x in line[1:] if x in ['+', '-']]
                l      = [st, en] + strand
            annotation[ide] = l
    return annotation


def gb2annotation(inFile):
    annotation = {}
    c = 1
    for rec in SeqIO.parse(inFile, "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
                    strand = '-' if feature.location.strand==-1 else '+'
                    start  = int(feature.location.start)
                    stop   = int(feature.location.end)
                    try:
                        genename = feature.qualifiers["locus_tag"][0]
                    except:
                        try:
                            genename = feature.qualifiers["gene"][0]
                        except:
                            genename = 'unknown'+str(c)
                            c += 1
                    annotation[genename] = [start, stop, strand]
    return annotation

def lists2dict(listA, listB):
    """ Given two lists of the same length, merge them in one dictionary """
    return dict(zip(listA, listB))


def remove_column(array, index):
    """ Remove the index column from a numpy array, index can be a list"""

    return np.delete(array, np.s_[index], axis=1)


def errorfill(x, y, yerr, color=None, linewidth=None, alpha_fill=0.3, ax=None, label=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, label=label, linewidth=linewidth)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)


# VENN 4

alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}

def get_labels(data, fill="number"):
    """
    to get a dict of labels for groups in data
    input
      data: data to get label for
      fill = ["number"|"logic"|"both"], fill with number, logic label, or both
    return
      labels: a dict of labels for different sets
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
    Out[12]:
    {'001': '001: 0',
     '010': '010: 5',
     '011': '011: 0',
     '100': '100: 3',
     '101': '101: 2',
     '110': '110: 2',
     '111': '111: 3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    if fill == "number":
        labels = {k: len(set_collections[k]) for k in set_collections}
    elif fill == "logic":
        labels = {k: k for k in set_collections}
    elif fill == "both":
        labels = {k: ("%s: %d" % (k, len(set_collections[k]))) for k in set_collections}
    else:  # invalid value
        raise Exception("invalid value for fill")
    return labels


def venn4(data=None, names=None, total=None, fill="number", show_names=True, show_plot=True, **kwds):

    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c', 'grey']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize)   # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45 , color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45 , color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))

    patches.append(Circle((200, 290), 20, color=colors[4], alpha=0.5))

    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 340); ax.set_ylim(80, 340)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, labels['1000'], fontsize=20, **alignment)
    pylab.text(280, 200, labels['0100'], fontsize=20, **alignment)
    pylab.text(155, 250, labels['0010'], fontsize=20, **alignment)
    pylab.text(245, 250, labels['0001'], fontsize=20, **alignment)
    # 2
    pylab.text(200, 115, labels['1100'], fontsize=20, **alignment)
    pylab.text(140, 225, labels['1010'], fontsize=20, **alignment)
    pylab.text(145, 155, labels['1001'], fontsize=20, **alignment)
    pylab.text(255, 155, labels['0110'], fontsize=20, **alignment)
    pylab.text(260, 225, labels['0101'], fontsize=20, **alignment)
    pylab.text(200, 240, labels['0011'], fontsize=20, **alignment)
    # 3
    pylab.text(235, 205, labels['0111'], fontsize=20, **alignment)
    pylab.text(165, 205, labels['1011'], fontsize=20, **alignment)
    pylab.text(225, 135, labels['1110'], fontsize=20, **alignment)
    pylab.text(175, 135, labels['1101'], fontsize=20, **alignment)
    # 4
    pylab.text(200, 175, labels['1111'], fontsize=20, **alignment)

    # Compute no classified
    pylab.text(200, 288, str(len(total.difference(data[0], data[1], data[2], data[3]))), fontsize=20, **alignment)
    pylab.text(200, 315, 'Undetected', fontsize=20, **alignment)

    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=20, **alignment)
        pylab.text(290, 110, names[1], fontsize=20, **alignment)
        pylab.text(130, 275, names[2], fontsize=20, **alignment)
        pylab.text(270, 275, names[3], fontsize=20, **alignment)

    # leg = ax.legend([names[0], names[2], names[3], names[1]] , loc='best', fancybox=True)
    # leg.get_frame().set_alpha(0.5)

    if show_plot:
        pylab.show()



##### FOR SEQUENCES

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    return ''.join([complement[k] if k in complement else 'N' for k in seq][::-1])
