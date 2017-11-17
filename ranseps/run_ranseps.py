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

def run_ranseps(genome              , cds                  , outDir='./'         , codon_table=0        , min_size=10, species_code=None,
                eval_thr=1e-8       , threads=12           ,
                eval_thr2=2e-8      , align_thr=0.0        , length_thr=0.0      , iden_thr=50.0        ,
                seps_percentage=0.25, positive_set_size=100, feature_set_size=100, negative_set_size=150,
                test_size=0.2       , folds=50            ):
    """
    RanSEPs provides a framework for bacterial genome re-annotation and novel small proteins (SEPs)
    detection adjusting the search to different genomic features that govern protein-coding capabilities.

    Only mandatory commands are:
    genome = Genome of reference in fasta format
    cds    = Coding DNA sequences of reference genome in fasta format.

    Optionally you can pass:
    outDir       = Output directory where all the results will be stored.
    codon_table  = Codon translation table. Accepted tables = [0, 4, 11]. 0 will use as START and STOP codons the ones observed for genes in --CDS
    min_size     = Minimum protein size considered.
    species_code = Identifier code used to name the putative proteins. By default RanSEPs uses the 5 first character of the genome name. Ex: mpneumoniae.fa > mpneu

    We use blast to define negative sets, wide-range criteria can be implemented
    in this step:
    eval_thr   = e-value threshold stablished to consider a hit as conserved in BlastP
    threads    = Number of threads to run BlastP
    eval_thr2  = e-value threshold stablished to consider a hit as conserved for the negative set definition
    align_thr  = Alignment threshold stablished to consider a hit as conserved for the negative set definition
    length_thr = Length threshold stablished to consider a hit as conserved for the negative set definition
    iden_thr   = Identity threshold stablished to consider a hit as conserved for the negative set definition

    Finally, to run the Random Forest Classifier you can adjust:
    seps_percentage   = Percentage of SEPs included in the positive and feature training sets
    positive_set_size = Training positive set size for each iteration
    feature_set_size  = Feature set size for each iteration
    negative_set_size = Training negative set size for each iteration
    test_size         = Percentage of sequences from the training sets used in testing
    folds             = Training positive set size for each iteration
    """

    # Create main and intermediaries folder
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    intDir = outDir+'intermediary_files/'
    if not os.path.exists(intDir):
        os.makedirs(intDir)

    # Define code
    if species_code==None:
        species_code = genome.split('/')[-1].split('.')[0][:5]

    # Genome size:
    genome  = u.load_multifasta(genome).values()[0]
    gl      = len(genome)

    # Generate DBs
    run_orfinder(genome=genome, cds=cds, outdir=intDir, ct=codon_table, min_size=min_size, species_code=species_code)
    run_gdbs(cds=cds, outdir=intDir, min_size=min_size, species_code=species_code, genome_length=gl)

    # Run Blast and find close species at the same time
    close_species = run_blaster(outdir=intDir, min_size=min_size, species_code=species_code, threshold=eval_thr, threads=threads)

    #############
    # Run RanSEPs

    # load all required info
    subcode        = intDir+species_code
    conservation   = sf.load_number_times_conserved(subcode, evalth=eval_thr2, alignth=align_thr, lenth=length_thr, identh=iden_thr, close=close_species, size=min_size)

    your_nt_seqs, your_aa_seqs, your_NCBI, your_annotation, your_lengths, your_frames, your_contra = sf.organism_info(subcode, gl, min_size)
    exclude         = [ide for ide in your_annotation if not sf.check_nt(your_nt_seqs[ide]) or not sf.check_aa(your_aa_seqs[ide])]
    your_nt_seqs    = sf.filter_dic(your_nt_seqs   , exclude)
    your_aa_seqs    = sf.filter_dic(your_aa_seqs   , exclude)
    your_annotation = sf.filter_dic(your_annotation, exclude)
    your_NCBI       = sf.filter_dic(your_NCBI      , exclude)
    print 'All required information loaded'

    # predefine sets
    big_prots, sep_prots = [], []
    for ide in your_NCBI.keys():
        if your_lengths[ide] > 100:
            big_prots.append(ide)
        else:
            sep_prots.append(ide)
    noconserved_SEPs_set = [ide for ide in your_aa_seqs.keys() if ide in conservation and conservation[ide]==0 and your_lengths[ide] <= 100 and your_frames[ide][3]=='NO']
    noconserved_BIGs_set = [ide for ide in your_aa_seqs.keys() if ide in conservation and conservation[ide]==0 and your_lengths[ide]  > 100 and your_frames[ide][3]=='NO']
    print 'Sets defined'

    sf.RanSEPs(genome=genome           , organism=species_code, nt_seqs=your_nt_seqs, aa_seqs=your_aa_seqs      , annotation=your_annotation,
               autoset=[seps_percentage, big_prots+sep_prots  , noconserved_SEPs_set, None]                     , set_sizes=[positive_set_size, feature_set_size, negative_set_size, None],
               positive_set=None       , feature_set=None     , negative_set=None   , to_exclude=[]             ,
               folds=folds             , sfolds=0             , test_size=test_size , random_state_test=None    ,
               n_estimators=100        , oob_score=1          , n_jobs=-1           , random_state=None         ,
               max_depth=None          , max_features="auto"  , min_samples_leaf=5  ,
               extension='rs'          , project_folder=intDir+'/rs_results')

    # Move and create files of interest
    print "Moving files of interest to "+outDir
    sf.handle_outputs(species_code, outDir)
    print "Thanks for using RanSEPs :D!\n||Stand and be True||"
    # All things serve the Beam.
