

import sys
import os
import argparse

from superimposer import Superimposer

def cmd():

    """Command line program processor."""

    superimposer=Superimposer()
    first_pdb = options.infile1
    second_pdb = options.infile2
    output_pdb = options.outfile

    if not first_pdb and not second_pdb:
        sys.stderr.write("\nThere is no input files selected. Please, select two PDB files\n\n")
        sys.exit()

    elif not output_pdb:
        sys.stderr.write("\nThere is no PDB output file selected.\n\n")
        sys.exit()

    elif not first_pdb.endswith(".pdb") or not second_pdb.endswith(".pdb"):
        sys.stderr.write("\nInput given not pdb files. Please, ensure the inputs are PDBs\n\n")
        sys.exit()

    elif not output_pdb.endswith(".pdb"):
        sys.stderr.write("\nThe output file has to be include the .pdb extension\n\n")
        sys.exit()

    else:
        print(">> Files accepted, starting Overlapy...")

        #Run the superimposition
        superimposer.parse(str(first_pdb), str(second_pdb))

        #Selecting chains (if they are set)
        if options.chains1 or options.chains2:
             print(">> Selecting chains {}, {}".format(options.chains1, options.chains2))
             chains1 = list(options.chains1) if options.chains1 else []
             chains2 = list(options.chains2) if options.chains2 else []
             superimposer.select_chains(chains1, chains2)

        #Superimpose:
        superimposer.superimpose()

        print(">> Overlapy has finished correctly. Processing outputs...")

        #Take the rmsd, matrix and alignment
        rmsd = superimposer.rmsd
        matrix = superimposer.rotation_matrix
        alignment = superimposer.get_multiple_sequence_alignment()

    #Output generation

    superimposer.save_superimposed_pdb(str(output_pdb))
    print(">> PDB with the structures superimposed generated.")

    #The matrix
    if options.matrix:
        numpy.savetxt(("matrix.tsv"), matrix, delimiter="\t")
        print(">> Rotation matrix generated.")

    else:
        print(">> Rotation matrix generated:\n\n{}\n".format(matrix))

    #RMSD
    if options.rmsd:
        with open("rmsd.txt", 'w') as o:
            o.write("PDB1:{}\nPDB2:{}\n\nRMSD:\n{}\n".format(first_pdb, second_pdb, rmsd))
        print(">> RMSD text file generated.")

    else:
        print(">> RMSD:\t{}\n".format(rmsd))

    #Alignment
    if options.align:
        with open("alignment.clustal", 'w') as o:
            o.write(alignment.format('clustal'))
        print(">> alignment.clustal generated.")

    else:
        print(">> Alignment:\n\n{}".format(alignment.format('clustal')))


#PARSER
parser = argparse.ArgumentParser(description = "Overlapy performs the superimposition of two protein structures in PDB")

parser.add_argument('-i1', '--input1',
                    dest="infile1",
                    action="store",
                    default=None,
                    help="First PDB input")

parser.add_argument('-i2', '--input2',
                    dest="infile2",
                    action="store",
                    help="Second PDB input")

parser.add_argument('-c1', '--chains1',
                    dest="chains1",
                    action="store",
                    default=None,
                    help="Chains selected from the first PDB, one letter code (H, L, A, B...). eg. -c1 HL")

parser.add_argument('-c2', '--chains2',
                    dest="chains2",
                    action="store",
                    default=None,
                    help="Chains selected from the second PDB, one letter code (H, L, A, B...). eg. -c2 HL")

parser.add_argument('-o', '--output',
                    dest="outfile",
                    action="store",
                    default=None,
                    help="PDB outputfile")

parser.add_argument('-m', '--matrix',
                    dest="matrix",
                    action="store_true",
                    default=False,
                    help="If -m, matrix stored in .tsv format")

parser.add_argument('-r', '--RMSD',
                    dest="rmsd",
                    action="store_true",
                    default=False,
                    help="If -r, RMSD value exported in .txt file")

parser.add_argument('-a', '--alignment',
                    dest="align",
                    action="store_true",
                    default=False,
                    help="If -a, alignment stored in .clustal format")

options = parser.parse_args()

if __name__=="__main__":
    cmd()

