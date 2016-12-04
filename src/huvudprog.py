#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import os
import tempfile
from reducenoise import * #Imports all functions from reducenoise.py
import dendropy
import csv
from dendropy import Tree
from dendropy.calculate import treecompare

infiles = sys.argv[1] #Input should be folder containing the fasta files
tns = dendropy.TaxonNamespace()


'''
                Head program for running the reduce noise program
'''


'''
Run FastPhylo fastprot

fastprot -> estimates the evolutionary distance between aligned protein sequences.
It implements two methods for calculating the distance between protein sequences: the maximum likelihood (ML) and the
expected distance.
'''
def fastprot(path):
    tmp = tempfile.mkstemp()
    os.system("fastprot -m -o " + tmp[1] + " " + path)
    os.close(tmp[0])

    return tmp[1]
	

'''
Run FastPhylo fnj

fnj -> The fnj program implements three tree reconstruction methods, and the default is FNJ
'''
def fnj(path):
    tmp = tempfile.mkstemp()
    os.system("fnj -O newick -m FNJ -o " + tmp[1] + " " + path)

    return tmp[0]


'''
Compare the given tree with the reference tree 
in regards of symmetric difference
'''
def compare_tree(recovered, reference_tree):
    # Open temporary stream and parse tree
    
    stream = os.fdopen(recovered)
    tree = Tree.get_from_stream(stream, schema='newick', taxon_namespace=tns)
	
    compared = treecompare.symmetric_difference(reference_tree, tree)

    stream.close()

    return compared

    

def find_reference_tree(directory):
    files = os.listdir(directory)

    for filename in files:
        if filename[0] != '.':
            if filename.endswith('.tree'):
                #If it is an .tree file then it is the reference tree
                return Tree.get_from_path(directory + "/" + filename, schema='newick', taxon_namespace=tns)

'''
Take the input files
'''

#Create a folder named "appbio11_noicereduced" to store the noise reduced alignments

input_file = sys.argv[1]
lfolder_name=len(input_file) #Takes out the lenght of the name of the folder

for directory in os.walk(input_file):
    reduced_directory = directory[0] #Takes the name of the folder
    reduced_directory=reduced_directory[lfolder_name+1:] #Takes all subdirectories
    if (reduced_directory != "") :
        path=input_file + '_noicereduced/' + reduced_directory
        createdirectory(path)
        createdirectory("../result")
        reference_tree = find_reference_tree(directory[0])

        files = os.listdir(directory[0])

        with open("../result/" + reduced_directory + ".csv", 'wb') as csvfile:

            for filename in files:
                if filename[0] != '.':
                   if filename.endswith('.msl'): #If it is an .msl file it should be noise reduced
                        noise_reduction(path, directory[0] + "/" +filename, filename)
                        '''
                        Make a tree of the alignment files
                        '''
                        alignment_fastprot_tree = fastprot(directory[0] + "/" +filename)
                        alignment_fnj_tree = fnj(alignment_fastprot_tree)
                        alignment_compare_tree = compare_tree(alignment_fnj_tree, reference_tree)

                        '''
                        Make a tree of the noise reduced files
                        '''
                        noise_reduced_fastprot_tree = fastprot(path + "/resultat_" + filename)
                        noise_reduced_fnj_tree = fnj(noise_reduced_fastprot_tree)
                        noise_reduced_compare_tree = compare_tree(noise_reduced_fnj_tree, reference_tree)
    
                        '''
                        Calculate the difference between the noise reduced tree and the alignment tree
                        '''
                        difference = noise_reduced_compare_tree - alignment_compare_tree
           
                    
                        writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        writer.writerow([filename, alignment_compare_tree, noise_reduced_compare_tree, difference])
                
                    #sys.stdout.write('%s\t%d\t%d\t%d\n' % (filename, alignment_compare_tree, noise_reduced_compare_tree, difference))

    

'''
DendroPy --> DendroPy is a Python library for phylogenetic computing.
It provides classes and functions for the simulation, processing, and manipulation of phylogenetic trees and character matrices
'''

