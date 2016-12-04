#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os

'''
Create a head folder that is going to contain all the noise reduced folders
'''
def foldercreator(dir):
    directory = dir + str('_noicereduced')
    if not os.path.exists(directory):
        os.makedirs(directory)

'''
Function for creating subdirectories (The noise reduced folders)
'''
def createdirectory(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

'''
A function that reads the inputfile and loads the sequences into a list
'''
sequences = []
name=[]
#Program som läser filen:
def readfile(filename):
    alignment = ''
    try:
        sequences = []
        alignment = AlignIO.read(filename, 'fasta')
        infile = open(filename, 'r')
        #for r in SeqIO.parse(infile, 'fasta'):
        for r in alignment:
            if len(r.seq) == 0:
                raise Exception
            empty = SeqRecord(Seq(""), r.id, r.name, r.description) #Creates an empty SeqRecord with the correct id, name and description
            sequences.append(empty) #Adds the SeqRecord to the sequence list
        infile.close()
    except:
        print 'Error! File did not open correctly. Please check your file.'
        sys.exit(0)  #closes the program
        
    return sequences, alignment


'''
A function which counts the occurences of each amino acid character and "-" in the given column, updates the value in the dictionary for them.
Returns the updated dictionary.
'''
def count_character(column):
    amino_acids = {'A':0,'R':0,'K':0,'M':0,'F':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0, '-':0}
    for char in column:
        if char in amino_acids:
            amino_acids[char] = amino_acids[char] + 1 #för varje gång vi hittar bokstaven i stringen i dictionaryt så räkna upp.
        else:    
            return amino_acids, True
        
    return amino_acids, False

'''
A function that calls for count_character to check the occurences of amino acids or indels and returns
True if the given conditions are fulfilled, otherwise False.
'''
def checknoise(column):
    amino_acids, outcome =count_character(column)
    if outcome == True:
        print 'Error: A column contained a character which was not an amino acid. Column was therefore removed.'
        return True
    else:
        lessthantwo = 0
        unique = 0
        indels = 0
        total_aa=0
        for acids in amino_acids:
            
            if acids != '-': #We don't want to check the '-' in the dictionary (blir fel då)
        #Check how many amino acids we have in our column:
                if amino_acids[acids] != 0:
                    total_aa += amino_acids[acids]
                
        #Check how many acids that are unique
                if amino_acids[acids] == 1:
                    unique += 1
        #Check how many acids appears less than twice
                if amino_acids[acids] <= 2:
                        lessthantwo += 1
        #More than 50% indels:        
            else:
                indels = amino_acids['-']

            #Conditions:

            #More than 50% indels:   tidigare hade vi >=0.5        
        if float(indels/len(column)) > 0.5:
            return True
        # Least 50% of amino acids are unique:
        if unique/total_aa >= 0.5:
            return True

        # No amino acid appears more than twice:
        if lessthantwo == 20:
            return True

    return False


'''
Function which takes the infile and sends it in to the other functions and gives the final output as a fine file saved on your wonderful computer.
'''
def noise_reduction(path, infile, filename):
    seqfile,alignment = readfile(infile)

    if seqfile ==[]:
        print 'Error: There must be a problem with your file, it is probably empty, please check it'
        #sys.exit('Program will be closed')

    else:
        allNoisy= True #Nothing as been inserted and all has been removed
        SIZE = alignment.get_alignment_length()

        for c in range(SIZE):
            column = alignment[:, c]
            if not checknoise(column): #if not true
                for k in range(len(seqfile)):
                    seqfile[k] += column[k]
                allNoisy=False #something has been inserted


        if allNoisy:
            print 'Error: All the columns have been removed, too much noise! Please take another file'    
            sys.exit()
        else:         
            outfile = path + '/resultat_' + str(filename)
            output = MultipleSeqAlignment(seqfile)
            AlignIO.write(output, outfile, 'fasta')

