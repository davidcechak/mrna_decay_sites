# Author: Anish Raju, NIH SIP (LGG)

"""This file contains various functions for processing data to produce inputs/outputs for models.
"""

from Bio import SeqIO
import argparse
import numpy as np
import pandas as pd
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import train_test_split

def parse_data(file, x, y, rbps, rbp_dict, type, rep, onehot, length=200):
    """Iterates through positive fasta file, pads sequences if necessary, and appends to corresponding arrays.
    @param file     Fasta file to read. Must be positive for decay.
    @param x        Array to append 200nt sequences to.
    @param y        Array to append pos/neg labels to. neg=0, pos=1
    @param rbps     Array to append RBP information to
    @param rbp_dict Dictionary matching a transcript to the start/ends of binding sites for an RBP
    @param type     Whether a file contains reads from 5' head (5), 3' tail (3), or the CDS
    @param rep      Whether to distinguish between soft-masked and hard-masked nucleotides.
    @param onehot   Whether to encode "x" to one-hot or leave nucleotides
    """
    fasta = SeqIO.parse(open(file), format="fasta")
    for fastaseq in fasta:
        sequence = fastaseq.seq

        # NOT USING NOW # Comparing read coordinates to match RBP binding locations
        # prelength = len(sequence)
        # rbp_array = compare_rbp(rbp_dict, fastaseq.id, prelength)
        # # If sequence is <200nt in length, pad RBP array to match below padding
        # if len(rbp_array)<length:
        #     if type=="5" or type==5:
        #         rbp_array = np.concatenate((np.zeros(length-prelength), rbp_array))
        #     elif type=="3" or type==3:
        #         rbp_array = np.concatenate((rbp_array, np.zeros(length-prelength)))
        # rbps.append(rbp_array)
           
        # Sequence padding for reads<200nt in length
        sequence = pad(sequence, length, type)
        
        # Shorten to desired len if longer
        if len(sequence) > length:
            longer_by = len(sequence) - length
            # sequence = sequence[longer_by//2:]
            # sequence = sequence[:(-1)*longer_by//2]
            sequence = sequence[longer_by//2-1:(len(sequence)-longer_by//2)]
 
        # x.append(encode(list(sequence), rep)) 
        if onehot:
            x.append(encode(list(sequence), rep))
        else:
            x.append(list(sequence))
        y.append(1)

def parse_neg(file, x, y, rbps, rbp_dict, rep, onehot, length=200):
    """Iterates through reference fasta file, selects random 200nt sequences, pads and appends to corresponding arrays.
    @param file     Fasta file to read. Must be positive for decay.
    @param x        Array to append 200nt sequences to.
    @param y        Array to append pos/neg labels to. neg=0, pos=1
    @param rbps     Array to append RBP information to
    @param rbp_dict Dictionary matching a transcript to the start/ends of binding sites for an RBP
    @param rep      Whether to distinguish between soft-masked and hard-masked nucleotides.
    @param onehot   Whether to encode "x" to one-hot or leave nucleotides
    """
    fasta = SeqIO.parse(open(file), format="fasta")
    for fastaseq in fasta:
        sequence = fastaseq.seq

        # Compare entire read to RBP binding locations
        # prelength = len(sequence)
        # rbp_array = compare_neg_rbp(rbp_dict, fastaseq.id, prelength)

        # Randomly select a 200nt section, with padding if the section is from the ends
        sequence = "N"*length + sequence + "N"*length
        center = np.random.randint(length, len(sequence)-length)
        sequence = sequence[center-length//2 : center+length//2]
        # CHECK
        if len(sequence) > length:
            print(len(sequence))
        # x.append(encode(list(sequence), rep))
        if onehot:
            x.append(encode(list(sequence), rep))
        else:
            x.append(list(sequence))
        y.append(0)

        # Extract the portion of the RBP array that corresponds to the randomly selected 200nt region
        # rbp_array = np.concatenate((np.zeros(length), rbp_array, np.zeros(length)))
        # rbp_array = rbp_array[center-length//2 : center+length//2]
        # rbps.append(rbp_array)

def compare_rbp(rbp_dict, fid, prelength):
    """Helper function for parse_data to compare reads to match RBP binding locations. Returns an array of same length as sequence, where 1 indicates
    that that position lies within the RBP binding site and 0 otherwise.
    @param rbp_dict     Dictionary matching a transcript to the starts/ends of binding sites for an RBP
    @param fid          ID from FASTA file, in format ENSTXX...XX:start-end
    @param prelength    Length of read
    """
    name, start, end = split_id(fid)
    if name in rbp_dict:
        rbp_start, rbp_end = rbp_dict[name][0], rbp_dict[name][1]
        rbp_array = np.zeros(prelength)
        if start<=rbp_start and rbp_start<=end:
            for i in range(rbp_start-start, rbp_end-start+1):
                if i<prelength:
                    rbp_array[i] = 1
    else:
        rbp_array = np.zeros(prelength)
    return rbp_array

def compare_neg_rbp(rbp_dict, fid, prelength):
    """Helper function for parse_neg to compare reads to match RBP binding locations. Returns an array of same length as sequence, where 1 indicates
    that that position lies within the RBP binding site and 0 otherwise.
    @param rbp_dict     Dictionary matching a transcript to the starts/ends of binding sites for an RBP
    @param fid          ID from FASTA file, in format ENSTXX...XX:start-end
    @param prelength    Length of read
    """
    name = fid
    if name in rbp_dict:
        rbp_start, rbp_end = rbp_dict[name][0], rbp_dict[name][1]
        rbp_array = np.zeros(prelength)
        if rbp_start<=prelength:
            for i in range(rbp_start, rbp_end+1):
                if i<prelength:
                    rbp_array[i] = 1
    else:
        rbp_array = np.zeros(prelength)
    return rbp_array

def split_id(id):
    """Helper function to split string id in format "ENSTXX...XX:start-end". Returns "ENSTXX...XX", start, end.
    If id only contains the transcript id, then simply returns it in an array.
    """
    if ":" not in id:
        return [id]
    s = id.split(":")
    name = s[0]
    coords = s[1].split("-")
    start, end = int(coords[0]), int(coords[1])
    return name, start, end

def encode(sequence, rep=True): # TODO make rep arugment
    """Performs one-hot encoding of nucleotide sequence.
    @param sequence Original nucleotide sequence to encode.
    @param rep      Whether to distinguish between soft-masked and hard-masked nucleotides
    """
    to_integer = {}
    total_nuc = ["A", "T", "C", "G", "N", "a", "t", "c", "g"]
    if rep==False:
        total_nuc = ["A", "T", "C", "G", "N"]
        sequence = [x.upper() for x in sequence]
        
    # Map each nucleotide to integer
    for i in range(len(total_nuc)):
        to_integer[total_nuc[i]] = i

    # Convert input sequence to corresponding integers
    for i in range(len(sequence)):
        sequence[i] = to_integer[sequence[i]]
    
    # One hot encode with keras backend
    onehot_encoding = to_categorical(sequence, num_classes=len(total_nuc))
    return onehot_encoding

def pad(sequence, cutlength, type):
    """Pads reads from 5' or 3' ends which are <200nt in length by calling helper functions. If not necessary, then returns without changing.
    @param sequence  Nucleotide sequence to pad
    @param cutlength Desired length of final sequence
    @param type      5 if from 5' head, 3 if from 3' tail
    """
    if type=="5" or type==5:
        return pad_5p(sequence, cutlength)
    elif type=="3" or type==3:
        return pad_3p(sequence, cutlength)
    else:
        return sequence

def pad_5p(sequence, cutlength):
    """Concatenates sequence with a string of "N"'s at beginning to bring total length to cutlength
    @param cutlength Desired length of final sequence
    """
    padded = "N"*(cutlength-len(sequence)) + sequence
    return padded

def pad_3p(sequence, cutlength):
    """Concatenates sequence with a string of "N"'s at end to bring total length to cutlength
    @param cutlength Desired length of final sequence
    """
    padded = sequence + "N"*(cutlength-len(sequence))
    return padded

def get_rbp_dict(rbp_bedfile):
    """Returns dictionary matching each transcript to the starts/ends of binding sites for an RBP.
    """
    rbp_list = pd.read_csv(rbp_bedfile, sep="\t", header=None, names=["transcript_id", "start", "end", "database", "length", "strand"])
    rbp_dict = {}
    for index, row in rbp_list.iterrows():
        name, start, end = row['transcript_id'], int(row['start']), int(row['end'])
        rbp_dict[name]=[start, end]
    return rbp_dict

def reformat(x_train, x_test, rbps_train, rbps_test, y_train, y_test, n_seqlength, n_features, onehot):
    """ Converts inputs and outputs into format model can use. Converts all to numpy arrays.
    """
    lentrain, lentest = len(x_train), len(x_test)
    lenrbptrain, lenrbptest = len(rbps_train), len(rbps_test)
    # Reshape input sequences
    if onehot:
        x_train = np.concatenate(x_train).reshape((lentrain, n_seqlength, n_features))
        x_test = np.concatenate(x_test).reshape((lentest, n_seqlength, n_features))
    else:
        x_train = np.array(x_train)
        x_test = np.array(x_test)

    rbps_train = np.concatenate(rbps_train).reshape((lenrbptrain, n_seqlength, 1))
    rbps_test = np.concatenate(rbps_test).reshape((lenrbptest, n_seqlength, 1))
    # Convert targets to numpy array
    y_train = np.array(y_train)
    y_test = np.array(y_test)
    return x_train, x_test, rbps_train, rbps_test, y_train, y_test

def motif_test():
    """Generates random sequences to test initial CNN, with artificial motifs inserted into half of them.
    """
    x, y = [], []
    for i in range(1000):
        # positive
        seq = [random.choice(['A', 'G', 'T', 'C']) for j in range(200)]
        seq = ["A"] * 9 + seq[9:100] + list("GATCCATG") + seq[108:192] + list("GATCCATG")
        x.append(encode(seq))
        y.append(1)
    for i in range(1000):
        # negative
        seq = [random.choice(['A', 'G', 'T', 'C']) for j in range(200)]
        x.append(encode(seq))
        y.append(0)

    return train_test_split(x, y)

