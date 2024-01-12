import re
import random

from enum import Enum, auto

class Variant_consequence(Enum):
    Stop_gained = 'Stop_gained'
    Stop_lost = 'Stop_lost'
    Start_lost = 'Start_lost'
    Missense_variant = 'Missense_variant'
    Amino_retained = 'Amino_retained'
    
class Clinical_significance(Enum):
    Likely_benign = 'Likely_benign'
    Benign = 'Benign'
    Conflicting_interpretations_of_pathogenicity = 'Conflicting_interpretations_of_pathogenicity'
    Uncertain_significance = 'Uncertain_significance'
    Benign_Likely_benign = 'Benign/Likely_benign'

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F',  # Phenylalanine
    'TTA': 'L', 'TTG': 'L',  # Leucine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',  # Serine
    'TAT': 'Y', 'TAC': 'Y',  # Tyrosine
    'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop codon
    'TGT': 'C', 'TGC': 'C',  # Cysteine
    'TGG': 'W',  # Tryptophan
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',  # Leucine
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline
    'CAT': 'H', 'CAC': 'H',  # Histidine
    'CAA': 'Q', 'CAG': 'Q',  # Glutamine
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',  # Arginine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',  # Isoleucine
    'ATG': 'M',  # Methionine
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine
    'AAT': 'N', 'AAC': 'N',  # Asparagine
    'AAA': 'K', 'AAG': 'K',  # Lysine
    'AGT': 'S', 'AGC': 'S',  # Serine
    'AGA': 'R', 'AGG': 'R',  # Arginine
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Valine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine
    'GAT': 'D', 'GAC': 'D',  # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'   # Glycine
}


def compare_sequence(seq_a, seq_b):
    len1 = len(seq_a)
    len2 = len(seq_b)
    print("Length of Sequence A: " + str(len1))    
    print("Length of Sequence B: " + str(len2))
    mismatches = []
    mis_indices = []
    for i, pos in enumerate(range(0, min(len1, len2))) :
        if seq_a[pos] == seq_b[pos]:
            mismatches.append('║')
            mis_indices.append(i)
        else:
            # mismatches.append('ˑ')
            mismatches.append('-')
    # print (mis_indices)
    print (seq_a)
    print ("".join(mismatches))
    print (seq_b, '\n')
    
    
def generate_random_dna_sequence(length):
    if length <= 0:
        return "Invalid input: Length must be a positive integer"
    return ''.join(random.choice('ACGT') for _ in range(length))


def find_codon_frame_index(distance):
    if distance < 0:
        return "Invalid input: Distance must be a non-negative integer"
    
    frame_index = distance % 3
    # if frame_index == 0:
    #     return 3
    # else:
    return frame_index


# Example usage
# print(translate_codon_dna("ATG"))  # Methionine (M)
# print(translate_codon_dna("TAA"))  # Stop codon (*)
def translate_codon_to_amino(codon):
    return GENETIC_CODE.get(codon.upper(), "Invalid codon")


def get_amino_change(codon_old, codon_new):
    amino_old = translate_codon_to_amino(codon_old)
    amino_new = translate_codon_to_amino(codon_new)
    if(amino_old == amino_new):
        change = amino_old
    else:
        change = amino_old + '/' + amino_new
    return change


# Example usage
# print(classify_mutation("ATG", "ATA"))  # Start_lost
# print(classify_mutation("TAA", "GAA"))  # Stop_lost
# print(classify_mutation("TTG", "TAG"))  # Stop_gained
# print(classify_mutation("GAG", "GAA"))  # Missense_variant
# print(classify_mutation("GGC", "GGT").value)  # Amino_retained
def classify_variant_consequence(original_codon, mutated_codon):
    original_aa = GENETIC_CODE.get(original_codon.upper(), None)
    mutated_aa = GENETIC_CODE.get(mutated_codon.upper(), None)

    if original_aa is None or mutated_aa is None:
        raise ValueError("Codon is None")
        
    if original_aa == 'M' and mutated_aa != 'M':
        return Variant_consequence.Start_lost
    elif original_aa != '*' and mutated_aa == '*':
        return Variant_consequence.Stop_gained
    elif original_aa == '*' and mutated_aa != '*':
        return Variant_consequence.Stop_lost
    elif original_aa == mutated_aa:
        return Variant_consequence.Amino_retained
    else:
        return Variant_consequence.Missense_variant
    
    
def convert_to_eqtl_format(hgvs_string):
    """
    Convert an HGVS string to eQTL variant format.

    Args:
    hgvs_string (str): An HGVS string, e.g., "NC_000013.11:g.114184638T>C".

    Returns:
    str: Converted string in https://www.ebi.ac.uk/eqtl/Data_access/ and https://elixir.ut.ee/eqtl format, e.g., "chr13_114184638_T_C".
    """
#     # Extract chromosome number, position, and alleles
#     parts = hgvs_string.split(':')
#     chrom_number = parts[0].split('.')[0][-2:]  # Gets the last two characters of the chromosome reference
#     position_alleles = parts[1].split('g.')[1]
#     position, alleles = position_alleles.split('>')
#     ref_allele, alt_allele = alleles[0], alleles[-1]

    # Extract chromosome number, position, and alleles
    parts = hgvs_string.split(':')
    chrom_number = parts[0].split('.')[0][-2:]  # Gets the last two characters of the chromosome reference
    position_alleles = parts[1].split('g.')[1]
    match = re.match(r"(\d+)([A-Za-z]>[A-Za-z])", position_alleles)
    if match:
        position, alleles = match.groups()
    else:
        position, alleles = s, None
    alleles = alleles.split('>')
    ref_allele, alt_allele = alleles[0], alleles[-1]

    # Construct the eQTL format string
    gtex_format = f"chr{int(chrom_number)}_{position}_{ref_allele}_{alt_allele}"
    return gtex_format