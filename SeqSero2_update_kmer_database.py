#!/usr/bin/env python3

import argparse
import os,subprocess
import pickle

### SeqSero Kmer
def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = argparse.ArgumentParser(usage='Just type "SeqSero2_update_kmer_database.py", it will update kmer database automatically')
    return parser.parse_args()

def reverse_complement(sequence):
    complement = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N',
        'M': 'K',
        'R': 'Y',
        'W': 'W',
        'S': 'S',
        'Y': 'R',
        'K': 'M',
        'V': 'B',
        'H': 'D',
        'D': 'H',
        'B': 'V'
    }
    return "".join(complement[base] for base in reversed(sequence))

def multifasta_dict(multifasta):
    multifasta_list = [
        line.strip() for line in open(multifasta, 'r') if len(line.strip()) > 0
    ]
    headers = [i for i in multifasta_list if i[0] == '>']
    multifasta_dict = {}
    for h in headers:
        start = multifasta_list.index(h)
        for element in multifasta_list[start + 1:]:
            if element[0] == '>':
                break
            else:
                if h[1:] in multifasta_dict:
                    multifasta_dict[h[1:]] += element
                else:
                    multifasta_dict[h[1:]] = element
    return multifasta_dict

def createKmerDict_reads(list_of_strings, kmer):
    kmer_table = {}
    for string in list_of_strings:
        sequence = string.strip('\n')
        for i in range(len(sequence) - kmer + 1):
            new_mer = sequence[i:i + kmer].upper()
            new_mer_rc = reverse_complement(new_mer)
            if new_mer in kmer_table:
                kmer_table[new_mer.upper()] += 1
            else:
                kmer_table[new_mer.upper()] = 1
            if new_mer_rc in kmer_table:
                kmer_table[new_mer_rc.upper()] += 1
            else:
                kmer_table[new_mer_rc.upper()] = 1
    return kmer_table

def multifasta_to_kmers_dict(multifasta):
    multi_seq_dict = multifasta_dict(multifasta)
    lib_dict = {}
    for h in multi_seq_dict:
        lib_dict[h] = set(
            [k for k in createKmerDict_reads([multi_seq_dict[h]], 27)])
    return lib_dict

def get_salmid_invA_database(ex_dir):
  # read invA kmer and return it
  a = open(ex_dir + '/invA_mers_dict', 'rb')
  invA_dict = pickle.load(a)
  try:
    del invA_dict['version']
  except:
    pass
  return invA_dict

def get_salmid_rpoB_database(ex_dir):
  # read invA kmer and return it
  a = open(ex_dir + '/rpoB_mers_dict', 'rb')
  rpoB_dict = pickle.load(a)
  try:
    del rpoB_dict['version']
  except:
    pass
  return rpoB_dict

def main():
  args = parse_args()
  ex_dir = os.path.dirname(os.path.realpath(__file__))
  lib_dict = multifasta_to_kmers_dict(ex_dir + '/H_and_O_and_specific_genes.fasta')
  invA_dict=get_salmid_invA_database(ex_dir)
  #rpoB_dict=get_salmid_rpoB_database(ex_dir)
  lib_dict_new = lib_dict.copy()
  #print(len(lib_dict_new))
  lib_dict_new.update(invA_dict)
  #print(len(lib_dict_new))
  #lib_dict_new.update(rpoB_dict)
  #print(len(lib_dict_new))
  f = open(ex_dir + '/antigens.pickle', "wb")
  pickle.dump(lib_dict_new, f)
  f.close()

if __name__ == '__main__':
  main()
