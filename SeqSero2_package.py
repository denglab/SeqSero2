#!/usr/bin/env python3

import sys
import time
import random
import os
import subprocess
import gzip
import io
import pickle
import argparse
import itertools
from distutils.version import LooseVersion

### SeqSero Kmer
def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = argparse.ArgumentParser(usage='SeqSero2_package.py -t <data_type> -m <mode> -i <input_data> [-p <number of threads>] [-b <BWA_algorithm>]\n\nDevelopper: Shaokang Zhang (zskzsk@uga.edu), Hendrik C Den-Bakker (Hendrik.DenBakker@uga.edu) and Xiangyu Deng (xdeng@uga.edu)\n\nContact email:seqsero@gmail.com')#add "-m <data_type>" in future
    parser.add_argument("-i",nargs="+",help="<string>: path/to/input_data")
    parser.add_argument("-t",choices=['1','2','3','4','5','6'],help="<int>: '1'(pair-end reads, interleaved),'2'(pair-end reads, seperated),'3'(single-end reads), '4'(assembly),'5'(nanopore fasta),'6'(nanopore fastq)")
    parser.add_argument("-b",choices=['sam','mem'],default="mem",help="<string>: mode for mapping, 'sam'(bwa samse/sampe), 'mem'(bwa mem), default=mem") 
    parser.add_argument("-p",default="1",help="<int>: threads used for mapping mode, if p>4, only 4 threads will be used for assembly, default=1") 
    parser.add_argument("-m",choices=['k','a'],default="k",help="<string>: 'k'(kmer mode), 'a'(allele mode), default=k") 
    parser.add_argument("-d",help="<string>: output directory name, if not set, the output directory would be 'SeqSero_result_'+time stamp+one random number")
    parser.add_argument("-c",action="store_true",help="<flag>: if '-c' was flagged, SeqSero2 will use clean mode and only output serotyping prediction, the directory containing log files will be deleted")
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


def multifasta_single_string(multifasta):
    multifasta_list = [
        line.strip() for line in open(multifasta, 'r')
        if (len(line.strip()) > 0) and (line.strip()[0] != '>')
    ]
    return ''.join(multifasta_list)


def chunk_a_long_sequence(long_sequence, chunk_size=60):
    chunk_list = []
    steps = len(long_sequence) // 60  #how many chunks
    for i in range(steps):
        chunk_list.append(long_sequence[i * chunk_size:(i + 1) * chunk_size])
    chunk_list.append(long_sequence[steps * chunk_size:len(long_sequence)])
    return chunk_list


def target_multifasta_kmerizer(multifasta, k, kmerDict):
    forward_length = 300  #if find the target, put forward 300 bases
    reverse_length = 2200  #if find the target, put backward 2200 bases
    chunk_size = 60  #it will firstly chunk the single long sequence to multiple smaller sequences, it controls the size of those smaller sequences
    target_mers = []
    long_single_string = multifasta_single_string(multifasta)
    multifasta_list = chunk_a_long_sequence(long_single_string, chunk_size)
    unit_length = len(multifasta_list[0])
    forward_lines = int(forward_length / unit_length) + 1
    reverse_lines = int(forward_length / unit_length) + 1
    start_num = 0
    end_num = 0
    for i in range(len(multifasta_list)):
        if i not in range(start_num, end_num):  #avoid computational repetition
            line = multifasta_list[i]
            start = int((len(line) - k) // 2)
            s1 = line[start:k + start]
            if s1 in kmerDict:  #detect it is a potential read or not (use the middle part)
                if i - forward_lines >= 0:
                    start_num = i - forward_lines
                else:
                    start_num = 0
                if i + reverse_lines <= len(multifasta_list) - 1:
                    end_num = i + reverse_lines
                else:
                    end_num = len(multifasta_list) - 1
                target_list = [
                    x.strip() for x in multifasta_list[start_num:end_num]
                ]
                target_line = "".join(target_list)
                target_mers += [
                    k1 for k1 in createKmerDict_reads([str(target_line)], k)
                ]  ##changed k to k1, just want to avoid the mixes of this "k" (kmer) to the "k" above (kmer length)
        else:
            pass
    return set(target_mers)


def target_read_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = []
    if file.endswith(".gz"):
        file_content = io.BufferedReader(gzip.open(file))
    else:
        file_content = open(file, "r").readlines()
    for line in file_content:
        start = int((len(line) - k) // 2)
        if i % 4 == 2:
            if file.endswith(".gz"):
                s1 = line[start:k + start].decode()
                line = line.decode()
            else:
                s1 = line[start:k + start]
            if s1 in kmerDict:  #detect it is a potential read or not (use the middle part)
                n_reads += 1
                total_coverage += len(line)
                target_mers += [
                    k1 for k1 in createKmerDict_reads([str(line)], k)
                ]  #changed k to k1, just want to avoid the mixes of this "k" (kmer) to the "k" above (kmer length)
        i += 1
        if total_coverage >= 4000000:
            break
    return set(target_mers)


def minion_fasta_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = {}
    for line in open(file):
        if i % 2 == 0:
            for kmer, rc_kmer in kmers(line.strip().upper(), k):
                if (kmer in kmerDict) or (rc_kmer in kmerDict):
                    if kmer in target_mers:
                        target_mers[kmer] += 1
                    else:
                        target_mers[kmer] = 1
                    if rc_kmer in target_mers:
                        target_mers[rc_kmer] += 1
                    else:
                        target_mers[rc_kmer] = 1
        i += 1
    return set([h for h in target_mers])


def minion_fastq_kmerizer(file, k, kmerDict):
    i = 1
    n_reads = 0
    total_coverage = 0
    target_mers = {}
    for line in open(file):
        if i % 4 == 2:
            for kmer, rc_kmer in kmers(line.strip().upper(), k):
                if (kmer in kmerDict) or (rc_kmer in kmerDict):
                    if kmer in target_mers:
                        target_mers[kmer] += 1
                    else:
                        target_mers[kmer] = 1
                    if rc_kmer in target_mers:
                        target_mers[rc_kmer] += 1
                    else:
                        target_mers[rc_kmer] = 1
        i += 1
    return set([h for h in target_mers])


def multifasta_single_string2(multifasta):
    single_string = ''
    with open(multifasta, 'r') as f:
        for line in f:
            if line.strip()[0] == '>':
                pass
            else:
                single_string += line.strip()
    return single_string


def kmers(seq, k):
    rev_comp = reverse_complement(seq)
    for start in range(1, len(seq) - k + 1):
        yield seq[start:start + k], rev_comp[-(start + k):-start]


def multifasta_to_kmers_dict(multifasta,k_size):#used to create database kmer set
    multi_seq_dict = multifasta_dict(multifasta)
    lib_dict = {}
    for h in multi_seq_dict:
        lib_dict[h] = set(
            [k for k in createKmerDict_reads([multi_seq_dict[h]], k_size)])
    return lib_dict


def Combine(b, c):
    fliC_combinations = []
    fliC_combinations.append(",".join(c))
    temp_combinations = []
    for i in range(len(b)):
        for x in itertools.combinations(b, i + 1):
            temp_combinations.append(",".join(x))
    for x in temp_combinations:
        temp = []
        for y in c:
            temp.append(y)
        temp.append(x)
        temp = ",".join(temp)
        temp = temp.split(",")
        temp.sort()
        temp = ",".join(temp)
        fliC_combinations.append(temp)
    return fliC_combinations


def seqsero_from_formula_to_serotypes(Otype, fliC, fljB, special_gene_list,subspecies):
    #like test_output_06012017.txt
    #can add more varialbles like sdf-type, sub-species-type in future (we can conclude it into a special-gene-list)
    from Initial_Conditions import phase1
    from Initial_Conditions import phase2
    from Initial_Conditions import phaseO
    from Initial_Conditions import sero
    from Initial_Conditions import subs
    seronames = []
    seronames_none_subspecies=[]
    for i in range(len(phase1)):
        fliC_combine = []
        fljB_combine = []
        if phaseO[i] == Otype: # no VII in KW, but it's there
            ### for fliC, detect every possible combinations to avoid the effect of "["
            if phase1[i].count("[") == 0:
                fliC_combine.append(phase1[i])
            elif phase1[i].count("[") >= 1:
                c = []
                b = []
                if phase1[i][0] == "[" and phase1[i][-1] == "]" and phase1[i].count(
                        "[") == 1:
                    content = phase1[i].replace("[", "").replace("]", "")
                    fliC_combine.append(content)
                    fliC_combine.append("-")
                else:
                    for x in phase1[i].split(","):
                        if "[" in x:
                            b.append(x.replace("[", "").replace("]", ""))
                        else:
                            c.append(x)
                    fliC_combine = Combine(
                        b, c
                    )  #Combine will offer every possible combinations of the formula, like f,[g],t: f,t  f,g,t
            ### end of fliC "[" detect
            ### for fljB, detect every possible combinations to avoid the effect of "["
            if phase2[i].count("[") == 0:
                fljB_combine.append(phase2[i])
            elif phase2[i].count("[") >= 1:
                d = []
                e = []
                if phase2[i][0] == "[" and phase2[i][-1] == "]" and phase2[i].count(
                        "[") == 1:
                    content = phase2[i].replace("[", "").replace("]", "")
                    fljB_combine.append(content)
                    fljB_combine.append("-")
                else:
                    for x in phase2[i].split(","):
                        if "[" in x:
                            d.append(x.replace("[", "").replace("]", ""))
                        else:
                            e.append(x)
                    fljB_combine = Combine(d, e)
            ### end of fljB "[" detect
            new_fliC = fliC.split(
                ","
            )  #because some antigen like r,[i] not follow alphabetical order, so use this one to judge and can avoid missings
            new_fliC.sort()
            new_fliC = ",".join(new_fliC)
            new_fljB = fljB.split(",")
            new_fljB.sort()
            new_fljB = ",".join(new_fljB)
            if (new_fliC in fliC_combine
                    or fliC in fliC_combine) and (new_fljB in fljB_combine
                                                  or fljB in fljB_combine):
                if subs[i] == subspecies:
                  seronames.append(sero[i])
                seronames_none_subspecies.append(sero[i])
    #analyze seronames
    subspecies_pointer=""
    if len(seronames) == 0 and len(seronames_none_subspecies)!=0:
      seronames=seronames_none_subspecies
      subspecies_pointer="1"
    if len(seronames) == 0:
        seronames = [
            "N/A (The predicted antigenic profile does not exist in the White-Kauffmann-Le Minor scheme)"
        ]
    star = ""
    star_line = ""
    if len(seronames) > 1:  #there are two possible predictions for serotypes
        star = "*"
        star_line = "The predicted serotypes share the same general formula:\t" + Otype + ":" + fliC + ":" + fljB + "\n"
    if subspecies_pointer=="1" and len(seronames_none_subspecies)!=0:
      star="*"
      star_line="The formula with this subspieces prediction can't get a serotype in KW manual, and the serotyping prediction was made without considering it."+star_line
    if  Otype=="":
      Otype="-"
    predict_form = Otype + ":" + fliC + ":" + fljB
    predict_sero = (" or ").join(seronames)
    ###special test for Enteritidis
    if predict_form == "9:g,m:-":
        sdf = "-"
        for x in special_gene_list:
            if x.startswith("sdf"):
                sdf = "+"
        predict_form = predict_form + " Sdf prediction:" + sdf
        if sdf == "-":
            star = "*"
            star_line = "Additional characterization is necessary to assign a serotype to this strain.  Commonly circulating strains of serotype Enteritidis are sdf+, although sdf- strains of serotype Enteritidis are known to exist. Serotype Gallinarum is typically sdf- but should be quite rare. Sdf- strains of serotype Enteritidis and serotype Gallinarum can be differentiated by phenotypic profile or genetic criteria.\n"
            predict_sero = "Gallinarum/Enteritidis sdf -"
    ###end of special test for Enteritidis
    elif predict_form == "4:i:-":
        predict_sero = "potential monophasic variant of Typhimurium"
    elif predict_form == "4:r:-":
        predict_sero = "potential monophasic variant of Heidelberg"
    elif predict_form == "4:b:-":
        predict_sero = "potential monophasic variant of Paratyphi B"
    elif predict_form == "8:e,h:1,2":
        predict_sero = "Newport"
        star = "*"
        star_line = "Serotype Bardo shares the same antigenic profile with Newport, but Bardo is exceedingly rare."
    claim = "The serotype(s) is/are the only serotype(s) with the indicated antigenic profile currently recognized in the Kauffmann White Scheme.  New serotypes can emerge and the possibility exists that this antigenic profile may emerge in a different subspecies.  Identification of strains to the subspecies level should accompany serotype determination; the same antigenic profile in different subspecies is considered different serotypes.\n"
    if "N/A" in predict_sero:
        claim = ""
    #special test for Typhimurium
    if "Typhimurium" in predict_sero or predict_form == "4:i:-":
        normal = 0
        mutation = 0
        for x in special_gene_list:
            if "oafA-O-4_full" in x:
                normal = float(special_gene_list[x])
            elif "oafA-O-4_5-" in x:
                mutation = float(special_gene_list[x])
        if normal > mutation:
            pass
        elif normal < mutation:
            predict_sero = predict_sero.strip() + "(O5-)"
            star = "*"
            star_line = "Detected the deletion of O5-."
        else:
            pass
    #special test for Paratyphi B
    if "Paratyphi B" in predict_sero or predict_form == "4:b:-":
        normal = 0
        mutation = 0
        for x in special_gene_list:
            if "gntR-family-regulatory-protein_dt-positive" in x:
                normal = float(special_gene_list[x])
            elif "gntR-family-regulatory-protein_dt-negative" in x:
                mutation = float(special_gene_list[x])
        #print(normal,mutation)
        if normal > mutation:
            predict_sero = predict_sero.strip() + "(dt+)"
            star = "*"
            star_line = "Didn't detect the SNP for dt- which means this isolate is a Paratyphi B variant L(+) tartrate(+)."
        elif normal < mutation:
            predict_sero = predict_sero.strip() + "(dt-)"
            star = "*"
            star_line = "Detected the SNP for dt- which means this isolate is a systemic pathovar of Paratyphi B."
        else:
            star = "*"
            star_line = "Failed to detect the SNP for dt-, can't decide it's a Paratyphi B variant L(+) tartrate(+) or not."
    #special test for O13,22 and O13,23
    if Otype=="13":
      ex_dir = os.path.dirname(os.path.realpath(__file__))
      f = open(ex_dir + '/special.pickle', 'rb')
      special = pickle.load(f)
      O22_O23=special['O22_O23']
      if predict_sero.split(" or ")[0] in O22_O23[-1]:
        O22_score=0
        O23_score=0
        for x in special_gene_list:
            if "O:22" in x:
                O22_score = O22_score+float(special_gene_list[x])
            elif "O:23" in x:
                O23_score = O23_score+float(special_gene_list[x])
        #print(O22_score,O23_score)
        for z in O22_O23[0]:
          if predict_sero.split(" or ")[0] in z:
            if O22_score > O23_score:
              star = "*"
              star_line = "Detected O22 specific genes to further differenciate '"+predict_sero+"'."
              predict_sero = z[0]
            elif O22_score < O23_score:
              star = "*"
              star_line = "Detected O23 specific genes to further differenciate '"+predict_sero+"'."
              predict_sero = z[1]
            else:
              star = "*"
              star_line = "Fail to detect O22 and O23 differences."
    #special test for O6,8 
    merge_O68_list=["Blockley","Bovismorbificans","Hadar","Litchfield","Manhattan","Muenchen"]
    for x in merge_O68_list:
      if x in predict_sero:
        predict_sero=x
        star=""
        star_line=""
    #special test for Montevideo; most of them are monophasic
    if "Montevideo" in predict_sero and "1,2,7" in predict_form:
      star="*"
      star_line="Montevideo is almost always monophasic, having an antigen called for the fljB position may be a result of Salmonella-Salmonella contamination."
    return predict_form, predict_sero, star, star_line, claim
### End of SeqSero Kmer part

### Begin of SeqSero2 allele prediction and output
def xml_parse_score_comparision_seqsero(xmlfile):
  #used to do seqsero xml analysis
  from Bio.Blast import NCBIXML
  handle=open(xmlfile)
  handle=NCBIXML.parse(handle)
  handle=list(handle)
  List=[]
  List_score=[]
  List_ids=[]
  List_query_region=[]
  for i in range(len(handle)):
    if len(handle[i].alignments)>0:
      for j in range(len(handle[i].alignments)):
        score=0
        ids=0
        cover_region=set() #fixed problem that repeated calculation leading percentage > 1
        List.append(handle[i].query.strip()+"___"+handle[i].alignments[j].hit_def)
        for z in range(len(handle[i].alignments[j].hsps)):
          hsp=handle[i].alignments[j].hsps[z]
          temp=set(range(hsp.query_start,hsp.query_end))
          if len(cover_region)==0:
            cover_region=cover_region|temp
            fraction=1
          else:
            fraction=1-len(cover_region&temp)/float(len(temp))
            cover_region=cover_region|temp
          if "last" in handle[i].query or "first" in handle[i].query:
            score+=hsp.bits*fraction
            ids+=float(hsp.identities)/handle[i].query_length*fraction
          else:
            score+=hsp.bits*fraction
            ids+=float(hsp.identities)/handle[i].query_length*fraction
        List_score.append(score)
        List_ids.append(ids)
        List_query_region.append(cover_region)
  temp=zip(List,List_score,List_ids,List_query_region)
  Final_list=sorted(temp, key=lambda d:d[1], reverse = True)
  return Final_list


def Uniq(L,sort_on_fre="none"): #return the uniq list and the count number
  Old=L
  L.sort()
  L = [L[i] for i in range(len(L)) if L[i] not in L[:i]]
  count=[]
  for j in range(len(L)):
    y=0
    for x in Old:
      if L[j]==x:
        y+=1
    count.append(y)
  if sort_on_fre!="none":
    d=zip(*sorted(zip(count, L)))
    L=d[1]
    count=d[0]
  return (L,count)

def judge_fliC_or_fljB_from_head_tail_for_one_contig(nodes_vs_score_list):
  #used to predict it's fliC or fljB for one contig, based on tail and head score, but output the score difference,if it is very small, then not reliable, use blast score for whole contig to test
  #this is mainly used for 
  a=nodes_vs_score_list
  fliC_score=0
  fljB_score=0
  for z in a:
    if "fliC" in z[0]:
      fliC_score+=z[1]
    elif "fljB" in z[0]:
      fljB_score+=z[1]
  if fliC_score>=fljB_score:
    role="fliC"
  else:
    role="fljB"
  return (role,abs(fliC_score-fljB_score))

def judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(node_name,Final_list,Final_list_passed):
  #used to predict contig is fliC or fljB, if the differnce score value on above head_and_tail is less than 10 (quite small)
  #also used when no head or tail got blasted score for the contig
  role=""
  for z in Final_list_passed:
    if node_name in z[0]:
      role=z[0].split("_")[0]
      break
  return role

def fliC_or_fljB_judge_from_head_tail_sequence(nodes_list,tail_head_list,Final_list,Final_list_passed):
  #nodes_list is the c created by c,d=Uniq(nodes) in below function
  first_target=""
  role_list=[]
  for x in nodes_list:
    a=[]
    role=""
    for y in tail_head_list:
      if x in y[0]:
        a.append(y)
    if len(a)==4:
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list,Final_list_passed)
    elif len(a)==3:
      ###however, if the one with highest score is the fewer one, compare their accumulation score
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list,Final_list_passed)
      ###end of above score comparison
    elif len(a)==2:
      #must on same node, if not, then decide with unit blast score, blast-score/length_of_special_sequence(30 or 37)
      temp=[]
      for z in a:
        temp.append(z[0].split("_")[0])
      m,n=Uniq(temp)#should only have one choice, but weird situation might occur too
      if len(m)==1:
        pass
      else:
        pass
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list,Final_list_passed)
        ###need to desgin a algorithm to guess most possible situation for nodes_list, See the situations of test evaluation
    elif len(a)==1:
      #that one
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list,Final_list_passed)
      #need to evaluate, in future, may set up a cut-off, if not met, then just find Final_list_passed best match,like when "a==0"
    else:#a==0
      #use Final_list_passed best match
      for z in Final_list_passed:
        if x in z[0]:
          role=z[0].split("_")[0]
          break
    #print x,role,len(a)
    role_list.append((role,x))
  if len(role_list)==2:
    if role_list[0][0]==role_list[1][0]:#this is the most cocmmon error, two antigen were assigned to same phase
      #just use score to do a final test
      role_list=[]
      for x in nodes_list: 
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list,Final_list_passed)
        role_list.append((role,x))
  return role_list

def decide_contig_roles_for_H_antigen(Final_list,Final_list_passed):
  #used to decide which contig is FliC and which one is fljB
  contigs=[]
  nodes=[]
  for x in Final_list_passed:
    if x[0].startswith("fl") and "last" not in x[0] and "first" not in x[0]:
      nodes.append(x[0].split("___")[1].strip())
  c,d=Uniq(nodes)#c is node_list
  #print c
  tail_head_list=[x for x in Final_list if ("last" in x[0] or "first" in x[0])]
  roles=fliC_or_fljB_judge_from_head_tail_sequence(c,tail_head_list,Final_list,Final_list_passed)
  return roles

def decide_O_type_and_get_special_genes(Final_list,Final_list_passed):
  #decide O based on Final_list
  O_choice="?"
  O_list=[]
  special_genes={}
  nodes=[]
  for x in Final_list_passed:
    if x[0].startswith("O-"):
      nodes.append(x[0].split("___")[1].strip())
    elif not x[0].startswith("fl"):
      special_genes[x[0]]=x[2]#08172018, x[2] changed from x[-1]
  #print "special_genes:",special_genes
  c,d=Uniq(nodes)
  #print "potential O antigen contig",c
  final_O=[]
  O_nodes_list=[]
  for x in c:#c is the list for contigs
    temp=0
    for y in Final_list_passed:
      if x in y[0] and y[0].startswith("O-"):
        final_O.append(y)
        break
  ### O contig has the problem of two genes on same contig, so do additional test
  potenial_new_gene=""
  for x in final_O:
    pointer=0 #for genes merged or not
    #not consider O-1,3,19_not_in_3,10, too short compared with others
    if "O-1,3,19_not_in_3,10" not in x[0] and int(x[0].split("__")[1].split("___")[0])*x[2]+850 <= int(x[0].split("length_")[1].split("_")[0]):#gene length << contig length; for now give 300*2 (for secureity can use 400*2) as flank region
      pointer=x[0].split("___")[1].strip()#store the contig name
      print(pointer)
    if pointer!=0:#it has potential merge event
      for y in Final_list:
        if pointer in y[0] and y not in final_O and (y[1]>=int(y[0].split("__")[1].split("___")[0])*1.5 or (y[1]>=int(y[0].split("__")[1].split("___")[0])*y[2] and y[1]>=400)):#that's a realtively strict filter now; if passed, it has merge event and add one more to final_O
          potenial_new_gene=y
          #print(potenial_new_gene)
          break
  if potenial_new_gene!="":
    print("two differnt genes in same contig, fix it for O antigen")
    print(potenial_new_gene[:3])
    final_O.append(potenial_new_gene)
  ### end of the two genes on same contig test
  final_O=sorted(final_O,key=lambda x: x[2], reverse=True)#sorted
  if len(final_O)==0 or (len(final_O)==1 and "O-1,3,19_not_in_3,10" in final_O[0][0]):
    #print "$$$No Otype, due to no hit"#may need to be changed
    O_choice="-"
  else:
    highest_O_coverage=max([float(x[0].split("_cov_")[-1]) for x in final_O if "O-1,3,19_not_in_3,10" not in x[0]])
    O_list=[]
    O_list_less_contamination=[]
    for x in final_O:
      if not "O-1,3,19_not_in_3,10__130" in x[0]:#O-1,3,19_not_in_3,10 is too small, which may affect further analysis; to avoid contamination affect, use 0.15 of highest coverage as cut-off
        O_list.append(x[0].split("__")[0])
        O_nodes_list.append(x[0].split("___")[1])
        if float(x[0].split("_cov_")[-1])>highest_O_coverage*0.15:
          O_list_less_contamination.append(x[0].split("__")[0])
    ### special test for O9,46 and O3,10 family
    if ("O-9,46_wbaV" in O_list or "O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254" in O_list) and O_list_less_contamination[0].startswith("O-9,"):#not sure should use and float(O9_wbaV)/float(num_1) > 0.1
      if "O-9,46_wzy" in O_list:#and float(O946_wzy)/float(num_1) > 0.1
        O_choice="O-9,46"
        #print "$$$Most possilble Otype:  O-9,46"
      elif "O-9,46,27_partial_wzy" in O_list:#and float(O94627)/float(num_1) > 0.1
        O_choice="O-9,46,27"
        #print "$$$Most possilble Otype:  O-9,46,27"
      else:
        O_choice="O-9"#next, detect O9 vs O2?
        O2=0
        O9=0
        for z in special_genes:
          if "tyr-O-9" in z:
            O9=special_genes[z]
          elif "tyr-O-2" in z:
            O2=special_genes[z]
        if O2>O9:
          O_choice="O-2"
        elif O2<O9:
          pass
        else:
          pass
          #print "$$$No suitable one, because can't distinct it's O-9 or O-2, but O-9 has a more possibility."
    elif ("O-3,10_wzx" in O_list) and ("O-9,46_wzy" in O_list) and (O_list[0].startswith("O-3,10") or O_list_less_contamination[0].startswith("O-9,46_wzy")):#and float(O310_wzx)/float(num_1) > 0.1 and float(O946_wzy)/float(num_1) > 0.1
      if "O-3,10_not_in_1,3,19" in O_list:#and float(O310_no_1319)/float(num_1) > 0.1
        O_choice="O-3,10"
        #print "$$$Most possilble Otype:  O-3,10 (contain O-3,10_not_in_1,3,19)"
      else:
        O_choice="O-1,3,19"
        #print "$$$Most possilble Otype:  O-1,3,19 (not contain O-3,10_not_in_1,3,19)"
    ### end of special test for O9,46 and O3,10 family
    else:
      try: 
        max_score=0
        for x in final_O:
          if x[2]>=max_score and float(x[0].split("_cov_")[-1])>highest_O_coverage*0.15:#use x[2],08172018, the "coverage identity = cover_length * identity"; also meet coverage threshold
            max_score=x[2]#change from x[-1] to x[2],08172018
            O_choice=x[0].split("_")[0]
        if O_choice=="O-1,3,19":
          O_choice=final_O[1][0].split("_")[0]
        #print "$$$Most possilble Otype: ",O_choice
      except:
        pass
        #print "$$$No suitable Otype, or failure of mapping (please check the quality of raw reads)"
  #print "O:",O_choice,O_nodes_list
  Otypes=[]
  for x in O_list:
    if x!="O-1,3,19_not_in_3,10":
      if "O-9,46_" not in x:
        Otypes.append(x.split("_")[0])
      else:
        Otypes.append(x.split("-from")[0])#O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254
  #Otypes=[x.split("_")[0] for x in O_list if x!="O-1,3,19_not_in_3,10"]
  Otypes_uniq,Otypes_fre=Uniq(Otypes)
  contamination_O=""
  if O_choice=="O-9,46,27" or O_choice=="O-3,10" or O_choice=="O-1,3,19":
    if len(Otypes_uniq)>2:
      contamination_O="potential contamination from O antigen signals"
  else:
    if len(Otypes_uniq)>1:
      if O_choice=="O-4" and len(Otypes_uniq)==2 and "O-9,46,27" in Otypes_uniq: #for special 4,12,27 case such as Bredeney and Schwarzengrund
        contamination_O=""
      elif O_choice=="O-9,46" and len(Otypes_uniq)==2 and "O-9,46_wbaV" in Otypes_uniq and "O-9,46_wzy" in Otypes_uniq: #for special 4,12,27 case such as Bredeney and Schwarzengrund
        contamination_O=""
      else:
        contamination_O="potential contamination from O antigen signals"
  return O_choice,O_nodes_list,special_genes,final_O,contamination_O
### End of SeqSero2 allele prediction and output

def get_input_files(make_dir,input_file,data_type,dirpath):
  #tell input files from datatype
  #"<int>: '1'(pair-end reads, interleaved),'2'(pair-end reads, seperated),'3'(single-end reads), '4'(assembly),'5'(nanopore fasta),'6'(nanopore fastq)"
  for_fq=""
  rev_fq=""
  os.chdir(make_dir)
  if data_type=="1":
    input_file=input_file[0].split("/")[-1]
    if input_file.endswith(".sra"):
      subprocess.check_call("fastq-dump --split-files "+input_file,shell=True)
      for_fq=input_file.replace(".sra","_1.fastq")
      rev_fq=input_file.replace(".sra","_2.fastq")
    else:
      core_id=input_file.split(".fastq")[0].split(".fq")[0]
      for_fq=core_id+"_1.fastq"
      rev_fq=core_id+"_2.fastq"
      if input_file.endswith(".gz"):
        subprocess.check_call("gzip -dc "+input_file+" | "+dirpath+"/deinterleave_fastq.sh "+for_fq+" "+rev_fq,shell=True)
      else:
        subprocess.check_call("cat "+input_file+" | "+dirpath+"/deinterleave_fastq.sh "+for_fq+" "+rev_fq,shell=True)
  elif data_type=="2":
    for_fq=input_file[0].split("/")[-1]
    rev_fq=input_file[1].split("/")[-1]
  elif data_type=="3":
    input_file=input_file[0].split("/")[-1]
    if input_file.endswith(".sra"):
      subprocess.check_call("fastq-dump --split-files "+input_file,shell=True)
      for_fq=input_file.replace(".sra","_1.fastq")
    else:
      for_fq=input_file
  elif data_type in ["4","5","6"]:
    for_fq=input_file[0].split("/")[-1]
  os.chdir("..")
  return for_fq,rev_fq

def predict_O_and_H_types(Final_list,Final_list_passed):
  #get O and H types from Final_list from blast parsing; allele mode
  fliC_choice="-"
  fljB_choice="-"
  fliC_contig="NA"
  fljB_contig="NA"
  fliC_region=set([0])
  fljB_region=set([0,])
  fliC_length=0 #can be changed to coverage in future
  fljB_length=0 #can be changed to coverage in future
  O_choice="-"#no need to decide O contig for now, should be only one
  O_choice,O_nodes,special_gene_list,O_nodes_roles,contamination_O=decide_O_type_and_get_special_genes(Final_list,Final_list_passed)#decide the O antigen type and also return special-gene-list for further identification
  O_choice=O_choice.split("-")[-1].strip()
  if (O_choice=="1,3,19" and len(O_nodes_roles)==1 and "1,3,19" in O_nodes_roles[0][0]) or O_choice=="":
    O_choice="-"
  H_contig_roles=decide_contig_roles_for_H_antigen(Final_list,Final_list_passed)#decide the H antigen contig is fliC or fljB
  log_file=open("SeqSero_log.txt","a")
  print("O_contigs:")
  log_file.write("O_contigs:\n")
  for x in O_nodes_roles:
    if "O-1,3,19_not_in_3,10" not in x[0]:#O-1,3,19_not_in_3,10 is just a small size marker
      print(x[0].split("___")[-1],x[0].split("__")[0],"blast score:",x[1],"identity%:",str(round(x[2]*100,2))+"%",str(min(x[-1]))+" to "+str(max(x[-1])))
      log_file.write(x[0].split("___")[-1]+" "+x[0].split("__")[0]+" "+"blast score: "+str(x[1])+"identity%:"+str(round(x[2]*100,2))+"% "+str(min(x[-1]))+" to "+str(max(x[-1]))+"\n")
  if len(H_contig_roles)!=0:
    highest_H_coverage=max([float(x[1].split("_cov_")[-1]) for x in H_contig_roles]) #less than highest*0.1 would be regarded as contamination and noises, they will still be considered in contamination detection and logs, but not used as final serotype output
  else:
    highest_H_coverage=0
  for x in H_contig_roles:
    #if multiple choices, temporately select the one with longest length for now, will revise in further change
    if "fliC" == x[0] and int(x[1].split("_")[3])>=fliC_length and x[1] not in O_nodes and float(x[1].split("_cov_")[-1])>highest_H_coverage*0.13:#remember to avoid the effect of O-type contig, so should not in O_node list
      fliC_contig=x[1]
      fliC_length=int(x[1].split("_")[3])
    elif "fljB" == x[0] and int(x[1].split("_")[3])>=fljB_length and x[1] not in O_nodes and float(x[1].split("_cov_")[-1])>highest_H_coverage*0.13:
      fljB_contig=x[1]
      fljB_length=int(x[1].split("_")[3])
  for x in Final_list_passed:
    if fliC_choice=="-" and "fliC_" in x[0] and fliC_contig in x[0]:
      fliC_choice=x[0].split("_")[1]
    elif fljB_choice=="-" and "fljB_" in x[0] and fljB_contig in x[0]:
      fljB_choice=x[0].split("_")[1]
    elif fliC_choice!="-" and fljB_choice!="-":
      break
  #now remove contigs not in middle core part
  first_allele="NA"
  first_allele_percentage=0
  for x in Final_list:
    if x[0].startswith("fliC") or x[0].startswith("fljB"):
      first_allele=x[0].split("__")[0] #used to filter those un-middle contigs
      first_allele_percentage=x[2]
      break 
  additional_contigs=[]
  for x in Final_list:
    if first_allele in x[0]:
      if (fliC_contig == x[0].split("___")[-1]): 
        fliC_region=x[3]
      elif fljB_contig!="NA" and (fljB_contig == x[0].split("___")[-1]):
        fljB_region=x[3]
      else:
        if x[1]*1.1>int(x[0].split("___")[1].split("_")[3]):#loose threshold by multiplying 1.1
          additional_contigs.append(x)
        #else:
          #print x[:3]
  #we can just use the fljB region (or fliC depends on size), no matter set() or contain a large locations (without middle part); however, if none of them is fully assembled, use 500 and 1200 as conservative cut-off
  if first_allele_percentage>0.9:
    if len(fliC_region)>len(fljB_region) and (max(fljB_region)-min(fljB_region))>1000:
      target_region=fljB_region|(fliC_region-set(range(min(fljB_region),max(fljB_region)))) #fljB_region|(fliC_region-set(range(min(fljB_region),max(fljB_region))))
    elif len(fliC_region)<len(fljB_region) and (max(fliC_region)-min(fliC_region))>1000:
      target_region=fliC_region|(fljB_region-set(range(min(fliC_region),max(fliC_region))))  #fljB_region|(fliC_region-set(range(min(fljB_region),max(fljB_region))))
    else:
      target_region=set()#doesn't do anything
  else:
    target_region=set()#doesn't do anything
  #print(target_region)
  #print(additional_contigs)
  target_region2=set(list(range(0,525))+list(range(1200,1700)))#I found to use 500 to 1200 as special region would be best
  target_region=target_region2|target_region
  for x in additional_contigs:
    removal=0
    contig_length=int(x[0].split("___")[1].split("length_")[-1].split("_")[0])
    if fljB_contig not in x[0] and fliC_contig not in x[0] and len(target_region&x[3])/float(len(x[3]))>0.65 and contig_length*0.5<len(x[3])<contig_length*1.5: #consider length and alignment length for now, but very loose,0.5 and 1.5 as cut-off
      removal=1
    else:
      if first_allele_percentage > 0.9 and float(x[0].split("__")[1].split("___")[0])*x[2]/len(x[-1])>0.96:#if high similiarity with middle part of first allele (first allele >0.9, already cover middle part)
        removal=1
      else:
        pass
    if removal==1:
      for y in H_contig_roles:
        if y[1] in x[0]:
          H_contig_roles.remove(y)
    else:
      pass
      #print(x[:3],contig_length,len(target_region&x[3])/float(len(x[3])),contig_length*0.5,len(x[3]),contig_length*1.5)
  #end of removing none-middle contigs
  print("H_contigs:")
  log_file.write("H_contigs:\n")
  H_contig_stat=[]
  H1_cont_stat={}
  H2_cont_stat={}
  for i in range(len(H_contig_roles)):
    x=H_contig_roles[i]
    a=0
    for y in Final_list_passed:
      if x[1] in y[0] and y[0].startswith(x[0]):
        if "first" in y[0] or "last" in y[0]: #this is the final filter to decide it's fliC or fljB, if can't pass, then can't decide
          for y in Final_list_passed: #it's impossible to has the "first" and "last" allele as prediction, so re-do it
            if x[1] in y[0]:#it's very possible to be third phase allele, so no need to make it must be fliC or fljB
              print(x[1],"can't_decide_fliC_or_fljB",y[0].split("_")[1],"blast_score:",y[1],"identity%:",str(round(y[2]*100,2))+"%",str(min(y[-1]))+" to "+str(max(y[-1])))
              log_file.write(x[1]+" "+x[0]+" "+y[0].split("_")[1]+" "+"blast_score: "+str(y[1])+" identity%:"+str(round(y[2]*100,2))+"% "+str(min(y[-1]))+" to "+str(max(y[-1]))+"\n")
              H_contig_roles[i]="can't decide fliC or fljB, may be third phase"
              break
        else:
          print(x[1],x[0],y[0].split("_")[1],"blast_score:",y[1],"identity%:",str(round(y[2]*100,2))+"%",str(min(y[-1]))+" to "+str(max(y[-1])))
          log_file.write(x[1]+" "+x[0]+" "+y[0].split("_")[1]+" "+"blast_score: "+str(y[1])+" identity%:"+str(round(y[2]*100,2))+"% "+str(min(y[-1]))+" to "+str(max(y[-1]))+"\n")
        if x[0]=="fliC":
          if y[0].split("_")[1] not in H1_cont_stat:
            H1_cont_stat[y[0].split("_")[1]]=y[2]
          else:
            H1_cont_stat[y[0].split("_")[1]]+=y[2]
        if x[0]=="fljB":
          if y[0].split("_")[1] not in H2_cont_stat:
            H2_cont_stat[y[0].split("_")[1]]=y[2]
          else:
            H2_cont_stat[y[0].split("_")[1]]+=y[2]
        break
  #detect contaminations
  #print(H1_cont_stat)
  #print(H2_cont_stat)
  H1_cont_stat_list=[x for x in H1_cont_stat if H1_cont_stat[x]>0.2]
  H2_cont_stat_list=[x for x in H2_cont_stat if H2_cont_stat[x]>0.2]
  contamination_H=""
  if len(H1_cont_stat_list)>1 or len(H2_cont_stat_list)>1:
    contamination_H="potential contamination from H antigen signals"
  elif len(H2_cont_stat_list)==1 and fljB_contig=="NA":
    contamination_H="potential contamination from H antigen signals, uncommon weak fljB signals detected"
  print(contamination_O)
  print(contamination_H)
  log_file.write(contamination_O+"\n")
  log_file.write(contamination_H+"\n")
  log_file.close()
  return O_choice,fliC_choice,fljB_choice,special_gene_list,contamination_O,contamination_H

def get_input_K(input_file,lib_dict,data_type,k_size):
  #kmer mode; get input_Ks from dict and data_type
  kmers = []
  for h in lib_dict:
      kmers += lib_dict[h]
  if data_type == '4':
      input_Ks = target_multifasta_kmerizer(input_file, k_size, set(kmers))
  elif data_type == '1' or data_type == '2' or data_type == '3':#set it for now, will change later
      input_Ks = target_read_kmerizer(input_file, k_size, set(kmers))
  elif data_type == '5':#minion_2d_fasta
      input_Ks = minion_fasta_kmerizer(input_file, k_size, set(kmers))
  if data_type == '6':#minion_2d_fastq
      input_Ks = minion_fastq_kmerizer(input_file, k_size, set(kmers))
  return input_Ks

def get_kmer_dict(lib_dict,input_Ks):
  #kmer mode; get predicted types
  O_dict = {}
  H_dict = {}
  Special_dict = {}
  for h in lib_dict:
      score = (len(lib_dict[h] & input_Ks) / len(lib_dict[h])) * 100
      if score > 1:  # Arbitrary cut-off for similarity score very low but seems necessary to detect O-3,10 in some cases
          if h.startswith('O-') and score > 25:
              O_dict[h] = score
          if h.startswith('fl') and score > 40:
              H_dict[h] = score
          if (h[:2] != 'fl') and (h[:2] != 'O-'):
              Special_dict[h] = score
  return O_dict,H_dict,Special_dict

def call_O_and_H_type(O_dict,H_dict,Special_dict,make_dir):
  log_file=open("SeqSero_log.txt","a")
  log_file.write("O_scores:\n")
  #call O:
  highest_O = '-'
  if len(O_dict) == 0:
      pass
  else:
      for x in O_dict:
          log_file.write(x+"\t"+str(O_dict[x])+"\n")
      if ('O-9,46_wbaV__1002' in O_dict and O_dict['O-9,46_wbaV__1002']>70) or ("O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254__1002" in O_dict and O_dict['O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254__1002']>70):  # not sure should use and float(O9_wbaV)/float(num_1) > 0.1
          if 'O-9,46_wzy__1191' in O_dict:  # and float(O946_wzy)/float(num_1) > 0.1
              highest_O = "O-9,46"
          elif "O-9,46,27_partial_wzy__1019" in O_dict:  # and float(O94627)/float(num_1) > 0.1
              highest_O = "O-9,46,27"
          else:
              highest_O = "O-9"  # next, detect O9 vs O2?
              O2 = 0
              O9 = 0
              for z in Special_dict:
                  if "tyr-O-9" in z:
                      O9 = float(Special_dict[z])
                  if "tyr-O-2" in z:
                      O2 = float(Special_dict[z])
              if O2 > O9:
                  highest_O = "O-2"
      elif ("O-3,10_wzx__1539" in O_dict) and (
              "O-9,46_wzy__1191" in O_dict
      ):  # and float(O310_wzx)/float(num_1) > 0.1 and float(O946_wzy)/float(num_1) > 0.1
          if "O-3,10_not_in_1,3,19__1519" in O_dict:  # and float(O310_no_1319)/float(num_1) > 0.1
              highest_O = "O-3,10"
          else:
              highest_O = "O-1,3,19"
      ### end of special test for O9,46 and O3,10 family
      else:
          try:
              max_score = 0
              for x in O_dict:
                  if float(O_dict[x]) >= max_score:
                      max_score = float(O_dict[x])
                      highest_O = x.split("_")[0]
              if highest_O == "O-1,3,19":
                  highest_O = '-'
                  max_score = 0
                  for x in O_dict:
                      if x == 'O-1,3,19_not_in_3,10__130':
                          pass
                      else:
                          if float(O_dict[x]) >= max_score:
                              max_score = float(O_dict[x])
                              highest_O = x.split("_")[0]
          except:
              pass
  #call_fliC:
  if len(H_dict)!=0:
    highest_H_score_both_BC=H_dict[max(H_dict.keys(), key=(lambda k: H_dict[k]))] #used to detect whether fljB existed or not
  else:
    highest_H_score_both_BC=0
  highest_fliC = '-'
  highest_fliC_raw = '-'
  highest_Score = 0
  log_file.write("\nH_scores:\n")
  for s in H_dict:
      log_file.write(s+"\t"+str(H_dict[s])+"\n")
      if s.startswith('fliC'):
          if float(H_dict[s]) > highest_Score:
              highest_fliC = s.split('_')[1]
              highest_fliC_raw = s
              highest_Score = float(H_dict[s])
  #call_fljB
  highest_fljB = '-'
  highest_fljB_raw = '-'
  highest_Score = 0
  for s in H_dict:
      if s.startswith('fljB'):
          if float(H_dict[s]) > highest_Score and float(H_dict[s]) > highest_H_score_both_BC * 0.65: #fljB is special, so use highest_H_score_both_BC to give a general estimate of coverage, currently 0.65 seems pretty good; the reason use a high (0.65) is some fliC and fljB shared with each other
              highest_fljB = s.split('_')[1]
              highest_fljB_raw = s
              highest_Score = float(H_dict[s])
  log_file.write("\nSpecial_scores:\n")
  for s in Special_dict:
    log_file.write(s+"\t"+str(Special_dict[s])+"\n")
  log_file.close()
  return highest_O,highest_fliC,highest_fljB

def get_temp_file_names(for_fq,rev_fq):
  #seqsero2 -a; get temp file names
  sam=for_fq+".sam"
  bam=for_fq+".bam"
  sorted_bam=for_fq+"_sorted.bam"
  mapped_fq1=for_fq+"_mapped.fq"
  mapped_fq2=rev_fq+"_mapped.fq"
  combined_fq=for_fq+"_combined.fq"
  for_sai=for_fq+".sai"
  rev_sai=rev_fq+".sai"
  return sam,bam,sorted_bam,mapped_fq1,mapped_fq2,combined_fq,for_sai,rev_sai

def map_and_sort(threads,database,fnameA,fnameB,sam,bam,for_sai,rev_sai,sorted_bam,mapping_mode):
  #seqsero2 -a; do mapping and sort
  print("building database...")
  subprocess.check_call("bwa index "+database+ " 2>> data_log.txt",shell=True)
  print("mapping...")
  if mapping_mode=="mem":
    subprocess.check_call("bwa mem -k 17 -t "+threads+" "+database+" "+fnameA+" "+fnameB+" > "+sam+ " 2>> data_log.txt",shell=True)
  elif mapping_mode=="sam":
    if fnameB!="":
      subprocess.check_call("bwa aln -t "+threads+" "+database+" "+fnameA+" > "+for_sai+ " 2>> data_log.txt",shell=True)
      subprocess.check_call("bwa aln -t "+threads+" "+database+" "+fnameB+" > "+rev_sai+ " 2>> data_log.txt",shell=True)
      subprocess.check_call("bwa sampe "+database+" "+for_sai+" "+ rev_sai+" "+fnameA+" "+fnameB+" > "+sam+ " 2>> data_log.txt",shell=True)
    else:
      subprocess.check_call("bwa aln -t "+threads+" "+database+" "+fnameA+" > "+for_sai+ " 2>> data_log.txt",shell=True)
      subprocess.check_call("bwa samse "+database+" "+for_sai+" "+for_fq+" > "+sam)
  subprocess.check_call("samtools view -@ "+threads+" -F 4 -Sh "+sam+" > "+bam,shell=True)
  ### check the version of samtools then use differnt commands
  samtools_version=subprocess.Popen(["samtools"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out, err = samtools_version.communicate()
  version = str(err).split("ersion:")[1].strip().split(" ")[0].strip()
  print("check samtools version:",version)
  ### end of samtools version check and its analysis
  if LooseVersion(version)<=LooseVersion("1.2"):
    subprocess.check_call("samtools sort -@ "+threads+" -n "+bam+" "+fnameA+"_sorted",shell=True)
  else:
    subprocess.check_call("samtools sort -@ "+threads+" -n "+bam+" >"+sorted_bam,shell=True)

def extract_mapped_reads_and_do_assembly_and_blast(current_time,sorted_bam,combined_fq,mapped_fq1,mapped_fq2,threads,fnameA,fnameB,database,mapping_mode):
  #seqsero2 -a; extract, assembly and blast
  subprocess.check_call("bamToFastq -i "+sorted_bam+" -fq "+combined_fq,shell=True)
  if fnameB!="":
    subprocess.check_call("bamToFastq -i "+sorted_bam+" -fq "+mapped_fq1+" -fq2 "+mapped_fq2 + " 2>> data_log.txt",shell=True)#2> /dev/null if want no output
  else:
    pass
  outdir=current_time+"_temp"
  print("assembling...")
  if int(threads)>4:
    t="4"
  else:
    t=threads
  if os.path.getsize(combined_fq)>100 and os.path.getsize(mapped_fq1)>100:#if not, then it's "-:-:-"
    if fnameB!="":
      subprocess.check_call("spades.py --careful --pe1-s "+combined_fq+" --pe1-1 "+mapped_fq1+" --pe1-2 "+mapped_fq2+" -t "+t+" -o "+outdir+ " >> data_log.txt 2>&1",shell=True)
    else:
      subprocess.check_call("spades.py --careful --pe1-s "+combined_fq+" -t "+t+" -o "+outdir+ " >> data_log.txt 2>&1",shell=True)
    new_fasta=fnameA+"_"+database+"_"+mapping_mode+".fasta"
    subprocess.check_call("mv "+outdir+"/contigs.fasta "+new_fasta+ " 2> /dev/null",shell=True)
    #os.system("mv "+outdir+"/scaffolds.fasta "+new_fasta+ " 2> /dev/null") contigs.fasta
    subprocess.check_call("rm -rf "+outdir+ " 2> /dev/null",shell=True)
    print("blasting...","\n")
    xmlfile="blasted_output.xml"#fnameA+"-extracted_vs_"+database+"_"+mapping_mode+".xml"
    subprocess.check_call('makeblastdb -in '+new_fasta+' -out '+new_fasta+'_db '+'-dbtype nucl >> data_log.txt 2>&1',shell=True) #temp.txt is to forbid the blast result interrupt the output of our program###1/27/2015
    subprocess.check_call("blastn -query "+database+" -db "+new_fasta+"_db -out "+xmlfile+" -outfmt 5 >> data_log.txt 2>&1",shell=True)###1/27/2015; 08272018, remove "-word_size 10"
  else:
    xmlfile="NA"
  return xmlfile

def judge_subspecies(fnameA):
  #seqsero2 -a; judge subspecies on just forward raw reads fastq
  salmID_output=subprocess.Popen("SalmID.py -i "+fnameA,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out, err = salmID_output.communicate()
  out=out.decode("utf-8")
  file=open("data_log.txt","a")
  file.write(out)
  file.close()
  salm_species_scores=out.split("\n")[1].split("\t")[6:]
  salm_species_results=out.split("\n")[0].split("\t")[6:]
  max_score=0
  max_score_index=1 #default is 1, means "I"
  for i in range(len(salm_species_scores)):
    if max_score<float(salm_species_scores[i]):
      max_score=float(salm_species_scores[i])
      max_score_index=i
  prediction=salm_species_results[max_score_index].split(".")[1].strip().split(" ")[0]
  if float(out.split("\n")[1].split("\t")[4]) > float(out.split("\n")[1].split("\t")[5]): #bongori and enterica compare
    prediction="bongori" #if not, the prediction would always be enterica, since they are located in the later part
  if max_score<10:
    prediction="-"
  return prediction

def judge_subspecies_Kmer(Special_dict):
  #seqsero2 -k;
  max_score=0
  prediction="-" #default should be I
  for x in Special_dict:
    if "mer" in x:
      if max_score<float(Special_dict[x]):
        max_score=float(Special_dict[x])
        prediction=x.split("_")[-1].strip()
      if x.split("_")[-1].strip()=="bongori" and float(Special_dict[x])>95:#if bongori already, then no need to test enterica
        prediction="bongori"
        break
  return prediction

def main():
  #combine SeqSeroK and SeqSero2, also with SalmID
  args = parse_args()
  input_file = args.i
  data_type = args.t
  analysis_mode = args.m
  mapping_mode=args.b
  threads=args.p
  make_dir=args.d
  clean_mode=args.c
  k_size=27 #will change for bug fixing
  database="H_and_O_and_specific_genes.fasta"
  dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
  if len(sys.argv)==1:
    subprocess.check_call(dirpath+"/SeqSero2_package.py -h",shell=True)#change name of python file
  else:
    request_id = time.strftime("%m_%d_%Y_%H_%M_%S", time.localtime())
    request_id += str(random.randint(1, 10000000))
    if make_dir is None:
      make_dir="SeqSero_result_"+request_id
    if os.path.isdir(make_dir):
      pass
    else:
      subprocess.check_call(["mkdir",make_dir])
    #subprocess.check_call("cp "+dirpath+"/"+database+" "+" ".join(input_file)+" "+make_dir,shell=True)
    subprocess.check_call("ln -sr "+dirpath+"/"+database+" "+" ".join(input_file)+" "+make_dir,shell=True)
  ############################begin the real analysis 
    if analysis_mode=="a":
      if data_type in ["1","2","3"]:#use allele mode
        for_fq,rev_fq=get_input_files(make_dir,input_file,data_type,dirpath)
        os.chdir(make_dir)
        ###add a function to tell input files
        fnameA=for_fq.split("/")[-1]
        fnameB=rev_fq.split("/")[-1]
        current_time=time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
        sam,bam,sorted_bam,mapped_fq1,mapped_fq2,combined_fq,for_sai,rev_sai=get_temp_file_names(fnameA,fnameB) #get temp files id
        map_and_sort(threads,database,fnameA,fnameB,sam,bam,for_sai,rev_sai,sorted_bam,mapping_mode) #do mapping and sort
        xmlfile=extract_mapped_reads_and_do_assembly_and_blast(current_time,sorted_bam,combined_fq,mapped_fq1,mapped_fq2,threads,fnameA,fnameB,database,mapping_mode) #extract the mapped reads and do micro assembly and blast
        if xmlfile=="NA":
          O_choice,fliC_choice,fljB_choice,special_gene_list,contamination_O,contamination_H=("-","-","-",[],"","")
        else:
          Final_list=xml_parse_score_comparision_seqsero(xmlfile) #analyze xml and get parsed results
          file=open("data_log.txt","a")
          for x in Final_list:
            file.write("\t".join(str(y) for y in x)+"\n")
          file.close()
          Final_list_passed=[x for x in Final_list if float(x[0].split("_cov_")[1])>=0.9 and (x[1]>=int(x[0].split("__")[1]) or x[1]>=int(x[0].split("___")[1].split("_")[3]) or x[1]>1000)]
          O_choice,fliC_choice,fljB_choice,special_gene_list,contamination_O,contamination_H=predict_O_and_H_types(Final_list,Final_list_passed) #predict O, fliC and fljB
        subspecies=judge_subspecies(fnameA) #predict subspecies
        ###output
        predict_form,predict_sero,star,star_line,claim=seqsero_from_formula_to_serotypes(O_choice,fliC_choice,fljB_choice,special_gene_list,subspecies)
        contamination_report=""
        if contamination_O!="" and contamination_H=="":
          contamination_report="#detected potential contamination of mixed serotypes from O antigen signals.\n"
        elif contamination_O=="" and contamination_H!="":
          contamination_report="#detected potential contamination of mixed serotypes or potential thrid H phase from H antigen signals.\n"
        elif contamination_O!="" and contamination_H!="":
          contamination_report="#detected potential contamination of mixed serotypes from both O and H antigen signals.\n"
        if clean_mode:
          subprocess.check_call("rm -rf ../"+make_dir,shell=True)
          make_dir="none-output-directory due to '-c' flag"
        else:
          new_file=open("Seqsero_result.txt","w")
          if O_choice=="":
            O_choice="-"
          new_file.write("Output_directory:"+make_dir+"\nInput files:\t"+for_fq+" "+rev_fq+"\n"+"O antigen prediction:\t"+O_choice+"\n"+"H1 antigen prediction(fliC):\t"+fliC_choice+"\n"+"H2 antigen prediction(fljB):\t"+fljB_choice+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted subspecies:\t"+subspecies+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+contamination_report+star+star_line+claim+"\n")#+##
          new_file.close()
          print("\n")
          #subprocess.check_call("cat Seqsero_result.txt",shell=True)
          #subprocess.call("rm H_and_O_and_specific_genes.fasta* *.sra *.bam *.sam *.fastq *.gz *.fq temp.txt *.xml "+fnameA+"*_db* 2> /dev/null",shell=True)
          subprocess.call("rm H_and_O_and_specific_genes.fasta* *.sra *.bam *.sam *.fastq *.gz *.fq temp.txt "+fnameA+"*_db* 2> /dev/null",shell=True)
        print("Output_directory:"+make_dir+"\nInput files:\t"+for_fq+" "+rev_fq+"\n"+"O antigen prediction:\t"+O_choice+"\n"+"H1 antigen prediction(fliC):\t"+fliC_choice+"\n"+"H2 antigen prediction(fljB):\t"+fljB_choice+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted subspecies:\t"+subspecies+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+contamination_report+star+star_line+claim+"\n")#+##
      else:
        print("Allele modes only support raw reads datatype, i.e. '-t 1 or 2 or 3'; please use '-m k'")
    elif analysis_mode=="k":
      ex_dir = os.path.dirname(os.path.realpath(__file__))
      #output_mode = args.mode
      for_fq,rev_fq=get_input_files(make_dir,input_file,data_type,dirpath)
      input_file = for_fq #-k will just use forward because not all reads were used
      os.chdir(make_dir)
      f = open(ex_dir + '/antigens.pickle', 'rb')
      lib_dict = pickle.load(f)
      f.close
      input_Ks=get_input_K(input_file,lib_dict,data_type,k_size)
      O_dict,H_dict,Special_dict=get_kmer_dict(lib_dict,input_Ks)
      highest_O,highest_fliC,highest_fljB=call_O_and_H_type(O_dict,H_dict,Special_dict,make_dir)
      subspecies=judge_subspecies_Kmer(Special_dict)
      if subspecies=="IIb" or subspecies=="IIa":
        subspecies="II"
      predict_form,predict_sero,star,star_line,claim = seqsero_from_formula_to_serotypes(
          highest_O.split('-')[1], highest_fliC, highest_fljB, Special_dict,subspecies)
      if clean_mode:
        subprocess.check_call("rm -rf ../"+make_dir,shell=True)
        make_dir="none-output-directory due to '-c' flag"
      else:
        if highest_O.split('-')[-1]=="":
          O_choice="-"
        else:
          O_choice=highest_O.split('-')[-1]
        #print("Output_directory:"+make_dir+"\tInput_file:"+input_file+"\tPredicted subpecies:"+subspecies + '\tPredicted antigenic profile:' + predict_form + '\tPredicted serotype(s):' + predict_sero)
        new_file=open("Seqsero_result.txt","w")
        new_file.write("Output_directory:"+make_dir+"\nInput files:\t"+input_file+"\n"+"O antigen prediction:\t"+O_choice+"\n"+"H1 antigen prediction(fliC):\t"+highest_fliC+"\n"+"H2 antigen prediction(fljB):\t"+highest_fljB+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted subspecies:\t"+subspecies+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+star+star_line+claim+"\n")#+##
        new_file.close()
        subprocess.call("rm *.fasta* *.fastq *.gz *.fq temp.txt *.sra 2> /dev/null",shell=True)
      print("Output_directory:"+make_dir+"\nInput files:\t"+input_file+"\n"+"O antigen prediction:\t"+O_choice+"\n"+"H1 antigen prediction(fliC):\t"+highest_fliC+"\n"+"H2 antigen prediction(fljB):\t"+highest_fljB+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted subspecies:\t"+subspecies+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+star+star_line+claim+"\n")#+##

if __name__ == '__main__':
  main()
