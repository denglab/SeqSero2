#!/usr/bin/env python

############################################################################
# Copyright (c) 2014-2015 University of Georgia
# Developer: Shaokang Zhang, zskzsk@uga.edu
# All Rights Reserved
############################################################################

import argparse,os,sys,time,random

def main():
  parser = argparse.ArgumentParser(usage='SeqSero2.py -i <input_data> [-p <number of threads>] [-b <BWA_algorithm>]\n\nDevelopper: Shaokang Zhang (zskzsk@uga.edu) and Xiangyu Deng (xdeng@uga.edu)\n\nContact email:seqsero@gmail.com')#add "-m <data_type>" in future
  parser.add_argument("-i",nargs="+", help="<string>: path/to/input_data")
  parser.add_argument("-b",choices=['sam','mem'],default="mem",help="<string>: 'sam'(bwa samse/sampe), 'mem'(bwa mem),default=mem") 
  parser.add_argument("-p",default="1",help="<int>: if p>4, only 4 threads will be used for assembly, default=1") 
  args=parser.parse_args()
  dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
  if len(sys.argv)==1:
    os.system(dirpath+"/SeqSero2.py -h")
  else:
    request_id = time.strftime("%m_%d_%Y_%H_%M_%S", time.localtime())
    request_id += str(random.randint(1, 10000000))
    make_dir="SeqSero_result_"+request_id
    os.system("mkdir "+make_dir)
    os.system("cp -rf "+dirpath+"/database/H_and_O_and_specific_genes.fasta "+make_dir)
    #mode_choice=args.m
    mapping_mode=args.b
    threads=args.p
    dataset=args.i
    os.system("cp "+dataset[0]+" "+make_dir)
    os.system("cp "+dataset[1]+" "+make_dir)
    fnameA=dataset[0].split("/")[-1]
    fnameB=dataset[1].split("/")[-1]
    os.chdir(make_dir)
    os.system("python2.7 "+dirpath+"/libs/mapping_and_assembly_hybrid.py H_and_O_and_specific_genes.fasta "+mapping_mode+" "+str(threads)+" "+fnameA+" "+fnameB)
    os.system("rm H_and_O_and_specific_genes.fasta* *.bam *.sam *.fastq *.fastq.gz *.fq temp.txt *.xml "+fnameA+"*_db* 2> /dev/null")
    print "Output_directory:",make_dir
    #print "\n\n\nResult:\n"
    #os.system("cat Seqsero_result.txt")

if __name__ == '__main__':
  main()
