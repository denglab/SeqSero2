import os,sys,glob,time,itertools,subprocess
from Initial_Conditions import phase1
from Initial_Conditions import phase2
from Initial_Conditions import phaseO
from Initial_Conditions import sero
from distutils.version import LooseVersion




def xml_parse_score_comparision_seqsero(xmlfile):
  #used to do seqsero xml analysis
  from Bio.Blast import NCBIXML
  handle=open(xmlfile)
  handle=NCBIXML.parse(handle)
  handle=list(handle)
  List=[]
  List_score=[]
  List_ids=[]
  for i in range(len(handle)):
    if len(handle[i].alignments)>0:
      for j in range(len(handle[i].alignments)):
        score=0
        ids=0
        List.append(handle[i].query.strip()+"___"+handle[i].alignments[j].hit_def)
        for z in range(len(handle[i].alignments[j].hsps)):
          if "last" in handle[i].query or "first" in handle[i].query:
            score+=handle[i].alignments[j].hsps[z].bits
            ids+=float(handle[i].alignments[j].hsps[z].identities)/handle[i].query_length
          else:
            if handle[i].alignments[j].hsps[z].align_length>=30:
              #for the long alleles, filter noise parts
              score+=handle[i].alignments[j].hsps[z].bits
              ids+=float(handle[i].alignments[j].hsps[z].identities)/handle[i].query_length
        List_score.append(score)
        List_ids.append(ids)
  temp=zip(List,List_score,List_ids)
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

def judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(node_name,Final_list_passed):
  #used to predict contig is fliC or fljB, if the differnce score value on above head_and_tail is less than 10 (quite small)
  #also used when no head or tail got blasted score for the contig
  role=""
  for z in Final_list_passed:
    if node_name in z[0]:
      role=z[0].split("_")[0]
      break
  return role


def fliC_or_fljB_judge_from_head_tail_sequence(nodes_list,tail_head_list,Final_list_passed):
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
      #compare two heads (37 > 30)
      #four contigs, most perfect assembly, high quality
      """
      for z in a:
        if "fliC_first_37" in z[0]:
          t1=z[1]
        elif "fljB_first_37" in z[0]:
          t2=z[1]
      if t1>=t2:
        role="fliC"
      else:
        role="fljB"
      """
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list_passed)
    elif len(a)==3:
      """
      #compare the number, because hybrid problem
      temp=[]
      for z in a:
        temp.append(z[0].split("_")[0])
      m,n=Uniq(temp)#only two choices in m or n
      if n[0]>n[1]:
        role=m[0]
      else:
        role=m[1]
      """
      ###however, if the one with highest score is the fewer one, compare their accumulation score
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list_passed)
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
        #print "head and tail not belong to same role, now let's guess based on maximum likelihood"
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list_passed)
        """
        max_unit_score=0
        for z in a:
          unit_score=z[-1]/int(z[0].split("__")[1])
          if unit_score>=max_unit_score:
            role=z[0].split("_")[0]
            max_unit_score=unit_score
        """
        ###need to desgin a algorithm to guess most possible situation for nodes_list, See the situations of test evaluation
    elif len(a)==1:
      #that one
      role,diff=judge_fliC_or_fljB_from_head_tail_for_one_contig(a)
      if diff<20:
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list_passed)
      #role=a[0][0].split("_")[0]
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
        role=judge_fliC_or_fljB_from_whole_contig_blast_score_ranking(x,Final_list_passed)
        role_list.append((role,x))
  return role_list

def decide_contig_roles_for_H_antigen(Final_list):
  #used to decide which contig is FliC and which one is fljB
  contigs=[]
  Final_list_passed=[x for x in Final_list if float(x[0].split("_cov_")[1])>=3.5 and (x[1]>=int(x[0].split("__")[1]) or x[1]>=int(x[0].split("___")[1].split("_")[3]))]
  nodes=[]
  for x in Final_list_passed:
    if x[0].startswith("fl") and "last" not in x[0] and "first" not in x[0]:
      nodes.append(x[0].split("___")[1].strip())
  c,d=Uniq(nodes)#c is node_list
  #print c
  tail_head_list=[x for x in Final_list if ("last" in x[0] or "first" in x[0])]
  roles=fliC_or_fljB_judge_from_head_tail_sequence(c,tail_head_list,Final_list_passed)
  return roles

def Combine(b,c):
  fliC_combinations=[]
  fliC_combinations.append(",".join(c))
  temp_combinations=[]
  for i in range(len(b)):
    for x in itertools.combinations(b,i+1):
      temp_combinations.append(",".join(x))
  for x in temp_combinations:
    temp=[]
    for y in c:
      temp.append(y)
    temp.append(x)
    temp=",".join(temp)
    temp=temp.split(",")
    temp.sort()
    temp=",".join(temp)
    fliC_combinations.append(temp)
  return fliC_combinations

def decide_O_type_and_get_special_genes(Final_list):
  #decide O based on Final_list
  O_choice="?"
  O_list=[]
  special_genes=[]
  Final_list_passed=[x for x in Final_list if float(x[0].split("_cov_")[1])>=3.5 and (x[1]>=int(x[0].split("__")[1]) or x[1]>=int(x[0].split("___")[1].split("_")[3]))]
  nodes=[]
  for x in Final_list_passed:
    if x[0].startswith("O-"):
      nodes.append(x[0].split("___")[1].strip())
    elif not x[0].startswith("fl"):
      special_genes.append(x)
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
    if "O-1,3,19_not_in_3,10" not in x[0] and int(x[0].split("__")[1].split("___")[0])+800 <= int(x[0].split("length_")[1].split("_")[0]):#gene length << contig length; for now give 300*2 (for secureity can use 400*2) as flank region
      pointer=x[0].split("___")[1].strip()#store the contig name
      print pointer
    if pointer!=0:#it has potential merge event
      for y in Final_list:
        if pointer in y[0] and y not in final_O and (y[1]>=int(y[0].split("__")[1].split("___")[0])*1.5 or (y[1]>=int(y[0].split("__")[1].split("___")[0])*y[2] and y[1]>=400)):#that's a realtively strict filter now; if passed, it has merge event and add one more to final_O
          potenial_new_gene=y
          print potenial_new_gene
          break
  if potenial_new_gene!="":
    print "two differnt genes in same contig, fix it for O antigen"
    final_O.append(potenial_new_gene)
  ### end of the two genes on same contig test
  if len(final_O)==0:
    #print "$$$No Otype, due to no hit"#may need to be changed
    O_choice="-"
  else:
    O_list=[]
    for x in final_O:
      O_list.append(x[0].split("__")[0])
      if not "O-1,3,19_not_in_3,10__130" in x[0]:#O-1,3,19_not_in_3,10 is too small, which may affect further analysis
        O_nodes_list.append(x[0].split("___")[1])
    ### special test for O9,46 and O3,10 family
    if "O-9,46_wbaV" in O_list:#not sure should use and float(O9_wbaV)/float(num_1) > 0.1
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
          if "tyr-O-9" in z[0]:
            O9=z[1]
          elif "tyr-O-2" in z[0]:
            O2=z[1]
        if O2>O9:
          O_choice="O-2"
        elif O2<O9:
          pass
        else:
          pass
          #print "$$$No suitable one, because can't distinct it's O-9 or O-2, but O-9 has a more possibility."
    elif ("O-3,10_wzx" in O_list) and ("O-9,46_wzy" in O_list):#and float(O310_wzx)/float(num_1) > 0.1 and float(O946_wzy)/float(num_1) > 0.1
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
          if x[1]>=max_score:
            max_score=x[1]
            O_choice=x[0].split("_")[0]
        if O_choice=="O-1,3,19":
          O_choice=final_O[1][0].split("_")[0]
        #print "$$$Most possilble Otype: ",O_choice
      except:
        pass
        #print "$$$No suitable Otype, or failure of mapping (please check the quality of raw reads)"
  #print "O:",O_choice,O_nodes_list
  return O_choice,O_nodes_list,special_genes,final_O

def seqsero_from_formula_to_serotypes(Otype,fliC,fljB,special_gene_list):
  #like test_output_06012017.txt
  #can add more varialbles like sdf-type, sub-species-type in future (we can conclude it into a special-gene-list)
  from Initial_Conditions import phase1
  from Initial_Conditions import phase2
  from Initial_Conditions import phaseO
  from Initial_Conditions import sero
  seronames=[]
  for i in range(len(phase1)):
    fliC_combine=[]
    fljB_combine=[]
    if phaseO[i]==Otype:
      ### for fliC, detect every possible combinations to avoid the effect of "["
      if phase1[i].count("[")==0:
        fliC_combine.append(phase1[i])
      elif phase1[i].count("[")>=1:
        c=[]
        b=[]
        if phase1[i][0]=="[" and phase1[i][-1]=="]" and phase1[i].count("[")==1:
          content=phase1[i].replace("[","").replace("]","")
          fliC_combine.append(content)
          fliC_combine.append("-")
        else:
          for x in phase1[i].split(","):
            if "[" in x:
              b.append(x.replace("[","").replace("]",""))
            else:
              c.append(x)
          fliC_combine=Combine(b,c) #Combine will offer every possible combinations of the formula, like f,[g],t: f,t  f,g,t
      ### end of fliC "[" detect
      ### for fljB, detect every possible combinations to avoid the effect of "["
      if phase2[i].count("[")==0:
        fljB_combine.append(phase2[i])
      elif phase2[i].count("[")>=1:
        d=[]
        e=[]
        if phase2[i][0]=="[" and phase2[i][-1]=="]" and phase2[i].count("[")==1:
          content=phase2[i].replace("[","").replace("]","")
          fljB_combine.append(content)
          fljB_combine.append("-")
        else:
          for x in phase2[i].split(","):
            if "[" in x:
              d.append(x.replace("[","").replace("]",""))
            else:
              e.append(x)
          fljB_combine=Combine(d,e)
      ### end of fljB "[" detect
      new_fliC=fliC.split(",") #because some antigen like r,[i] not follow alphabetical order, so use this one to judge and can avoid missings
      new_fliC.sort()
      new_fliC=",".join(new_fliC)
      new_fljB=fljB.split(",")
      new_fljB.sort()
      new_fljB=",".join(new_fljB)
      if (new_fliC in fliC_combine or fliC in fliC_combine) and (new_fljB in fljB_combine or fljB in fljB_combine):
        seronames.append(sero[i])
  #analyze seronames
  if len(seronames)==0:
    seronames=["N/A (The predicted antigenic profile does not exist in the White-Kauffmann-Le Minor scheme)"]
  star=""
  star_line=""
  if len(seronames)>1:#there are two possible predictions for serotypes
    star="*"
    star_line="The predicted serotypes share the same general formula:\t"+Otype+":"+fliC+":"+fljB+"\n"##
  print "\n"
  predict_form=Otype+":"+fliC+":"+fljB#
  predict_sero=(" or ").join(seronames)
  ###special test for Enteritidis
  if predict_form=="9:g,m:-":
    sdf="-"
    for x in special_gene_list:
      if x[0].startswith("sdf"):
        sdf="+"
    predict_form=predict_form+"\nSdf prediction:"+sdf
    if sdf=="-":
      star="*"
      star_line="Additional characterization is necessary to assign a serotype to this strain.  Commonly circulating strains of serotype Enteritidis are sdf+, although sdf- strains of serotype Enteritidis are known to exist. Serotype Gallinarum is typically sdf- but should be quite rare. Sdf- strains of serotype Enteritidis and serotype Gallinarum can be differentiated by phenotypic profile or genetic criteria.\n"#+##
      predict_sero="See comments below"
  ###end of special test for Enteritidis
  elif predict_form=="4:i:-":
    predict_sero="potential monophasic variant of Typhimurium"
  elif predict_form=="4:r:-":
    predict_sero="potential monophasic variant of Heidelberg"
  elif predict_form=="4:b:-":
    predict_sero="potential monophasic variant of Paratyphi B"
  elif predict_form=="8:e,h:1,2":
    predict_sero="Newport"
    star="*"
    star_line="Serotype Bardo shares the same antigenic profile with Newport, but Bardo is exceedingly rare."
  claim="The serotype(s) is/are the only serotype(s) with the indicated antigenic profile currently recognized in the Kauffmann White Scheme.  New serotypes can emerge and the possibility exists that this antigenic profile may emerge in a different subspecies.  Identification of strains to the subspecies level should accompany serotype determination; the same antigenic profile in different subspecies is considered different serotypes."##
  if "N/A" in predict_sero:
    claim=""
  if "Typhimurium" in predict_sero or predict_form=="4:i:-":
    normal=0
    mutation=0
    for x in special_gene_list:
      if "oafA-O-4_full" in x[0]:
        normal=x[1]
      elif "oafA-O-4_5-" in x[0]:
        mutation=x[1]
    if normal>mutation:
      #print "$$$Typhimurium"
      pass
    elif normal<mutation:
      predict_sero=predict_sero.strip()+"(O5-)"
      star="*"#
      star_line="Detected the deletion of O5-."
      #print "$$$Typhimurium_O5-"
    else:
      #print "$$$Typhimurium, even no 7 bases difference"
      pass
  return predict_form,predict_sero,star,star_line,claim

def main():
  database=sys.argv[1]#used to extract reads
  mapping_mode=sys.argv[2]#mem or sampe
  threads=sys.argv[3]
  for_fq=sys.argv[4]
  rev_fq=sys.argv[5]
  current_time=time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
  sam=for_fq+".sam"
  bam=for_fq+".bam"
  sorted_bam=for_fq+"_sorted.bam"
  mapped_fq1=for_fq+"_mapped.fq"
  mapped_fq2=rev_fq+"_mapped.fq"
  combined_fq=for_fq+"_combined.fq"
  for_sai=for_fq+".sai"
  rev_sai=rev_fq+".sai"
  print "building database..."
  #os.system("bwa index "+database+ " 2> /dev/null")
  os.system("bwa index "+database+ " 2>> data_log.txt ")
  print "mapping..."
  if mapping_mode=="mem":
    os.system("bwa mem -t "+threads+" "+database+" "+for_fq+" "+rev_fq+" > "+sam+ " 2>> data_log.txt")
  elif mapping_mode=="sam":
    os.system("bwa aln -t "+threads+" "+database+" "+for_fq+" > "+for_sai+ " 2>> data_log.txt")
    os.system("bwa aln -t "+threads+" "+database+" "+rev_fq+" > "+rev_sai+ " 2>> data_log.txt")
    os.system("bwa sampe "+database+" "+for_sai+" "+ rev_sai+" "+for_fq+" "+rev_fq+" > "+sam+ " 2>> data_log.txt")
  os.system("samtools view -@ "+threads+" -F 4 -Sbh "+sam+" > "+bam)
  os.system("samtools view -@ "+threads+" -h -o "+sam+" "+bam)
  ### check the version of samtools then use differnt commands
  samtools_version=subprocess.Popen(["samtools"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out, err = samtools_version.communicate()
  version = err.split("ersion:")[1].strip().split(" ")[0].strip()
  print "check samtools version:",version
  if LooseVersion(version)<=LooseVersion("1.2"):
    os.system("samtools sort -@ "+threads+" -n "+bam+" "+for_fq+"_sorted")
  else:
    os.system("samtools sort -@ "+threads+" -n "+bam+" >"+sorted_bam)
  ### end of samtools version check and its analysis
  os.system("bamToFastq -i "+sorted_bam+" -fq "+combined_fq)
  os.system("bamToFastq -i "+sorted_bam+" -fq "+mapped_fq1+" -fq2 "+mapped_fq2 + " 2>> data_log.txt")#2> /dev/null if want no output
  outdir=current_time+"_temp"
  print "assembling..."
  if int(threads)>4:
    t="4"
  else:
    t=threads
  os.system("spades.py --careful --pe1-s "+combined_fq+" --pe1-1 "+mapped_fq1+" --pe1-2 "+mapped_fq2+" -t "+t+" -o "+outdir+ " >> data_log.txt 2>&1")
  new_fasta=for_fq+"_"+database+"_"+mapping_mode+".fasta"
  os.system("mv "+outdir+"/scaffolds.fasta "+new_fasta+ " 2> /dev/null")
  os.system("rm -rf "+outdir+ " 2> /dev/null")
  ### begin blast
  print "blasting..."
  print "\n"
  xmlfile=for_fq+"-extracted_vs_"+database+"_"+mapping_mode+".xml"
  os.system('makeblastdb -in '+new_fasta+' -out '+new_fasta+'_db '+'-dbtype nucl >> data_log.txt 2>&1') #temp.txt is to forbid the blast result interrupt the output of our program###1/27/2015
  os.system("blastn -word_size 10 -query "+database+" -db "+new_fasta+"_db -out "+xmlfile+" -outfmt 5 >> data_log.txt 2>&1")###1/27/2015
  Final_list=xml_parse_score_comparision_seqsero(xmlfile)
  Final_list_passed=[x for x in Final_list if float(x[0].split("_cov_")[1])>=3.5 and (x[1]>=int(x[0].split("__")[1]) or x[1]>=int(x[0].split("___")[1].split("_")[3]))]
  fliC_choice="-"
  fljB_choice="-"
  fliC_contig="NA"
  fljB_contig="NA"
  fliC_length=0 #can be changed to coverage in future
  fljB_length=0 #can be changed to coverage in future
  O_choice=""#no need to decide O contig for now, should be only one
  O_choice,O_nodes,special_gene_list,O_nodes_roles=decide_O_type_and_get_special_genes(Final_list)#decide the O antigen type and also return special-gene-list for further identification
  O_choice=O_choice.split("-")[-1].strip()
  H_contig_roles=decide_contig_roles_for_H_antigen(Final_list)#decide the H antigen contig is fliC or fljB
  log_file=open("SeqSero_hybrid_assembly_log.txt","a")
  print "O_contigs:"
  log_file.write("O_contigs:\n")
  for x in O_nodes_roles:
    if "O-1,3,19_not_in_3,10" not in x[0]:#O-1,3,19_not_in_3,10 is just a small size marker
      print x[0].split("___")[-1],x[0].split("__")[0],"blast score:",x[1],"identity%:",str(round(x[2]*100,2))+"%"
      log_file.write(x[0].split("___")[-1]+" "+x[0].split("__")[0]+" "+"blast score: "+str(x[1])+"identity%:"+str(round(x[2]*100,2))+"%"+"\n")
  print "H_contigs:"
  log_file.write("H_contigs:\n")
  H_contig_stat=[]
  for x in H_contig_roles:
    a=0
    for y in Final_list_passed:
        if x[1] in y[0] and y[0].startswith(x[0]):
            print x[1],x[0],y[0].split("_")[1],"blast_score:",y[1],"identity%:",str(round(y[2]*100,2))+"%"
            log_file.write(x[1]+" "+x[0]+" "+y[0].split("_")[1]+" "+"blast_score: "+str(y[1])+" identity%:"+str(round(y[2]*100,2))+"%"+"\n")
            break
  for x in H_contig_roles:
    #if multiple choices, temporately select the one with longest length for now, will revise in further change
    if "fliC" == x[0] and int(x[1].split("_")[3])>=fliC_length and x[1] not in O_nodes:#remember to avoid the effect of O-type contig, so should not in O_node list
      fliC_contig=x[1]
      fliC_length=int(x[1].split("_")[3])
    elif "fljB" == x[0] and int(x[1].split("_")[3])>=fljB_length and x[1] not in O_nodes:
      fljB_contig=x[1]
      fljB_length=int(x[1].split("_")[3])
  for x in Final_list_passed:
    if fliC_choice=="-" and "fliC_" in x[0] and fliC_contig in x[0] :
      fliC_choice=x[0].split("_")[1]
    elif fljB_choice=="-" and "fljB_" in x[0] and fljB_contig in x[0]:
      fljB_choice=x[0].split("_")[1]
    elif fliC_choice!="-" and fljB_choice!="-":
      break
  print "\n"
  print "SeqSero Input files:",for_fq,rev_fq
  print "Most possible O antigen:",O_choice
  print "Most possible H1 antigen:",fliC_choice
  print "Most possible H2 antigen:",fljB_choice
  #print Final_list
  ###output
  predict_form,predict_sero,star,star_line,claim=seqsero_from_formula_to_serotypes(O_choice,fliC_choice,fljB_choice,special_gene_list)
  new_file=open("Seqsero_result.txt","w")
  new_file.write("Input files:\t"+for_fq+" "+rev_fq+"\n"+"O antigen prediction:\t"+"O-"+O_choice+"\n"+"H1 antigen prediction(fliC):\t"+fliC_choice+"\n"+"H2 antigen prediction(fljB):\t"+fljB_choice+"\n"+"Predicted antigenic profile:\t"+predict_form+"\n"+"Predicted serotype(s):\t"+predict_sero+star+"\n"+star+star_line+claim+"\n")#+##
  new_file.close()
  os.system("cat Seqsero_result.txt")

if __name__ == '__main__':
  main()
