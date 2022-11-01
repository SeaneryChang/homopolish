import sys
import os
import numpy as np
import time
import gzip
import pysam
import pandas as pds
import multiprocessing
import glob
import modules.download as dl
import modules.polish_interface as mlp
import modules.alignment as ma
import modules.getCSV as CSV
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.utils.FileManager import FileManager
from modules.utils.TextColor import TextColor
from modules import ani
from modules.VAtypeClass import FixSNP


#save fasta when modpolish not work
def saveNofixSeq(fixData):
  filePath = fixData.output_dir+"debug/"+fixData.contig_id
  record = SeqRecord(
         fixData.seq,
         id=fixData.contig_id,
         description = fixData.gen_desc
  )
  contig_output_filePath=write_for_new_fasta(record,filePath,fixData.contig_id+"_modpolish")
  return contig_output_filePath

def saveDraftContig(fixData):
  filePath = fixData.output_dir+"debug/"+fixData.contig_id
  record = SeqRecord(
    fixData.seq,
    id=fixData.contig_id,
    description = fixData.gen_desc
  )
  contig_output_filePath=write_for_new_fasta(record,filePath,fixData.contig_id+"_contig")
  return contig_output_filePath

def startModpolsh(fixData,debug_mod):
  modpolish_filePath = ""
  contig_fix_ary = []
  fileName = fixData.draft_genome_file.split('/')[-1].split('.')[0]	
  
 

  #file path end must exist "/"
  if(fixData.output_dir[-1] != "/"):
       fixData.output_dir = fixData.output_dir+"/"


  #create filePath
  if not os.path.isdir(fixData.output_dir+"debug"):
       os.makedirs(fixData.output_dir+"debug")
  
  
  #create bam file
  if(fixData.bamFile == ""):
      fixData.bamFile = getBamPileUp(fixData,fileName,fixData.thread,fixData.draft_genome_file,fixData.reads_file)

  for contig in SeqIO.parse(fixData.draft_genome_file, 'fasta'):
      fixData.contig_id = contig.name#get fa Contig_id
      fixData.seq = contig.seq#get fa seq
      fixData.gen_desc = contig.description
      
      #create each contig debug filePath
      if not os.path.isdir(fixData.output_dir+"debug/"+ fixData.contig_id):
        os.makedirs(fixData.output_dir+"debug/"+ fixData.contig_id)

      #target contig 
      fixData.contig_path = saveDraftContig(fixData)
      ####
      filePath,mod_fix_flag = getPos(fixData,debug_mod)
      contig_fix_ary.append(mod_fix_flag)
      
      if(filePath != ""):
        modpolish_filePath = modpolish_filePath +" "+filePath
  
  #check each contig without modpolish
  if(not any(contig_fix_ary)):
     timestr = time.strftime("[%Y/%m/%d %H:%M]")
     sys.stderr.write(TextColor.PURPLE + str(timestr) + " INFO:Methyl position less than threshold!"+ "\n" + TextColor.END)  
  else:
     fileName = fileName+"_modpolish"
     
  catAllfasta(fixData,modpolish_filePath,fileName)
  
  if(debug_mod):
     shutil.rmtree(fixData.output_dir+"debug")



def getHomoFile(fixData):
   file_path = ""
   flag = True

   #use genomes from user side
   if(fixData.homo_files != ""):
     for homoFile in fixData.homo_files:
        file_path = file_path +" "+homoFile
     
     print(file_path)
     db_path = fixData.output_dir+"debug/"+fixData.contig_id+"/All_homologous_sequences.fna.gz"
     os.system('cat {} > {}'.format(file_path, db_path))
   
   else:
    #download homogenome file from NCBI
    flag = dlHomoFile(fixData)  
   
   return flag


def getPos(fixData,debug_mod):
    

    fileName = fixData.draft_genome_file.split('/')[-1].split('.')[0]
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star modpolish with sequence length: "+ str(len(fixData.seq))  + "\n" + TextColor.END)
    
    #download homologous genomes files
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star download Homologous files"+ "\n" + TextColor.END)
    
    homo_flag = getHomoFile(fixData)
    if(homo_flag == False):
      return saveNofixSeq(fixData),False
    #homologous  array  
    homoPath = fixData.output_dir+"debug/"+fixData.contig_id
    H_misAry,H_AllAry = getPileUpAry(fixData,homoPath,homoPath+"/All_homologous_sequences.fna.gz")   
   
    #Reads array   
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star get Reads File data"+ "\n" + TextColor.END)
    
    #get read mis array
    R_misAry_bam,R_AllAry_bam,totalCovergae = MismatchPileup_read_bam(fixData.bamFile,len(fixData.seq),fixData.seq,fixData.contig_id)

    #check quality and special pattern
    fixary = fix_pos_homo_read_pattern(fixData,H_AllAry,R_misAry_bam,R_AllAry_bam)
    

    #star fix
    modpolish_filePath= fixProcess(fixary,fixData,fileName)
    CSV.getFixPosCSV(fixData.contig_id,fixary,fixData.output_dir)

    return modpolish_filePath,True



  


def catAllfasta(fixData,genFilePath,fileName):
    subFilePath = fixData.draft_genome_file.split('/')[-1].split('.')[0]
    timestr =time.strftime("[%Y/%m/%d %H:%M]")
    # create a directory for each contig
    contig_output_dir = mlp.make_output_dir("contig", fixData.output_dir, subFilePath)
    outPutFilePath = contig_output_dir+"/"+fileName+".fa"
    os.system('cat {} > {}'.format(genFilePath,outPutFilePath)) 
    
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: output dir: " + contig_output_dir + "\n" + TextColor.END)




def getBamPileUp(fixData:FixSNP,fileName,threads,fasta,fastq):
    filePath = fixData.output_dir+'debug/'+fixData.contig_id
    bam = filePath + '/reads_sorted.bam'
    os.system('minimap2 -ax asm5 --cs=long -t {thread} {draft} {reference} | '
              'samtools view -@ {thread} -S -b - | '
              'samtools sort -@ {thread} -o {bam} -'.format(thread=threads, draft=fasta, reference=fastq, bam=bam))
    os.system('samtools index {}'.format(bam))
    bam = filePath+'/reads_sorted.bam'
    return bam
    
    
def getPileUpAry(fixData:FixSNP,pafPath,asemberlyFile):
   filePath = fixData.output_dir+'debug/'+fixData.contig_id
   ma.align(fixData.contig_path,"asm5",fixData.thread,"",filePath,asemberlyFile)
   misAry,posData_ary = MismatchPileup(filePath+"/truth.paf",len(fixData.seq))#get SNP position ATCG array
   
   return misAry,posData_ary   
   
   

def dlHomoFile(fixData:FixSNP):
   fixFlag = True
   filePath = fixData.output_dir+'debug/'+fixData.contig_id
   ncbi_id =  mlp.mash_select_closely_related(fixData.sketch_path,False,fixData.thread,filePath,fixData.mash_threshold,fixData.dl_contig_nums,fixData.contig_path,fixData.contig_id)
   if(len(ncbi_id) == 0):
     return False 
   
   if(len(ncbi_id)<2):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.PURPLE + str(timestr) + "Closely-related genomes less than 2, not to modpolish...\n" + TextColor.END) 
    fixFlag = False
   
   url_list =  dl.parser_url(ncbi_id)
   dl_path = dl.download(filePath,ncbi_id,url_list,fixData.contig_path,99,5)

   return fixFlag



def getHomoVal(S_PosAry,S_percentRate):  # return Homogenome pos value
    Homo_val = ""
    decision_Fix = False

    Homo_A = S_PosAry[0]+S_PosAry[5]
    Homo_T = S_PosAry[1]+S_PosAry[6]
    Homo_C = S_PosAry[2]+S_PosAry[7]
    Homo_G = S_PosAry[3]+S_PosAry[8]
    Homo_Ins_Del = S_PosAry[4]

    S_total = Homo_A+Homo_T+Homo_C+Homo_G
    if(Homo_A >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_T >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_C >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_G >= (S_total*S_percentRate)):
      decision_Fix = True 

    if(decision_Fix):
      Homo_ATCG = [Homo_A,Homo_T,Homo_C,Homo_G,Homo_Ins_Del] #是否連ins del判斷比較好?    
      pos =  Homo_ATCG.index(max(Homo_ATCG))  
      Homo_val = getATCG_pos_use(pos)
      if(pos == 0 and Homo_A == 0):
      	Homo_val = ""

    return Homo_val

def pattern(ST,local_actg):
    patt = True
    if(ST == ""):
        return '',False
    if(ST == 'CCAGC'):
        if(local_actg[0:5] == "CCAAG"):
            value = 'C'
        elif(local_actg[1:6] == "CCGAC" or local_actg[1:6] == "CCAAC" or local_actg[1:6] == "CCAAG"):
            value = 'G'
        elif(local_actg[2:7] == "CCGGC" or local_actg[2:7] == "CCGAC"):
            value = 'A'
        elif(local_actg[2:7] == "GTCGG" or local_actg[2:7] == "GCCGG"):
            value = 'T'
        elif(local_actg[3:8] == "GTCGG" or local_actg[3:8] == "GTTGG" or local_actg[3:8] == "TTTGG"):
            value = 'C'
        elif(local_actg[4:9] == "TTTGG"):
            value = 'G'
        else:
            patt = False
            value = ''    
    elif(ST == 'GCAGC'):
        if(local_actg[1:6] == "GCGAC"):
            value = 'G'
        elif(local_actg[2:7] == "GCGAC" or local_actg[2:7] == "GCGGC"):
            value = 'A'
        elif(local_actg[2:7] == "GTCGC" or local_actg[2:7] == "GCCGC"):
            value = 'T'
        elif(local_actg[3:8] == "GTCGC"):
            value = 'C'
        else:
            patt = False
            value = ''
    else:
        patt = False
        value = ''
    return value,patt

def fix_pos_homo_read_pattern(fixData,S_AllAry,R_MisAry,R_Ary_bam):
    pos_ary = []
    fix_pass = False
    count = 0
    for reads_all in range(0,len(fixData.seq)):
        if(R_MisAry[reads_all][0]!=0):
            S_val=""
            S_percentRate = 1
            diff_D_T = False
            R_precentFlag = True
            misPos = R_MisAry[reads_all]

            R_Pos_ATCG = R_Ary_bam[reads_all]
            R_ATCG_Total = R_Pos_ATCG[0][0]+R_Pos_ATCG[1][0]+R_Pos_ATCG[2][0]+R_Pos_ATCG[3][0]+R_Pos_ATCG[4][0]+R_Pos_ATCG[5][0]+R_Pos_ATCG[6][0]+R_Pos_ATCG[7][0]

            S_val = getHomoVal(S_AllAry[misPos[0]],1)

            if((fixData.seq[misPos[0]] != S_val) and S_val !=""):
               diff_D_T = True
            if((R_Ary_bam[reads_all][0][0]+R_Ary_bam[reads_all][4][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][1][0]+R_Ary_bam[reads_all][5][0])>(R_ATCG_Total*0.95)):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][2][0]+R_Ary_bam[reads_all][6][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][3][0]+R_Ary_bam[reads_all][7][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
 
            if(S_val!="" and diff_D_T  and R_precentFlag):
               ary = [misPos[0],S_val]
               count += 1
               pos_ary.append(ary)

            elif(reads_all>=4 or reads_all <=(len(fixData.seq)-5)):    
                q_all = 0
                n_sum = 0
                local_actg = ""

                for i in range(reads_all-4,reads_all+5):
                    local_actg += fixData.seq[i]
                for j in range(8):
                    n_sum += R_Ary_bam[reads_all][j][0]
                    q_all += R_Ary_bam[reads_all][j][1]
                if(q_all/n_sum < 15):
                    P_val,patt = pattern(fixData.spPattern,local_actg)
                    if(patt):
                        ary = [reads_all,P_val]
                        pos_ary.append(ary)
                    else:
                        H_val = getHomoVal(S_AllAry[reads_all],0.7)
                        ary = [reads_all,H_val]
                        pos_ary.append(ary)
        else:
            q_all = 0
            n_sum = 0
            local_actg = ""
            
            if(reads_all<4 or reads_all>(len(fixData.seq)-5)):
             continue
            for i in range(reads_all-4,reads_all+5):
                local_actg += fixData.seq[i]
            for j in range(8):
                n_sum += R_Ary_bam[reads_all][j][0]
                q_all += R_Ary_bam[reads_all][j][1] 
            if(n_sum == 0):
              continue
            if(q_all/n_sum < 15):
                P_val,patt = pattern(fixData.spPattern,local_actg)
                if(patt):
                    ary = [reads_all,P_val]
                    pos_ary.append(ary)
                else:
                    H_val = getHomoVal(S_AllAry[reads_all],0.7)
                    ary = [reads_all,H_val]
                    pos_ary.append(ary)
    return pos_ary

          
      
def getATCG_pos_use(index:int):
   if(index == 0):
   	return "A"
   elif(index == 1):
   	return "T"
   elif(index == 2):
   	return "C"
   elif(index == 3):
   	return "G"
   else :
   	return ""

def write_for_new_fasta(contig, output_dir_debug,fileName):
    
    contig_name = output_dir_debug + '/' + fileName + '.fasta'
    SeqIO.write(contig, contig_name, "fasta")
    return contig_name

def fixGem(fasta,posAry):
  for ary in  posAry:
     if(ary[1] != ""):
       str1 = str(fasta[:int(ary[0])])
       str2 = str(fasta[int(ary[0])+1:])
       fasta = str1+str(ary[1])+str2

  if(not isinstance(fasta, str)):
      return str(fasta)
  return fasta


def getSeq(fixAry,fixData):
  if(len(fixAry)>1):
    fixSeq = fixGem(fixData.seq,fixAry)   
     #save fasta
    record = SeqRecord(
         Seq(fixSeq),
         id=fixData.contig_id,
         description = fixData.gen_desc
    )
    return record
  else:
    record = SeqRecord(
         fixData.seq,
         id=fixData.contig_id,
         description = fixData.gen_desc
    )
    return record


def fixProcess(fixAry,fixData,fileName):#fix the draft
    
    record = getSeq(fixAry,fixData)
    contig_output_filePath=write_for_new_fasta(record,fixData.output_dir+"debug",fixData.contig_id+"_modpolish")
    return contig_output_filePath
def check_ATCG( char ):
    if char == 'A' or char == 'a':  # A:0, T:1, C:2, G:3
        return 0
    elif char == 'T' or char == 't':
        return 1
    elif char == 'C' or char == 'c':
        return 2
    elif char == 'G' or char == 'g':
        return 3
    elif char == 'N' or char == 'n':
        return 4

def MismatchPileup(file_name, genome_size):

    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star Homo files pileup with sequence length: "+ str(genome_size)  + "\n" + TextColor.END)


    LONG_DELETION_LENGTH = 50
    misAry = np.array([np.array([0,np.zeros(8)])for i in range(genome_size)])
    ins_len = 7
    arr = np.zeros((genome_size, 9), dtype=np.int)
    coverage = np.zeros(genome_size, dtype=np.int)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int)
    over_ins = []

    with open(file_name, 'r') as f:
      for line in f:
        line = line.split()
        t_start = line[7] #reference

        if line[11] != '0': #mapping quality != 0
          cigar = line[-1]
          start_pos = int(t_start)
          flag = 0
          longdel_count = 0
          longdel_status = 0
          if line[4] == "-":#取正反股
            reverse = 5
          else:
            reverse = 0

          i = 5
          cigar_len = len( cigar )
          if( start_pos < genome_size ):
            while i < cigar_len:
              if(start_pos>=genome_size):
                break

              if cigar[i] == "=":
                i += 1
                while i < cigar_len and cigar[i] != "=" and cigar[i] != "*" and cigar[i] != "-" and cigar[i] != "+":
                  base = check_ATCG( cigar[i] )
                  if base != 4: # cigar[i] != 'N' or 'n'
                    arr[start_pos][base + reverse] += 1

                  start_pos += 1
                  i += 1

              elif cigar[i] == "*":

                #base = check_ATCG( cigar[i+1] ) #正解
                #misAry[start_pos][1][base + reverse-1] += 1

                base = check_ATCG( cigar[i+2] ) #原本的
                if base != 4: # cigar[i] != 'N' or 'n'
                  arr[start_pos][base + reverse] += 1

                start_pos +=1
                i += 3

              elif cigar[i] == "+": #start_pos don't need to plus
                i += 1
                while i < cigar_len and cigar[i] != "=" and cigar[i] != "*" and cigar[i] != "-" and cigar[i] != "+":
                  i += 1

              elif cigar[i] == "-":
                i += 1

              else:      #ATCG behind "-"
                i += 1
                start_pos += 1

    return misAry,arr

def MismatchPileup_read_bam(bam, genome_size, ref_seq , contig_name):
    misAry = np.array([np.array([0,np.zeros(8)])for i in range(genome_size)])
    arr = np.zeros((genome_size, 9,3), dtype=np.int)
    totalCovergae = 0
    bamData = pysam.AlignmentFile(bam,'rb')
    ref_name = bamData.get_reference_name(0)
    for each_read in bamData.fetch(until_eof=True):
        x = each_read.mapping_quality
        if(x == 0 or each_read.reference_name != contig_name):
          continue

        #start_pos = each_read.pos
        base_quality = each_read.query_qualities

        qseq = each_read.query_sequence

        if qseq == None:
          continue

        for qpos, refpos in each_read.get_aligned_pairs(True):

          if qseq[qpos] == 'A':
            if ref_seq[refpos] != 'A':
              misAry[refpos][0] = refpos

            arr[refpos][0][0] += 1
            arr[refpos][0][1] += base_quality[qpos]
          elif qseq[qpos] == 'T':
            if ref_seq[refpos] != 'T':
              misAry[refpos][0] = refpos

            arr[refpos][1][0] += 1
            arr[refpos][1][1] += base_quality[qpos]
          elif qseq[qpos] == 'C':
            if ref_seq[refpos] != 'C':
              misAry[refpos][0] = refpos

            arr[refpos][2][0] += 1
            arr[refpos][2][1] += base_quality[qpos]
          elif qseq[qpos] == 'G':
            if ref_seq[refpos] != 'G':
              misAry[refpos][0] = refpos

            arr[refpos][3][0] += 1
            arr[refpos][3][1] += base_quality[qpos]

    #totalCovergae = bamData.count_coverage( ref_name, quality_threshold = 0 )
    return misAry,arr,totalCovergae
