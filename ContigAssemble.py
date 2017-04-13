#!/usr/bin/python

##FASTQ Reads data structure
fq_reads_1st = {} #tag line
fq_reads_2nd = {} #seq line
fq_reads_3rd = {} #description line
fq_reads_4th = {} #quality line

#####################################################################################

RULE={'A':'T','T':'A','C':'G','G':'C','N':'N','a':'T','t':'A','c':'G','g':'C','n':'N'}
def REVCOMP(seq):
    return "".join(map(lambda x:RULE[x],seq))[::-1] 

def read_fastq(fastq_path):
    print "Reading File: "+fastq_path
    read_tag = 0
    name_i = ""
    for line in open(fastq_path, 'r'):
        if read_tag == 0:
            ## name_i: common part in pe1 read and correspongding pe2 read 
            name_i = line.strip().split(" ")[0].split("#")[0]
            fq_reads_1st[name_i] = line.strip()
            read_tag = 1
            continue
        if read_tag == 1:
            fq_reads_2nd[name_i] = line.strip()
            read_tag = 2
            continue
        if read_tag == 2:
            fq_reads_3rd[name_i] = line.strip()
            read_tag = 3
            continue
        if read_tag == 3:
            fq_reads_4th[name_i] = line.strip()
            read_tag = 0
            continue
    print "Reading Complete !"
    return 0
    

def read_fasta(fasta_path):
    print "Reading File: "+fasta_path
    name_i = ""
    sequence = ""
    for line in open(fasta_path, 'r'):
        if line[0]==">":
            if sequence != "": 
                fq_reads_2nd[name_i] = sequence
                sequence = ""
            name_i = line.strip().split(" ")[0].split("#")[0]
        else:
            sequence += line.strip()
    fq_reads_2nd[name_i] = sequence
    print "Reading Complete !"
    return 0
    
######################################################################################

def write_fasta(fp,sequence):
    set_len = 60
    while True:
        if len(sequence) > set_len:
            sub_seq = sequence[0:set_len]
            fp.write(sub_seq+"\n")
            sequence=sequence[set_len:]
        else:
            sub_seq = sequence
            fp.write(sub_seq+"\n")
            break
    return 0 


#######################################################################################
def getOverlapSeq(seq_s,seq_q,index=0,errorN=1):
    seq_s_list = list(seq_s[index:])
    seq_q_list = list(seq_q)
    s_l = len(seq_s_list)-1
    q_l = len(seq_q_list)-1
    error_n = 0
    q_pos = 0
    s_pos = 0
    while s_pos < s_l and q_pos < q_l:
        #print seq_s_list[s_pos],seq_q_list[q_pos]
        if seq_s_list[s_pos] != seq_q_list[q_pos]:
            error_n += 1
            if seq_s_list[s_pos+1] == seq_q_list[q_pos+1]:
                ## mismatch
                s_pos += 1
                q_pos += 1
                continue
            else:
                error_n += 1
            if seq_s_list[s_pos] == seq_q_list[q_pos+1]:
                s_pos += 1
                q_pos += 2
                continue
            else:
                error_n += 1
            if seq_s_list[s_pos+1] == seq_q_list[q_pos]:
                s_pos += 2
                q_pos += 1
                continue
            else:
                error_n += 1

            if error_n > errorN:
                break
        
        s_pos += 1
        q_pos += 1
        
    return False if error_n > errorN else True
#######################################################################################

def Contig_Assemble(contig_1,contig_2,miniOverlap = 100):
    len_1 = len(contig_1)
    len_2 = len(contig_2)
    
    scaffolds = []
    
    if len_1 > len_2:
        max_len = len_1
        min_len = len_2
        seq_s = contig_1 
        seq_q = contig_2
    else:
        max_len = len_2
        min_len = len_1
        seq_s = contig_2
        seq_q = contig_1
      
    index = miniOverlap
    while index < min_len:
        if index % 100 ==0:
            print "Searching... "+str(index/100)+"/"+str(min_len/100)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = checkError(seq_s[-index:],seq_q[0:index],0,50)
        if check_res:
            seq_s += seq_q[index-1:]
            scaffolds.append(seq_s)
            break
        else:
            index += 1        
    
    
    
    if index < min_len:
        return scaffolds
     
    index = miniOverlap
    seq_q_rev = REVCOMP(seq_q)
    
    while index < min_len:
        if index % 100 ==0:
            print "Searching... "+str(index/100)+"/"+str(min_len/100)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = checkError(seq_s[-index:],seq_q_rev[0:index],0,50)
        if check_res:
            seq_s += seq_q_rev[index-1:]
            scaffolds.append(seq_s)
            break
        else:
            index += 1  
    
    return scaffolds
########################################################################################

def checkError(seq_s,seq_q,index=0,errorN=1):
    seq_s_list = list(seq_s[index:])
    seq_q_list = list(seq_q)
    s_l = len(seq_s_list)-1
    q_l = len(seq_q_list)-1
    error_n = 0
    q_pos = 0
    s_pos = 0
    while s_pos < s_l and q_pos < q_l:
        if error_n > errorN:
            break
        if seq_s_list[s_pos] != seq_q_list[q_pos]:
            error_n += 1
            if seq_s_list[s_pos+1] == seq_q_list[q_pos+1]:
                ## mismatch
                s_pos += 1
                q_pos += 1
                continue

            if seq_s_list[s_pos] == seq_q_list[q_pos+1]:
                s_pos += 1
                q_pos += 2
                continue
            
            if seq_s_list[s_pos+1] == seq_q_list[q_pos]:
                s_pos += 2
                q_pos += 1
                continue
            error_n += 1 
        
        s_pos += 1
        q_pos += 1
        
    return 0 if error_n > errorN else q_pos

########################################################################################

fq_reads_2nd.clear()
read_fasta("contig_trimed.fa")

contigs = []
scaffold_res = {}
for id_i in fq_reads_2nd:
    print id_i
    read_i = fq_reads_2nd[id_i]
    contigs.append(read_i)


contig_num = len(contigs)
print "Contig Numbers:",contig_num

n = 0
for i in range(contig_num):
    for j in range(i+1,contig_num):
        n += 1
        if len(contigs[j]) == len(contigs[i]):
            continue
        
        print "Processing Contig_"+str(i+1)+" : Contig_"+str(j+1)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        scaffold_tmp = Contig_Assemble(contigs[i],contigs[j])
        for scaffold_tmp_i in scaffold_tmp:
            if scaffold_tmp_i not in scaffold_res:
                scaffold_res[scaffold_tmp_i]=len(scaffold_tmp_i) 

scaffold_sorted = sorted(scaffold_res.iteritems(), key=lambda d:d[1],reverse = True)

fp_result = open("scaffolds_draft.fa",'w')
i = 0
for scaffold_i in scaffold_sorted:
    i += 1
    fp_result.write(">Scaffold_"+str(i)+"_Length:"+str(scaffold_i[1])+"\n")
    write_fasta(fp_result,scaffold_i[0])
fp_result.close()

print "Contig Assembly Running End !!!"

   

        
