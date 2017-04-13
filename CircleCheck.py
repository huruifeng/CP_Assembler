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
######################################################################################
def CircleCheck(scaffold):
        
    cp_genome = scaffold
    cp_len = len(scaffold)

    cp_head = cp_genome
    cp_tail = cp_genome
    
    min_overlap = 100
    
    index = 0
    while index < cp_len/2:
        if index % 1000 == 0:
            print "Checking..."+str(index)+"/"+str(cp_len/2)+"["+str(i)+"/"+str(scaffold_n)+"]"
        if getOverlapSeq(cp_head[0:index+min_overlap],cp_tail[-index-min_overlap:],0,50):
            break
        index+=1
        
    if index >= cp_len / 2:
        return 0
    else:
        return index

#######################################################################################
##Check Circle or Region
fq_reads_2nd.clear()
read_fasta("scaffolds_draft.fa")

scaffolds = {}
for id_i in fq_reads_2nd:
    read_i = fq_reads_2nd[id_i]
    scaffolds[read_i] = len(read_i)
    
scaffold_sorted = sorted(scaffolds.iteritems(), key=lambda d:d[1],reverse = True)

fp_final = open("scaffolds_final.fa",'w')
i = 0
n = 0
scaffold_n = len(scaffolds)
for scaffold_i in scaffold_sorted:
    i+=1
    Circle_check = CircleCheck(scaffold_i[0])
    if Circle_check:
        n+=1
        ##scaffold_temp = scaffold_i[0][Circle_check:]
        scaffold_temp = scaffold_i[0]
        fp_final.write(">Scaffold_"+str(i)+"|Length:"+str(len(scaffold_temp))+"|OverLap:"+str(Circle_check)+"\n")
        write_fasta(fp_final,scaffold_temp)
fp_final.close()

print "CircleCheck Running End !!![Num = "+str(n)+"]"

        















