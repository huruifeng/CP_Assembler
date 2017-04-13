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
    
#######################################################################################
def contigTrim(contig,miniOverlap = 100):
    
    contig_len = len(contig)
    
    seq_s = contig 
    seq_q = REVCOMP(contig)
    index = miniOverlap
    
    contig_trimed = seq_s
    
    while index < contig_len:
        if index % 100 == 0:
            print "Search..."+str(index/100)+"/"+str(contig_len/100)+" ["+str(i+1)+"/"+str(contig_num)+"]"
        check_res = getOverlapSeq(seq_s[0:index],seq_q[0:index],0,5)
        if not check_res:
            contig_trimed = seq_s[:-index+1]
            break
        else:
            index += 1
    print index
    return contig_trimed
########################################################################################


fq_reads_2nd.clear()
read_fasta("contig_extend.fa")

contigs = []
scaffold_res = {}
for id_i in fq_reads_2nd:
    read_i = fq_reads_2nd[id_i]
    contigs.append(read_i)


contig_num = len(contigs)

fp_result = open("contig_trimed.fa",'w')
for i in range(contig_num):
    print "Triming Contig_"+str(i+1)+" ["+str(i+1)+"/"+str(contig_num)+"]"
    contig_trimed = contigTrim(contigs[i],20)
    fp_result.write(">Contig_"+str(i+1)+"_Length:"+str(len(contig_trimed))+"\n")
    write_fasta(fp_result,contig_trimed)
fp_result.close()
    
print "Contig Trim Running End !!!"

































    

        

