#!/usr/bin/python
import common
arith = common.arithmetic()

##FASTQ Reads data structure
fq_reads_1st = {} #tag line
fq_reads_2nd = {} #seq line
fq_reads_3rd = {} #description line
fq_reads_4th = {} #quality line

#####################################################################################

RULE={'A':'T','T':'A','C':'G','G':'C','N':'N','a':'T','t':'A','c':'G','g':'C','n':'N'}
def REVCOMP(seq):
    return "".join(map(lambda x:RULE[x],seq))[::-1] 

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
def contigTrim(contig,miniOverlap = 100):
    
    contig_len = len(contig)
    
    seq_s = contig 
    seq_q = REVCOMP(contig)
    
    index = miniOverlap
    
    contig_trimed = seq_s
    
    while index < contig_len/2:
        if index % 100 == 0:
            print "Search..."+str(index/100)+"/"+str(contig_len/100)+" ["+str(i+1)+"/"+str(contig_num)+"]"
        check_res = arith.checkErr(seq_s[0:index],seq_q[0:index],2)
        if not check_res:
            contig_trimed = seq_s[:-index+1]
            break
        else:
            index += 1
    print index
    return contig_trimed
########################################################################################


fq_reads_2nd.clear()
read_fasta("contig_extend.fasta")

scaffold_res = {}

contig_num = len(fq_reads_2nd)

fp_result = open("contig_trimed.fasta",'w')
i = 0
for id_i in fq_reads_2nd:
    i += 1
    print "Triming Contig_"+str(i+1)+" ["+str(i+1)+"/"+str(contig_num)+"]"
    #contig_trimed = contigTrim(fq_reads_2nd[id_i],20)
    contig_trimed = fq_reads_2nd[id_i]
    fp_result.write(id_i.split(":")[0]+":"+str(len(contig_trimed))+"\n")
    write_fasta(fp_result,contig_trimed)
fp_result.close()
    
print "Contig Trim Running End !!!"

































    

        

