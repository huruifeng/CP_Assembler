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


#####################################################################################
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
        if arith.checkErr(cp_head[0:index+min_overlap],cp_tail[-index-min_overlap:]):
            break
        index+=1
        
    if index >= cp_len / 2:
        return 0
    else:
        return index

#######################################################################################
##Check Circle or Region
fq_reads_2nd.clear()
read_fasta("scaffolds_draft.fasta")

scaffolds = {}
for id_i in fq_reads_2nd:
    read_i = fq_reads_2nd[id_i]
    scaffolds[read_i] = len(read_i)
    
scaffold_sorted = sorted(scaffolds.iteritems(), key=lambda d:d[1],reverse = True)

fp_final = open("scaffolds_final.fasta",'w')
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

        















