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

def Contig_Assemble(contig_1,contig_2,miniOverlap = 100):
    len_1 = len(contig_1)
    len_2 = len(contig_2)
    
    scaffolds = []
    
    if len_1 > len_2:
        min_len = len_2
        seq_s = contig_1 
        seq_q = contig_2
    else:
        min_len = len_1
        seq_s = contig_2
        seq_q = contig_1
      
    index = miniOverlap
    while index < min_len:
        if index % 100 ==0:
            print "Searching 1/2... "+str(index/100)+"/"+str(min_len/100)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = arith.checkErr(seq_s[-index:],seq_q[0:index])
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
            print "Searching 2/2... "+str(index/100)+"/"+str(min_len/100)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = arith.checkErr(seq_s[-index:],seq_q_rev[0:index])
        if check_res:
            seq_s += seq_q_rev[index-1:]
            scaffolds.append(seq_s)
            break
        else:
            index += 1  
    
    return scaffolds

########################################################################################

fq_reads_2nd.clear()
read_fasta("contig_trimed.fasta")

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

fp_result = open("scaffolds_draft.fasta",'w')
i = 0
for scaffold_i in scaffold_sorted:
    i += 1
    fp_result.write(">Scaffold_"+str(i)+"_Length:"+str(scaffold_i[1])+"\n")
    write_fasta(fp_result,scaffold_i[0])
fp_result.close()

print "Contig Assembly Running End !!!"

   

        
