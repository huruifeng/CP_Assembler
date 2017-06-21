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
    
    scaffolds = {}

    if len_1 > len_2:
        min_len = len_2
    else:
        min_len = len_1
    
    
    ############################################################################          
    seq_s = contig_1 
    seq_q = contig_2
    
    index = miniOverlap
    while index < min_len:
        if index % 1000 ==0:
            print "Searching 1/2... "+str(index/1000)+"/"+str(min_len/1000)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = arith.checkError_hrf(seq_s[-index:],seq_q[0:index],2)
        if check_res:
            #seq_s += seq_q[index-1:]
            scaffolds[seq_s]=[]
            scaffolds[seq_s].append(seq_q)
            scaffolds[seq_s].append(index)
            break
        else:
            index += 1         
    
    if index < min_len:
        ## assembly is successful!
        return scaffolds
    
    ##############################################################################
    seq_s = contig_2 
    seq_q = contig_1
    
    index = miniOverlap
    while index < min_len:
        if index % 1000 ==0:
            print "Searching 2/2... "+str(index/1000)+"/"+str(min_len/1000)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        check_res = arith.checkError_hrf(seq_s[-index:],seq_q[0:index],2)
        if check_res:
            #seq_s += seq_q[index-1:]
            scaffolds[seq_s]=[]
            scaffolds[seq_s].append(seq_q)
            scaffolds[seq_s].append(index)
            break
        else:
            index += 1         
            
    return scaffolds

########################################################################################

fq_reads_2nd.clear()
read_fasta("contig_trimed.fasta")

read_id = {}
contigs = []
for id_i in fq_reads_2nd:
    print id_i
    read_i = fq_reads_2nd[id_i]
    read_id[read_i] = id_i
    contigs.append(read_i)


contig_num = len(contigs)
print "Contig Numbers:",contig_num

n = 0
Ass_Path = {}
scaffold_res = ""
for i in range(contig_num):
    for j in range(i+1,contig_num):
        n += 1
        if len(contigs[j]) == len(contigs[i]):
            continue
        
        print "Processing Contig_"+str(i+1)+" : Contig_"+str(j+1)+" ["+str(n)+"/"+str(contig_num*(contig_num-1)/2)+"]"
        scaffold_tmp = Contig_Assemble(contigs[i],contigs[j])
        
        for contig_i in scaffold_tmp:
            Ass_Path[contig_i]=scaffold_tmp[contig_i]

##print Ass_Path

Path_X = 0
Ass_Map_Name={}
Ass_Map_Name[Path_X]=[]

Ass_Map={}
Ass_Map[Path_X]=[]
for contig_i in Ass_Path:
    flag = 0
    next_seq = Ass_Path[contig_i][0]
    for Path_i in Ass_Map:
        if contig_i in Ass_Map[Path_i]:
            flag = 1
            pos_i = Ass_Map[Path_i].index(contig_i)
            Ass_Map[Path_i].insert(pos_i+1,next_seq)
            Ass_Map_Name[Path_i].insert(pos_i+1,read_id[next_seq])
        elif next_seq in Ass_Map[Path_i]:
            flag = 2
            pos_i = Ass_Map[Path_i].index(next_seq)
            Ass_Map[Path_i].insert(pos_i,contig_i)
            Ass_Map_Name[Path_i].insert(pos_i,read_id[contig_i])
    if flag == 0:
        Path_X += 1
        Ass_Map[Path_X]=[]
        Ass_Map[Path_X].append(contig_i)
        Ass_Map[Path_X].append(next_seq)
        
        Ass_Map_Name[Path_X]=[]
        Ass_Map_Name[Path_X].append(read_id[contig_i])
        Ass_Map_Name[Path_X].append(read_id[next_seq])

        
del Ass_Map[0]
del Ass_Map_Name[0]

#print Ass_Map
print Ass_Map_Name

scaffolds =[]

for path_i in Ass_Map:
    index = 0
    i = 1
    s_temp = ""
    for contig_i in Ass_Map[path_i]:
        if i == len(Ass_Map[path_i]):
            break
        s_temp += contig_i[index:]
        index = Ass_Path[contig_i][1]
        i += 1
    s_temp += contig_i[index:]
    scaffolds.append(s_temp)

fp_result = open("scaffolds_draft.fasta",'w')
i = 0
for scaffold_i in scaffolds:
    i += 1
    fp_result.write(">Scaffold_"+str(i)+"|Length:"+str(len(scaffold_i))+"\n")
    write_fasta(fp_result,scaffold_i)
fp_result.close()

print "Contig Assembly Running End !!!"

   
      
