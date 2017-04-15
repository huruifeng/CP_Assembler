#!/usr/bin/python
import sys
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

def EXTEND(seq,length):

    min_repeat = 50
    repeat_domain_len = 100
    
    ## cut off the 200bp at the end for the low quality
    seq = seq[:-100]
    
    origin_length = len(seq)
    
    while True:
        subseq = seq[-length:]
        ##index_loc = len(seq)-length
        queryReads = []
        queryIndexs = []

        index = 1
        while index < length - baitLength:
            baitseq = subseq[index:index+baitLength]
            if baitseq not in ReadsBin:
                index += 1
                continue
            for q_read in ReadsBin[baitseq]:
                res_chk = arith.checkErr(subseq[index:],q_read[:-index],ErrorN)
                if res_chk:
                    queryReads.append(q_read)
                    queryIndexs.append(index)
                    #print baitseq,ReadsBin[baitseq],q_read
                    #ReadsBin[baitseq].remove(q_read)
                    if q_read in BoneReads:
                        BoneReads[q_read] = 0
            index += 1
        

        voteCount={"A":0,"T":0,"G":0,"C":0}
        if len(queryReads)==0:
            print "No reads for extending!"
            seq = seq[:-300]
            break
        else:
            for i in range(len(queryReads)):
                seq2process = queryReads[i]
                seq2pIndex = queryIndexs[i]
                vote = seq2process[readLength-seq2pIndex:readLength-seq2pIndex+1]
                if vote == "A":
                    voteCount["A"] += Unique_Reads_Count[seq2process]
                elif vote == "T":
                    voteCount["T"] += Unique_Reads_Count[seq2process]
                elif vote == "G":
                    voteCount["G"] += Unique_Reads_Count[seq2process]
                elif vote == "C":
                    voteCount["C"] += Unique_Reads_Count[seq2process]
                else:
                    pass
            voteCount_sorted = sorted(voteCount.items(), key=lambda d:d[1], reverse = True)
            selectBase = voteCount_sorted[0][0]
            seq += selectBase
            
            current_len = len(seq)
            if len(seq)%1000 == 0:
                print "Current Length: "+str(current_len)
                        
            delta_l = current_len - origin_length
            
            if extend_length == 0:
                pass
            else:
                if delta_l >= extend_length * 1000:
                    print "Entended "+str(delta_l)+" bp"
                    seq = seq[:-300]
                    break
            
            if current_len > 100000 and arith.checkErr(seq[-100:],seq[0:100],ErrorN):
                print "Circle Complete!"
                break
            
            ## Check Repeat
            repeat_str = seq[-min_repeat:]
            repeat_domain = seq[-repeat_domain_len-min_repeat:-min_repeat]
            if repeat_domain.find(repeat_str) !=-1 and len(seq)>100000:
                print "Repeat Domain"
                break
            
            if len(seq)>150000: ##150K
                print "longer than 150K!"
                break

    return seq
#####################################################################################

def printUsage():
    print "Usage:"
    print "    python ContigExtend.py -b bone_file -1 source_file_1 -2 source_file_2 -l readLength -m minOverlaplength -e extendLength\n"
    print "Note:"
    print "The bone_file could be the result file of Assembly programs with long assembled plastid reads,"\
          "or, the initial reads file that is intended to assembly."
    print "\nOptions:"
    for option_i in option_dict:
        print '    '+option_i + "    " + option_dict[option_i]
###############################################################

## init

baitLength = 100 #min overlap length
readLength = 300
extend_length = 5

ErrorN = 1

Bone_File = "top3.fasta"
file_path_1 = "302_R1_select.fastq"
file_path_2 = "302_R2_select.fastq"


option_dict = {"-b":"The path to input file that is intended to extend.",\
               "-1":"For paired-end reads, one of the source reads file.",\
               "-2":"For paired-end reads, another source reads file corresponding to the -1 option.",\
               "-l":"The length of source reads.",\
               "-m":"The min overlap length between bone_reads and source_reads,default:100 bp",\
               "-e":"Set the Extended Length that you want(default:5, means 5Kbp,0 means as long as it can extend).",\
               "-v":"Print the version of this program.",\
               "-h":"Print the Usage options of this program."}

if len(sys.argv) >= 2:
    option_list = sys.argv[1:]
    #print option_list
    running_dict = {}
    if len(option_list)==1:
        if option_list[0] =="-h":
            printUsage()
            sys.exit(1)
        elif option_list[0] =="-v":
            print "Version 17.04.06"
            sys.exit(1)
        else:
            print "ERROR!!!"
            printUsage()
            sys.exit(1)
    for option_i in range(0,len(option_list),2):
        running_dict[option_list[option_i]] = option_list[option_i+1]
    for running_i in running_dict:
        if running_i not in option_dict:
            print "Option "+ running_i +" is invalid. It is not one of the Options."
            sys.exit(1)
        else:
            ##TODO
            
            if running_i =="-b":
                Bone_File = running_dict[running_i]
            if running_i =="-1":
                file_path_1 = running_dict[running_i]
            if running_i =="-2":
                file_path_2 = running_dict[running_i]
            if running_i =="-l":
                readLength = int(running_dict[running_i])
            if running_i =="-e":
                extend_length = int(running_dict[running_i])
            if running_i =="-m":
                baitLength = int(running_dict[running_i])
               
else:
    printUsage()
    #sys.exit(1)

##########################################################

fq_reads_1st.clear()
fq_reads_2nd.clear()
fq_reads_3rd.clear()
fq_reads_4th.clear()
if file_path_1.split(".")[-1].strip() == "fa" or file_path_1.split(".")[-1].strip() == "fasta":
    read_fasta(file_path_1)
    print "Reads Number in File_1:",len(fq_reads_2nd)
    read_fasta(file_path_2)
    print "Reads Number in File_2:", len(fq_reads_2nd)
else:
    read_fastq(file_path_1)
    print "Reads Number in File_1:", len(fq_reads_2nd)
    read_fastq(file_path_2)
    print "Reads Number in File_2:", len(fq_reads_2nd)

print "Reads Processing ..."
Unique_Reads_Count = {}

for id_i in fq_reads_2nd:
    #print id_i
    read_i = fq_reads_2nd[id_i]
    if read_i in Unique_Reads_Count:
        Unique_Reads_Count[read_i] += 1
    else:
        Unique_Reads_Count[read_i] = 1
        
    read_i_revcomp = REVCOMP(read_i)
    if read_i_revcomp in Unique_Reads_Count:
        Unique_Reads_Count[read_i_revcomp] += 1
    else:
        Unique_Reads_Count[read_i_revcomp] = Unique_Reads_Count[read_i]
    

ReadsBin = {}
for read_i in Unique_Reads_Count:
    #print read_i
    baitSeq = read_i[0:baitLength]
    if baitSeq in ReadsBin:
        ReadsBin[baitSeq].append(read_i)
    else:
        ReadsBin[baitSeq] = []
        ReadsBin[baitSeq].append(read_i)

fq_reads_1st.clear()
fq_reads_2nd.clear()
fq_reads_3rd.clear()
fq_reads_4th.clear()

BoneReads = {}

read_fasta(Bone_File)
for id_i in fq_reads_2nd:
    read_i = fq_reads_2nd[id_i]
    if read_i in BoneReads:
        BoneReads[read_i]+=1
    else:
        BoneReads[read_i] = 1

i = 0L

fp_contig_result = open("contig_extend.fasta",'w')
for read_i in BoneReads:
    if BoneReads[read_i] == 0:
        continue
    i += 1
    print "Extending Contig "+str(i)+"..."
    seq3pExtend = EXTEND(read_i,readLength)
    
    ## length is longer than 10kbp and first 100bp == last 100bp, means that it is a circle.
    if len(seq3pExtend)>10000 and arith.checkErr(seq3pExtend[-100:],seq3pExtend[0:100],ErrorN):
        fp_contig_result.write(">Contig_"+str(i)+"|Length:"+str(len(seq3pExtend))+"\n")
        write_fasta(fp_contig_result,seq3pExtend)
        ######################
        fp_contig_result.write(">Contig_"+str(i)+"_RC|Length:"+str(len(seq3pExtend))+"\n")
        write_fasta(fp_contig_result,REVCOMP(seq3pExtend))
        break
    ###############################################################################################
    print "RC Extending..."
    seq3pExtend_revcomp = REVCOMP(seq3pExtend)
    #seq3p5pExtend = seq3pExtend_revcomp
    seq3p5pExtend = EXTEND(seq3pExtend_revcomp,readLength)
    
    print "Current Contig Length: "+ str(len(seq3p5pExtend))
    
    fp_contig_result.write(">Contig_"+str(i)+"|Length:"+str(len(seq3p5pExtend))+"\n")
    write_fasta(fp_contig_result,seq3p5pExtend)
    ###############
    fp_contig_result.write(">Contig_"+str(i)+"_RC|Length:"+str(len(seq3p5pExtend))+"\n")
    write_fasta(fp_contig_result,REVCOMP(seq3p5pExtend))

fp_contig_result.close()
  
print "Contig Extend Running End !!!"

    

        
