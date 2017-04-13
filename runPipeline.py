import sys
import os

def printUsage():
    print "Usage:"
    print "    python runPipeline.py -b bone_file -1 source_file_1 -2 source_file_2 -l readLength -m minOverlaplength\n"
    print "Note:"
    print "The bone_file could be the result file of Assembly programs with long assembled plastid reads,"\
          "or, the initial reads file that is intended to assembly."
    print "\nOptions:"
    for option_i in option_dict:
        print '    '+option_i + "    " + option_dict[option_i]
###############################################################

## init

baitLength = 100 #min overlap length
readLength = 150

Bone_File = "scaffolds_top5.fa"
file_path_1 = "SR_filtered_1.fa"
file_path_2 = "SR_filtered_2.fa"


option_dict = {"-b":"The path to input file that is intended to extend.",\
               "-1":"For paired-end reads, one of the source reads file.",\
               "-2":"For paired-end reads, another source reads file corresponding to the -1 option.",\
               "-l":"The length of source reads.",\
               "-m":"The min overlap length between bone_reads and source_reads",\
               "-v":"Print the version of this program.",\
               "-h":"Print the Usage options of this program."}

if len(sys.argv) >= 2:
    option_list = sys.argv[1:]
    print option_list
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
                readLength = running_dict[running_i]
            if running_i =="-m":
                baitLength = running_dict[running_i]
               
else:
    printUsage()
    sys.exit(1)

print #####################ContigExtend 1/4###########################
command = "python ContigExtend.py -b "+Bone_File+" -1 "+file_path_1+" -2 "+file_path_2+" -l "+readLength+" -m "+baitLength
os.system(command);

print #####################ContigExtend 2/4###########################
command = "python ContigTrim.py"
os.system(command);

print #####################ContigExtend 3/4###########################
command = "python ContigAssemble.py"
os.system(command);

print #####################ContigExtend 4/4###########################
command = "python CircleCheck.py"
os.system(command);

print ########################-Finidhed-##############################
