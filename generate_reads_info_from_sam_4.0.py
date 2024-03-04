#!/usr/bin/env python3

'''
#-------------------------------------------------------------------------------
#useage: python script.py snp_list bam_file reads_length output_file results_file
#-------------------------------------------------------------------------------
'''

import sys
import re
import os
import getopt
def usage():
    print('''Useage: python script.py [option] [parameter]
    -s/--snp_file        input the snp file
    -b/--bam_file        input the bam/sam file
    -l/--reads_length    set the reads length (default:125)
    -o/--output          the output results file
    -h/--help            show possible options''')
#######################default
reads_length = 125

opts, args = getopt.getopt(sys.argv[1:], "hs:b:l:o:",["help","snp_file=","bam_file=","reads_length=","output="])
for op, value in opts:
    if op == "-s" or op == "--snp_file":
        snp_file = value
    elif op == "-b" or op == "--bam_file":
        bam_file = value
    elif op == "-l" or op == "--reads_length":
        reads_length = int(value)
    elif op == "-o" or op == "--output":
        output = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) != 9:
    usage()
    sys.exit(1)

f1=open(snp_file)
f2 = os.popen('samtools view '+bam_file)
f3=open(output,'w')
#load snp dictionary########################
'''
1       612     T       C
1       638     A       C
1       681     G       C
1       1596    T       C
'''

ref_dict={}
alt_dict={}
genome_dict={}
for snp in f1:
    snp=snp.split()
    index=snp[0]+'-'+snp[1]
    ref_dict[index]=snp[2]
    alt_dict[index]=snp[3]
    genome_dict[snp[0]]=genome_dict.get(snp[0],0)+1
########################################

def get_region():
    '''
    reads_region_ref: the reads in which chromosome
    reads_region_start: the reads's start in chromosome
    reads_region_end: the reads's end in chromosome
    '''
    reads_region_ref[reads[0]]=reads[2]
    if re.search('-',reads[8]):
        reads_region_end[reads[0]]=int(reads[3])+reads_length
    else:
        reads_region_start[reads[0]]=int(reads[3])

def decide_ref_or_alt_sub(i,sub_reads_info):
    snp_index=reads[2]+'-'+str(int(reads[3])+i)
    if snp_index in ref_dict.keys():
        if sub_reads_info==ref_dict[snp_index]:
            reads_ref[reads[0]]=reads_ref.get(reads[0],0)+1
        elif sub_reads_info==alt_dict[snp_index]:
            reads_alt[reads[0]]=reads_alt.get(reads[0],0)+1
        else:
            pass
    else:
        pass
    
def decide_ref_or_alt():
    '''
    in this function, we match 7 different modle.
    for example:
        125M
        10M5N115M
        10M5N5M10N110M
        107M2D18M
        112M1I12M
        6M1I63M2I55M
        31M1D54M1D40M
        
    '''
    if re.search('^(\d+)M$',reads[5]):  ###125M###
        for i in range(reads_length):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)N(\d+)M$',reads[5]): ###10M5N115M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$',reads[5]):  ###10M5N5M10N110M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3])+int(pos_number[4])):
            sub_reads_info=reads[9][i-int(pos_number[1])-int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)D(\d+)M$',reads[5]):     ###107M2D18M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)I(\d+)M$',reads[5]): ###112M1I12M##
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[2])):
            sub_reads_info=reads[9][i+int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M$',reads[5]): ###31M1D54M1D40M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])):
            sub_reads_info=reads[9][i-int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3]),int(pos_number[0])+int(pos_number[1])+int(pos_number[2])+int(pos_number[3])+int(pos_number[4])):
            sub_reads_info=reads[9][i-int(pos_number[1])-int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)I(\d+)M(\d+)I(\d+)M$',reads[5]): ###31M1I54M1I40M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[2])):
            sub_reads_info=reads[9][i+int(pos_number[1])]
            decide_ref_or_alt_sub(i,sub_reads_info)
        for i in range(int(pos_number[0])+int(pos_number[2]),int(pos_number[0])+int(pos_number[2])+int(pos_number[4])):
            sub_reads_info=reads[9][i+int(pos_number[1])+int(pos_number[3])]
            decide_ref_or_alt_sub(i,sub_reads_info)
    else:
        pass                       
################## define dictionary #############
reads_ref={}
reads_alt={}
reads_name={}
reads_region_ref={}
reads_region_start={}
reads_region_end={}
##################################################


for reads in f2.readlines():
    reads=reads.split()
    try:
        if int(reads[4]) >= 0:
            if reads[2] in genome_dict and reads[6]=='=':
                decide_ref_or_alt()
                get_region()
            else:
                pass
        else:
            pass
    except:
        pass  
#f3.write('reads_name\tchr\tstart\tend\tref\talt\n')


#################### put all keys into reads_ref ###################
for key,value in reads_alt.items():
    if key in reads_ref.keys():
        pass
    else:
        reads_ref[key]=0
####################################################################

############################## output what we want ################
for key,value in reads_ref.items():
    try:
        f3.write(str(key)+'\t'\
                +str(reads_region_ref[key])+'\t'\
                +str(reads_region_start.get(key,(reads_region_end.get(key,0)-reads_length)))+'\t'\
                +str(reads_region_end.get(key,(reads_region_start.get(key,0)+reads_length)))+'\t'\
                +str(reads_ref.get(key,0))+'\t'\
                +str(reads_alt.get(key,0))+'\n')
    except KeyError:
        pass

###################### close file ###############################
f1.close()
f2.close()
f3.close()
