#!/usr/bin/env python3.4

#######################modle####
import sys
import re
import os
import getopt
def usage():
    print('''Useage: python script.py [option] [parameter]
    -b/--bam_file        input the bam/sam file
    -i/--reads_info      get reads species info files 
    -r/--ref             output ref sam results file
    -a/--alt             output alt sam results file
    -h/--help            show possible options''')
####################### argv ######
opts, args = getopt.getopt(sys.argv[1:], "hb:i:r:a:",["help","bam_file=","reads_info=","ref=","alt="])
for op, value in opts:
    if op == "-b" or op == "--bam_file":
        bam_file = value
    elif op == "-i" or op == "--reads_info":
        reads_info = value    
    elif op == "-r" or op == "--ref":
        ref = value
    elif op == "-a" or op == "--alt":
        alt = value   
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 9:
    usage()
    sys.exit(1)
f2=open(reads_info)
'''
####load sam f1###
HWI-7001457:246:C79W4ANXX:2:1101:14832:2309     339     27      37368206        1       125M    =       37368047        -284    GCGAGGGCGGTGAGTTTTATGTAAGTAGGTATGGTTATTTCTGGGACGGTTGTGGGGTAGATATTGTTGGAGATGAAGAATCCGGCAAAAATGCTGCCAATT
HWI-70014
CGGCAGTGGCAGCGGGAGCAGCTGCGCCGCCAGCGTGGAGCGGGAGGCGTGTGGGTGTGGAGAGAAGGGAGCGCCGGACCGGAACGGATGTGTACATCTCATGCTTGATGTGAAGAAA
HWI-7001457:246:C79W4ANXX:2:1101:14793:2291     99      MT      13152   1       125M    =       13271   244     TCTTAATTGGCAGCATTTTTGCCGGATTCTTCATCTCCAACAATATCTACCCCACAACCGTCCCAGAAATAACCATACCTACTTACATAAAACTCACCGCCCTCGCAGTAACCATCCT
#### load reads count file f2###
HWI-7001457:246:C79W4ANXX:2:2216:1528:44977     26      35498819        35499050        4       0
HWI-7001457:246:C79W4ANXX:2:2309:2519:89849     5       9257700 9257928 1       0
HWI-7001457:246:C79W4ANXX:2:2114:3383:52075     12      33034824        33035075        1       0
HWI-7001457:246:C79W4ANXX:2:1103:10384:70864    1       92340920        92341159        0       2
HWI-7001457:246:C79W4ANXX:2:1203:8835:98972     2       56472223        56472508        0       1
HWI-7001457:246:C79W4ANXX:2:2304:15605:9808     30      28465488        28467111        2       0
HWI-7001457:246:C79W4ANXX:2:1215:14218:79456    22      34800418        34800777        2       0
'''
refdic={}
altdic={}
reads_ref={}
reads_alt={}
all_row = 0
ref_row = 0
alt_row = 0
for readscount in f2:
    readscount=readscount.split()
    if int(readscount[5])==0:
        refdic[readscount[0]]=''
    elif int(readscount[4])==0:
        altdic[readscount[0]]=''
    else:
        pass
os.system("samtools view -H "+bam_file+">"+ref)
os.system("samtools view -H "+bam_file+">"+alt)
f1=os.popen('samtools view '+bam_file)
for reads in f1:
    all_row += 1
    reads=reads.split()
    try:
        if reads[0] in refdic:
            info='\t'.join(reads)
            reads_ref[reads[0]]=reads_ref.get(reads[0],'')+'\n'+info
        elif reads[0] in altdic:
            info='\t'.join(reads)
            reads_alt[reads[0]]=reads_alt.get(reads[0],'')+'\n'+info
        else:
            pass
    except IndexError:
        pass
f3=open(ref,'a')
f4=open(alt,'a')
for key,value in reads_ref.items():
    f3.write(reads_ref[key].strip('^\n')+'\n')
    ref_row += 1
for key,value in reads_alt.items():
    f4.write(reads_alt[key].strip('^\n')+'\n')
    alt_row += 1
f1.close()
f2.close()
f3.close()
f4.close()

f5 = open('RefAlt_Allbam.txt','a+')
f5.write(bam_file.split('/')[-1] + '\t' + str(all_row) + '\t' + str(ref_row) + '\t' + str(alt_row) + '\t' + str((ref_row +alt_row)/all_row) + '\n')
f5.close()
