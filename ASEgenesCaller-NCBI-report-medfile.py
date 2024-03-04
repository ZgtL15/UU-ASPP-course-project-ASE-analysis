#!/usr/bin/env python3

'''

#-------------------------------------------------------------------------------
'''
#######################modle####
import sys
import re
import getopt
import os
import gc
def usage():
    print('''Useage: python script.py [option] [parameter]
    -s/--snp_file        input the snp file
    -b/--bam_file        input the bam/sam file
    -g/--gff_file        input the NCBI genes annotation files
    -o/--output          the output results file
    -r/reads_info        output the med file reads info 
    -m/med_file          output the med file reads number
    -h/--help            show possible options''')
#######################default#####
####################### argv ######
opts, args = getopt.getopt(sys.argv[1:], "hs:b:g:o:r:m:",["help","snp_file=","bam_file=","gff_file=","output=","reads_info=","med_file="])
for op, value in opts:
    if op == "-s" or op == "--snp_file":
        snp_file = value
    elif op == "-b" or op == "--bam_file":
        bam_file = value
    elif op == "-g" or op == "--gff_file":
        gff_file = value    
    elif op == "-o" or op == "--output":
        output = value
    elif op == "-r" or op == "--reads_info":
        reads_info = value
    elif op == "-m" or op == "--med_file":
        med_file = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 9:
    usage()
    sys.exit(1)
f1=open(snp_file)
f2=os.popen('samtools view '+bam_file)
f3=open(gff_file)
f4=open(output,'w')
f5=open(reads_info,'w')
f6=open(med_file,'w')

#load snp dictionary########################
'''
1       612     T       C
1       638     A       C
1       681     G       C
1       1596    T       C
'''
#load sam file ########################
'''
HWI-D00524:32:C4D07ACXX:6:1101:2199:1984        83      2       39256853        255     100M    =       39256777        -176    GTTTCTAAGGTCAAGCTGGTCTAGATCCATACATTTCCTCCAAGGCAAAATGTACCCATGACTTTTTTGTG
HWI-D00524:32:C4D07ACXX:6:1101:2199:1984        163     2       39256777        255     100M    =       39256853        176     AGCGGCCCTAGGAGTCACCTCAGCACCTTCTCCCCACCCTCCGACGCCACCTCCCTCTGGGGCTTGCTTGC
'''
#load gtf file ########################
'''
1       ensembl gene    11193   15975   .       +       .       gene_id "ENSECAG00000012421"; gene_version "1"; gene_name "SYCE1"; gene_source "ensembl"; gene_biotype "protein_coding";
1       ensembl transcript      11193   15975   .       +       .       gene_id "ENSECAG00000012421"; gene_version "1"; transcript_id "ENSECAT00000013004"; transcript_version "1"; gene_name "SYCE1";
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
###############define ########################

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
    elif re.search('^(\d+)S(\d+)M$',reads[5]):    ###32S118M###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0]),reads_length):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)M(\d+)S$',reads[5]):  ####136M14S##
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0])):
            sub_reads_info=reads[9][i]
            decide_ref_or_alt_sub(i,sub_reads_info)
    elif re.search('^(\d+)S(\d+)M(\d+)S$',reads[5]):###7S127M16S###
        pos_number=re.findall('\d+',reads[5])
        for i in range(int(pos_number[0]),int(pos_number[0])+int(pos_number[1])):
            sub_reads_info=reads[9][i]
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
readsdic={}
genelength={}
ASEgenes={}
genesitedic={}
CDSsitedic={}
UTRsitedic={}
exonsitedic={}
chrdic={}
sitedicalt={}
sitedicref={}
reads_ref={}
reads_alt={}
gene_ref={}
gene_alt={}
exon_ref={}
exon_alt={}
UTR_ref={}
UTR_alt={}
CDS_ref={}
CDS_alt={}
genedic={}
genedic1={}
exondic={}
exondic1={}
CDSdic={}
CDSdic1={}
UTRdic={}
UTRdic1={}
geneinfo={}
genename={}
exoninfo={}
exonname={}
CDSinfo={}
CDSname={}
UTRinfo={}
UTRname={}
geneinfodic={}
genelength={}
CDSlength={}
exonsetdic={}
CDSsetdic={}
##################################################
for reads in f2:
    reads=reads.split()
    try:
        if int(reads[4]) >= 50:
            if reads[2] in genome_dict:
                reads_length=len(reads[9])
                decide_ref_or_alt()
                get_region()
            else:
                pass
        else:
            pass
    except:
        pass  
#################### put all keys into reads_ref ###################
del ref_dict
gc.collect()
del alt_dict
gc.collect()
for key,value in reads_alt.items():
    if key in reads_ref.keys():
        pass
    else:
        reads_ref[key]=0
############################## genarate distinguish species reads dictionary ################
for key,value in reads_ref.items():
    try:
        readsdic[key]=str(key)+'\t'\
                +str(reads_region_ref[key])+'\t'\
                +str(reads_region_start.get(key,(reads_region_end.get(key,0)-125)))+'\t'\
                +str(reads_region_end.get(key,(reads_region_start.get(key,0)+125)))+'\t'\
                +str(reads_ref.get(key,0))+'\t'\
                +str(reads_alt.get(key,0))
        f5.write(str(key)+'\t'\
                +str(reads_region_ref[key])+'\t'\
                +str(reads_region_start.get(key,(reads_region_end.get(key,0)-125)))+'\t'\
                +str(reads_region_end.get(key,(reads_region_start.get(key,0)+125)))+'\t'\
                +str(reads_ref.get(key,0))+'\t'\
                +str(reads_alt.get(key,0))+'\n')
    except KeyError:
        pass
f5.close()       
###################### calculate distingguish species reads number ###############################        
species_info_reads_number=int(len(readsdic))
del reads_alt
gc.collect()
del reads_ref
gc.collect()
############################## gene select reads count ################

for key,value in readsdic.items():
    reads=value.split()
    if reads[4]=='0':
        sitedicalt[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]]=''
    elif reads[5]=='0':
        sitedicref[reads[0]+'\t'+reads[1]+'\t'+reads[2]+'\t'+reads[3]]=''
    else:
        pass
for gene in f3:
    if gene[0]=='#':
        pass
    else:
        gene=gene.split('\t')
        if gene[2]=='region':
            x=re.findall('Name=\w*',gene[8])
            if x==[]:
                pass
            else:
                i=x[0].split('=')
                chrdic[gene[0]]=i[1]
        if gene[2]=='gene' or gene[2]=='ncRNA' or gene[2]=='region':
            x=re.findall('GeneID:\w*',gene[8])
            if x==[]:
                pass
            else:
                i=x[0].split(':')
                z=i[1]
                cdsc=[]
                exonsetdic[z]=set()
                CDSsetdic[z]=set()
                x=re.findall('gene=\w*',gene[8])
                if x==[]:
                    name='NA'
                else:
                    i=x[0].split('=')
                    name=i[1]
                x=re.findall('gene_biotype=\w*',gene[8])
                if x==[]:
                    biotype='NA'
                else:
                    i=x[0].split('=')
                    biotype=i[1]
                genesitedic[z]=[z,int(gene[3]),int(gene[4])]
                geneinfo[z]=chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]
                geneinfodic[z]=name+'\t'+biotype+'\t'+gene[6]+'\t'+str(gene[3])+'\t'+str(gene[4])
                gene_alt[chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]]=0
                gene_ref[chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]]=0
                a=[]
                b=[]
                a=[gene[3],gene[4]]
                genename[gene[3]+'\t'+gene[4]]=[]
                if chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000)) in genedic:
                    genedic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))].append(a)
                else:
                    genedic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))]=[]
                    genedic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))].append(a)
                if chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000)) in genedic1:
                    genedic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))].append(a)
                else:
                    genedic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))]=[]
                    genedic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))].append(a)
        elif gene[2]=='exon':
            x=re.findall('GeneID:\w*',gene[8])
            if x==[]:
                pass
            else:
                i=x[0].split(':')
                z=i[1]
            exonsitedic[z+'\t'+gene[3]+'\t'+gene[4]]=[chrdic[gene[0]],gene[3],gene[4]]
            exon_alt[chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]]=0
            exon_ref[chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]]=0
            exoninfo[chrdic[gene[0]]+'\t'+gene[3]+'\t'+gene[4]]=z
            d=[]
            c=[gene[3],gene[4]]
            wy=gene[3]+'-'+gene[4]
            exonsetdic[z].add(wy)
            exonname[gene[3]+'\t'+gene[4]]=[]
            if chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000)) in exondic:
                exondic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))].append(c)
            else:
                exondic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))]=[]
                exondic[chrdic[gene[0]]+'\t'+str(round(int(gene[3])/1000000))].append(c)
            if chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000)) in exondic1:
                exondic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))].append(c)
            else:
                exondic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))]=[]
                exondic1[chrdic[gene[0]]+'\t'+str(round(int(gene[4])/1000000))].append(c)
        elif gene[2]=='CDS':
            x=re.findall('GeneID:\w*',gene[8])
            i=x[0].split(':')
            z=i[1]
            hj=gene[3]+'-'+gene[4]
            CDSsetdic[z].add(hj)
            cdsa=int(gene[3])
            cdsb=int(gene[4])
            cdsc.append(cdsa)
            cdsc.append(cdsb)
            cdsmin=str(min(cdsc))
            cdsmax=str(max(cdsc))
            CDSsitedic[z]=[chrdic[gene[0]],cdsmin,cdsmax]
            CDS_alt[chrdic[gene[0]]+'\t'+cdsmin+'\t'+cdsmax]=0
            CDS_ref[chrdic[gene[0]]+'\t'+cdsmin+'\t'+cdsmax]=0
            CDSinfo[z]=chrdic[gene[0]]+'\t'+cdsmin+'\t'+cdsmax
        else:
            pass
for key,value in exonsetdic.items():
    genelength[key]=0
    for x in value:
        x=x.split('-')
        genelength[key]+=(int(x[1])-int(x[0]))
for key,value in CDSsetdic.items():
    CDSlength[key]=0
    for x in value:
        x=x.split('-')
        CDSlength[key]+=(int(x[1])-int(x[0]))
for key,value in exonsitedic.items():
    try:
        x=key.split('\t')
        a=CDSsitedic[x[0]]
        UTRsitedic[x[0]+'\tstart']=[str(a[0]),str(genesitedic[x[0]][1]),str(a[1])]
        UTRsitedic[x[0]+'\tend']=[str(a[0]),str(a[2]),str(genesitedic[x[0]][2])]
        UTRinfo[x[0]+'\tstart']=str(a[0])+'\t'+str(genesitedic[x[0]][1])+'\t'+str(a[1])
        UTRinfo[x[0]+'\tend']=str(a[0])+'\t'+str(a[2])+'\t'+str(genesitedic[x[0]][2])
    except KeyError:
        pass
for key,value in CDSsitedic.items():
    a=[]
    b=[]
    a=[str(value[1]),str(value[2])]
    CDSname[str(value[1])+'\t'+str(value[2])]=[]
    if value[0]+'\t'+str(round(int(value[1])/1000000)) in CDSdic:
        CDSdic[value[0]+'\t'+str(round(int(value[1])/1000000))].append(a)
    else:
        CDSdic[value[0]+'\t'+str(round(int(value[1])/1000000))]=[]
        CDSdic[value[0]+'\t'+str(round(int(value[1])/1000000))].append(a)
    if value[0]+'\t'+str(round(int(value[2])/1000000)) in CDSdic1:
        CDSdic1[value[0]+'\t'+str(round(int(value[2])/1000000))].append(a)
    else:
        CDSdic1[value[0]+'\t'+str(round(int(value[2])/1000000))]=[]
        CDSdic1[value[0]+'\t'+str(round(int(value[2])/1000000))].append(a)
for key,value in UTRsitedic.items():
    a=[]
    b=[]
    a=[value[1],value[2]]
    UTRname[str(value[1])+'\t'+str(value[2])]=[]
    if value[0]+'\t'+str(round(int(value[1])/1000000)) in UTRdic:
        UTRdic[value[0]+'\t'+str(round(int(value[1])/1000000))].append(a)
    else:
        UTRdic[value[0]+'\t'+str(round(int(value[1])/1000000))]=[]
        UTRdic[value[0]+'\t'+str(round(int(value[1])/1000000))].append(a)
    if value[0]+'\t'+str(round(int(value[2])/1000000)) in UTRdic1:
        UTRdic1[value[0]+'\t'+str(round(int(value[2])/1000000))].append(a)
    else:
        UTRdic1[value[0]+'\t'+str(round(int(value[2])/1000000))]=[]
        UTRdic1[value[0]+'\t'+str(round(int(value[2])/1000000))].append(a)
for key,value in sitedicref.items():
    key=key.split()
    try:
        a=genedic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=genedic1[key[1]+'\t'+str(round(int(key[3])/1000000))]        
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicalt.items():
    key=key.split()
    try:
        a=genedic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=genedic1[key[1]+'\t'+str(round(int(key[3])/1000000))]
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in genename[x[0]+'\t'+x[1]]:
                    pass
                else:
                    genename[x[0]+'\t'+x[1]].append(key[0])
                    gene_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=gene_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicref.items():
    key=key.split()
    try:
        a=exondic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=exondic1[key[1]+'\t'+str(round(int(key[3])/1000000))]        
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicalt.items():
    key=key.split()
    try:
        a=exondic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=exondic1[key[1]+'\t'+str(round(int(key[3])/1000000))]
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in exonname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    exonname[x[0]+'\t'+x[1]].append(key[0])
                    exon_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=exon_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
for key,value in sitedicref.items():
    key=key.split()
    try:
        a=CDSdic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in CDSname[str(x[0])+'\t'+str(x[1])]:
                    pass
                else:
                    CDSname[str(x[0])+'\t'+str(x[1])].append(key[0])
                    CDS_ref[key[1]+'\t'+str(x[0])+'\t'+str(x[1])]=CDS_ref.get(key[1]+'\t'+str(x[0])+'\t'+str(x[1]),0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in CDSname[str(x[0])+'\t'+str(x[1])]:
                    pass
                else:
                    CDSname[str(x[0])+'\t'+str(x[1])].append(key[0])
                    CDS_ref[key[1]+'\t'+str(x[0])+'\t'+str(x[1])]=CDS_ref.get(key[1]+'\t'+str(x[0])+'\t'+str(x[1]),0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=CDSdic1[key[1]+'\t'+str(round(int(key[3])/1000000))]        
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in CDSname[str(x[0])+'\t'+str(x[1])]:
                    pass
                else:
                    CDSname[str(x[0])+'\t'+str(x[1])].append(key[0])
                    CDS_ref[key[1]+'\t'+str(x[0])+'\t'+str(x[1])]=CDS_ref.get(key[1]+'\t'+str(x[0])+'\t'+str(x[1]),0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in CDSname[str(x[0])+'\t'+str(x[1])]:
                    pass
                else:
                    CDSname[str(x[0])+'\t'+str(x[1])].append(key[0])
                    CDS_ref[key[1]+'\t'+str(x[0])+'\t'+str(x[1])]=CDS_ref.get(key[1]+'\t'+str(x[0])+'\t'+str(x[1]),0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicalt.items():
    key=key.split()
    try:
        a=CDSdic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in CDSname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    CDSname[str(x[0])+'\t'+str(x[1])].append(key[0])
                    CDS_alt[key[1]+'\t'+str(x[0])+'\t'+str(x[1])]=CDS_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in CDSname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    CDSname[x[0]+'\t'+x[1]].append(key[0])
                    CDS_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=CDS_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=CDSdic1[key[1]+'\t'+str(round(int(key[3])/1000000))]
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in CDSname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    CDSname[x[0]+'\t'+x[1]].append(key[0])
                    CDS_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=CDS_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in CDSname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    CDSname[x[0]+'\t'+x[1]].append(key[0])
                    CDS_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=CDS_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
for key,value in sitedicref.items():
    key=key.split()
    try:
        a=UTRdic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=UTRdic1[key[1]+'\t'+str(round(int(key[3])/1000000))]        
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_ref[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_ref.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass 
    except KeyError:
        pass
for key,value in sitedicalt.items():
    key=key.split()
    try:
        a=UTRdic[key[1]+'\t'+str(round(int(key[2])/1000000))]
        for x in a:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
    try:
        b=UTRdic1[key[1]+'\t'+str(round(int(key[3])/1000000))]
        for x in b:
            if int(key[2])>=int(x[0]) and int(key[2])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            elif int(key[3])>=int(x[0]) and int(key[3])<int(x[1]):
                if key[0] in UTRname[x[0]+'\t'+x[1]]:
                    pass
                else:
                    UTRname[x[0]+'\t'+x[1]].append(key[0])
                    UTR_alt[key[1]+'\t'+x[0]+'\t'+x[1]]=UTR_alt.get(key[1]+'\t'+x[0]+'\t'+x[1],0)+1
            else:
                pass
    except KeyError:
        pass
for key,value in exon_alt.items():
    if key in exon_ref.keys():
        pass
    else:
        exon_ref[key]=0
hj=[]
for key,value in exon_ref.items():
    a=key.split('\t')
    chr=a[0]
    x=exoninfo[key]
    genekey=geneinfo[x]
    UTRskey=UTRinfo.get(x+'\tstart','NA')
    UTRekey=UTRinfo.get(x+'\tend','NA')
    CDSkey=CDSinfo.get(x,'0\t0\t0')
    wy=CDSkey.split('\t')
    w=geneinfodic[x].split('\t')
    genename=w[0]
    biotype=w[1]
    strand=w[2]
    genes=w[3]
    genee=w[4]
    if x not in hj:
        if strand=='+':
            f6.write('1'+'\t'+'gene'+'\t'+chr+'\t'+genes+'\t'+genee+'\t'+str(genelength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(gene_ref.get(genekey,0))+'\t'\
                    +str(gene_alt.get(genekey,0))+'\t'\
                    +str(gene_ref.get(genekey,0)+gene_alt.get(genekey,0))+'\n'\
                    +'2'+'\t'+'CDS'+'\t'+chr+'\t'+wy[1]+'\t'+wy[2]+'\t'+str(CDSlength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(CDS_ref.get(CDSkey,0))+'\t'\
                    +str(CDS_alt.get(CDSkey,0))+'\t'\
                    +str(CDS_ref.get(CDSkey,0)+CDS_alt.get(CDSkey,0))+'\n'\
                    +'3'+'\t'+'3UTR'+'\t'+chr+'\t'+genes+'\t'+wy[1]+'\t'+str(int(wy[1])-int(genes))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0)+UTR_ref.get(UTRskey,0))+'\n'\
                    +'4'+'\t'+'5UTR'+'\t'+chr+'\t'+wy[2]+'\t'+genee+'\t'+str(int(genee)-int(wy[2]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0)+UTR_ref.get(UTRekey,0))+'\n'\
                    +'5'+'\t'+'exon'+'\t'+str(key)+'\t'+str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(exon_ref.get(key,0))+'\t'\
                    +str(exon_alt.get(key,0))+'\t'\
                    +str(exon_ref.get(key,0)+exon_alt.get(key,0))+'\n')
            ASEgenes['1'+'\t'+'gene'+'\t'+chr+'\t'+genes+'\t'+genee]=str(genelength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(gene_ref.get(genekey,0))+'\t'\
                    +str(gene_alt.get(genekey,0))+'\t'\
                    +str(gene_ref.get(genekey,0)+gene_alt.get(genekey,0))
            ASEgenes['2'+'\t'+'CDS'+'\t'+chr+'\t'+wy[1]+'\t'+wy[2]]=str(CDSlength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(CDS_ref.get(CDSkey,0))+'\t'\
                    +str(CDS_alt.get(CDSkey,0))+'\t'\
                    +str(CDS_ref.get(CDSkey,0)+CDS_alt.get(CDSkey,0))
            ASEgenes['3'+'\t'+'3UTR'+'\t'+chr+'\t'+genes+'\t'+wy[1]]=str(int(wy[1])-int(genes))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0)+UTR_ref.get(UTRskey,0))
            ASEgenes['4'+'\t'+'5UTR'+'\t'+chr+'\t'+wy[2]+'\t'+genee]=str(int(genee)-int(wy[2]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0)+UTR_ref.get(UTRekey,0))
            ASEgenes['5'+'\t'+'exon'+'\t'+str(key)]=str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(exon_ref.get(key,0))+'\t'\
                    +str(exon_alt.get(key,0))+'\t'\
                    +str(exon_ref.get(key,0)+exon_alt.get(key,0))
        elif strand=='-':
            f6.write('1'+'\t'+'gene'+'\t'+chr+'\t'+genes+'\t'+genee+'\t'+str(genelength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(gene_ref.get(genekey,0))+'\t'\
                    +str(gene_alt.get(genekey,0))+'\t'\
                    +str(gene_ref.get(genekey,0)+gene_alt.get(genekey,0))+'\n'\
                    +'2'+'\t'+'CDS'+'\t'+chr+'\t'+wy[1]+'\t'+wy[2]+'\t'+str(CDSlength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(CDS_ref.get(CDSkey,0))+'\t'\
                    +str(CDS_alt.get(CDSkey,0))+'\t'\
                    +str(CDS_ref.get(CDSkey,0)+CDS_alt.get(CDSkey,0))+'\n'\
                    +'3'+'\t'+'3UTR'+'\t'+chr+'\t'+wy[2]+'\t'+genee+'\t'+str(int(genee)-int(wy[2]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0)+UTR_ref.get(UTRekey,0))+'\n'\
                    +'4'+'\t'+'5UTR'+'\t'+chr+'\t'+genes+'\t'+wy[1]+'\t'+str(int(wy[1])-int(genes))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0)+UTR_ref.get(UTRskey,0))+'\n'\
                    +'5'+'\t'+'exon'+'\t'+str(key)+'\t'+str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(exon_ref.get(key,0))+'\t'\
                    +str(exon_alt.get(key,0))+'\t'\
                    +str(exon_ref.get(key,0)+exon_alt.get(key,0))+'\n')
            ASEgenes['1'+'\t'+'gene'+'\t'+chr+'\t'+genes+'\t'+genee]=str(genelength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(gene_ref.get(genekey,0))+'\t'\
                    +str(gene_alt.get(genekey,0))+'\t'\
                    +str(gene_ref.get(genekey,0)+gene_alt.get(genekey,0))
            ASEgenes['2'+'\t'+'CDS'+'\t'+chr+'\t'+wy[1]+'\t'+wy[2]]=str(CDSlength.get(x,0))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(CDS_ref.get(CDSkey,0))+'\t'\
                    +str(CDS_alt.get(CDSkey,0))+'\t'\
                    +str(CDS_ref.get(CDSkey,0)+CDS_alt.get(CDSkey,0))
            ASEgenes['4'+'\t'+'5UTR'+'\t'+chr+'\t'+genes+'\t'+wy[1]]=str(int(wy[1])-int(genes))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0))+'\t'+str(UTR_alt.get(UTRskey,0)+UTR_ref.get(UTRskey,0))
            ASEgenes['3'+'\t'+'3UTR'+'\t'+chr+'\t'+wy[2]+'\t'+genee]=str(int(genee)-int(wy[2]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(UTR_ref.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0))+'\t'+str(UTR_alt.get(UTRekey,0)+UTR_ref.get(UTRekey,0))
            ASEgenes['5'+'\t'+'exon'+'\t'+str(key)]=str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(exon_ref.get(key,0))+'\t'\
                    +str(exon_alt.get(key,0))+'\t'\
                    +str(exon_ref.get(key,0)+exon_alt.get(key,0))
        hj.append(x)
    else:
        f6.write('5'+'\t'+'exon'+'\t'+str(key)+'\t'+str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                 +str(exon_ref.get(key,0))+'\t'\
                 +str(exon_alt.get(key,0))+'\t'\
                 +str(exon_ref.get(key,0)+exon_alt.get(key,0))+'\n')
        ASEgenes['5'+'\t'+'exon'+'\t'+str(key)]=str(int(a[2])-int(a[1]))+'\t'+x+';'+genename+';'+biotype+';'+strand+'\t'\
                    +str(exon_ref.get(key,0))+'\t'\
                    +str(exon_alt.get(key,0))+'\t'\
                    +str(exon_ref.get(key,0)+exon_alt.get(key,0))
            
        
f1.close()
f2.close()
f3.close()
f6.close()
del gene_ref
gc.collect()
del gene_alt
gc.collect()
for key,value in ASEgenes.items():
    reads=value.split()
    if int(reads[0])==0:
        ref=0
        alt=0
        total=0
    else:
        ref=int(reads[2])*1000000000/(int(species_info_reads_number)*float(int(reads[0])))
        alt=int(reads[3])*1000000000/(int(species_info_reads_number)*float(int(reads[0])))
        total=int(reads[4])*1000000000/(int(species_info_reads_number)*float(int(reads[0])))
    f4.write(key+'\t'+reads[0]+'\t'+reads[1]+'\t'+str(round(ref,3))+'\t'+str(round(alt,3))+'\t'+str(round(total,3))+'\n')
f4.close()
