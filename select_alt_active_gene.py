#----------------------------------------------------------------------------------------------------------
import re,os,sys,logging,time,datetime,getopt
def usage():
    print('''
    #-------------------------------------------------------------------------------
    #Author:Guoxiang Xie(zgtl15@foxmail.com)
    #Time:  2022/04/12
    #Version: 1.0
    ###本脚本对得到的进行SNP位点信息根据基因进行筛选
    Useage: python script.py [option] [parameter]
    -i/--input           input  SNPsplit list
    -o/--out             output result file
    -h/--help            show possible options''')


opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["help","input=","output=",])
for op, value in opts:
    if op == "-i" or op == "---input":
        infile1 = value
    elif op == "-o" or op == "--out":
        out = value
    elif op == "-h" or op == "--help":
        usage()
        sys.exit(1)
if len(sys.argv) < 5:
    usage()
    sys.exit(1)


###infile
'''
1       42888308        42888308        T       C       78      47      0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888578        42888578        A       G       23      0       0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888612        42888612        G       T       17      7       0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888668        42888668        G       C       22      9       0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888676        42888676        C       A       18      13      0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888689        42888689        G       A       26      12      0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888713        42888713        T       C       21      18      0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
1       42888744        42888744        G       A       16      15      0       1       42888156        42888878        ENSBTAG00000004124      CPOX    7       CPOX
'''

###根据ensembl gene name对所得到的位点进行字典的写入，key为ensembl_gene_name
def genedic(infile1):
    f1=open(infile1)
    genedic={}
    for xgx in f1:
        xgx_spl = xgx.split()
        key = xgx_spl[11]
        value = xgx
        if key in genedic and eval(xgx_spl[6]) >5 and eval(xgx_spl[6]) > eval(xgx_spl[7]):  ###对alt进行筛选且other<alt###
            genedic[key].append(value)
        elif eval(xgx_spl[6]) >5 :
            genedic[key] = []
            genedic[key].append(value)   
    f1.close()
    return genedic

genedic = genedic(infile1)
f2 = open(out,'w')
for key,value in genedic.items():
    alt = []
    other = []
    exon = []
    for aaa in value:
        aaa = aaa.split()
        alt.append(eval(aaa[6]))
        other.append(eval(aaa[7]))
        exon.append(eval(aaa[13]))
    exon_dedup=[]  
    [exon_dedup.append(i) for i in exon if i not in exon_dedup]    ###对exon的列表进行去重，这样可以避免多个在同一条exon上的SNP位点影响###
    alt_fil=list(filter(lambda x:x>=5,alt))         ###对alt的位点进行筛选，筛选掉小于counts数目5的点###
    if len(exon_dedup) >=2 and len(alt)>=3:
    	f2.write(f'{key}\n')

f2.close()
