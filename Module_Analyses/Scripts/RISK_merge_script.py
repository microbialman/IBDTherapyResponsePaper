#this script was used to generate a merged version of the individual raw RISK data from GEO (GSE57945_RAW)
#This was used to get the correct gene symbols as the pre-merged rpkm version on GEO has the gene symbols correupted by date conversion in Excel
import sys,gzip

files=sys.argv[2:]
output=sys.argv[1]

rows=[]

gname=False
sampcount=0
for i in files:
    sampcount+=1
    print(sampcount)
    fopen=gzip.open(i,"rU")
    header=False
    rowindex=0
    for j in fopen:
        if gname==False:
            rows.append([])
        if header==False:
            sample=j.strip("\n").split()[4]
            if gname==False:
                rows[rowindex].append("Location\tGene.Symbol\t{}".format(sample))
            else:
                rows[rowindex].append(sample)
            header=True
        else:
            dat=j.strip("\n").strip("\r").split("\t")
            loc=dat[0]
            gene=dat[1]
            count=dat[2]
            if gname==False:
                rows[rowindex].append("{}\t{}\t{}".format(loc,gene,count))
            else:
                rows[rowindex].append(count)
        rowindex+=1
    gname=True

ofile=open(output,"w")
for i in rows:
    ofile.write("{}\n".format("\t".join(i)))
ofile.close()

    
    
