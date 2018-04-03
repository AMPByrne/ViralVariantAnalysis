import os,sys,csv

def writecsv(fname,matrix):
    with open(fname,"wb") as fileOut:
        writer=csv.writer(fileOut)
        writer.writerows(matrix)

args=sys.argv

fin=args[1]
thqual=int(args[2])
softPath=args[3]
errorths=args[4].split("_")


fpileup=open(fin,"r")
fOut=os.path.join(fin[:-7]+"_callQth_"+str(thqual)+".csv")
thmapq=fin.split("_")[-1].split(".")[0]
stats = [["seq","pos","cov","afterFilt","A","C","G","T","First","Second","Third","Four"]]
cont=0
for line in fpileup:
    if cont<10: print line
    cont=cont+1
    stat = [0]*10
    line=line.strip().split("\t")
    genome=line[0]
    pos=line[1]
    ref=line[2]
    cov = int(line[3])
    if cov>0:
        bases = list(line[4])
        for i in range(len(bases)):
            if bases[i]=="." or bases[i]==",":
                bases[i]=ref
        qual = map(ord,list(line[5]))
        qual = [x-33 for x in qual]
        basesFin = [bases[i] for i in range(len(bases)) if qual[i]>=thqual]
        acgt = map(basesFin.count,["A","C","G","T"])
        totalIn = sum([acgt[0],acgt[1],acgt[2],acgt[3]]) 
        acgtsort=sorted([acgt[0],acgt[1],acgt[2],acgt[3]])
        acgtsort = acgtsort[::-1]
        if totalIn>0:
            acgtsort=[round(float(x)/totalIn,4) for x in acgtsort]
        stat = [genome,pos,cov,totalIn,acgt[0],acgt[1],acgt[2],acgt[3]]+acgtsort
        stats = stats + [stat]
    else:
        stats = stats + [[genome,pos,cov]+["NA"]*9]
            
writecsv(fOut,stats)

RoutFileName=fOut.split(os.sep)[-1][:-4]

for errorth in errorths:
    cmd = "Rscript "+os.sep+os.path.join(softPath,"codePlot.r") + " " + fOut + " " + str(errorth)+" "+RoutFileName
    print cmd
    os.system(cmd)	

