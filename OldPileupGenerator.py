import os, sys

args=sys.argv

if len(args)<2:
    fname ="/media/javier/4TB1/MoreWork/Virology/Sara/NGS_sub420_Dec17/VariantAnalysis/RVF01-PoolsAB__Versus__RVF01-RVFV_ZH501_singleREF_11979bp_consensus_callQth_30__mapQth_60/RVF01-PoolsAB.filtered.mapped.sorted.sam"
    reffile = "/media/javier/4TB1/MoreWork/Virology/Sara/NGS_sub420_Dec17/VariantAnalysis/RVF01-PoolsAB__Versus__RVF01-RVFV_ZH501_singleREF_11979bp_consensus_callQth_30__mapQth_60/RVF01-RVFV_ZH501_singleREF_11979bp_consensus.fa"
    mapqth = 0
else:    
    fname =sys.argv[1]
    reffile=sys.argv[2]
    mapqth=int(sys.argv[3])  


revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

def readRef(fname):
    fileIn = open(fname, 'rb')
    dataOut = fileIn.readlines()
    dataOut=[x.rstrip() for x in dataOut]
    dataOut="".join(dataOut[1:])
    fileIn.close()
    return dataOut
    
def writeListToTextFile(lista,fileName):
    f = open(fileName, 'w')
    f.write("\n".join(lista))
    f.close()
    print "file "+fileName+" saved."
def isNumber(ch):
    if ch in "0123456789":
        return 1
    else:
        return 0
def indChar(strin):
     stin = [isNumber(x) for x in strin]
     return [i for i, x in enumerate(stin) if x == 0]
def cigar(strin):
    inds = indChar(strin)    
    return [strin[:inds[0]]]+[strin[inds[i]+1:inds[i+1]] for i in range(len(inds)-1)],[strin[inds[i]] for i in range(len(inds)-1)]+[strin[-1]]
    

refn=reffile.split(os.sep)[-1].split(".")[0]
ref=readRef(reffile).upper()

pileup = [[refn,str(i+1),ref[i],str(0),"","","",""] for i in range(len(ref))]

fileIn = open(fname, 'r')
cont1=0
for line in fileIn:
    if not line[0]=="@":
        cont1=cont1+1
fileIn.close()

cont=0
with open(fname, "r") as sam:
    for line in sam:
        cont=cont+1
        #if cont%1000000==0:
        #    print str(cont)+" of "+str(cont1)+" rows done" 
        if not line[0]=="@":
            line = line.strip()
            line=line.split("\t")
            flag=int(line[1])
            mapq=int(line[4])
            mapqChr=chr(mapq+33)
            pos=int(line[3])
            if pos>0 and mapq>=mapqth: #line[3]!="0"
                seq = line[9].upper()
                qual =line[10]
                cig = line[5]
                cigaro = cigar(cig)
                char="."
                if flag & 16:
                    #seq=revcompl(seq)
                    char=","
                #for i in range(len(cigaro[0])):
                #    print cigaro
                seqPos=0
                refPos=0
                for i in range(len(cigaro[0])):
                    if cigaro[1][i]=="M":
                        lenMatch=int(cigaro[0][i])
                        for p in range(lenMatch):
                            if pileup[pos-1+refPos+p][2]==seq[seqPos+p]:  
                                pileup[pos-1+refPos+p][4]=pileup[pos-1+refPos+p][4]+char
                            else:
                                pileup[pos-1+refPos+p][4]=pileup[pos-1+refPos+p][4]+seq[seqPos+p]
                            pileup[pos-1+refPos+p][5]=pileup[pos-1+refPos+p][5]+qual[seqPos+p]
                            pileup[pos-1+refPos+p][6]=pileup[pos-1+refPos+p][6]+mapqChr
                    if cigaro[1][i]=="I":
                        seqPos=seqPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["S","H"]:
                        seqPos=seqPos+int(cigaro[0][i])
                        #refPos=refPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["M"]:
                        seqPos=seqPos+int(cigaro[0][i])
                        refPos=refPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["D","N","P"]:
                        refPos=refPos+int(cigaro[0][i])
                
sam.close()

for i in range(len(pileup)): pileup[i][3]=str(len(pileup[i][4]))
pileup=["\t".join(x) for x in pileup]
newfname=fname.split(".")[0]+"_mapQth_"+str(mapqth)+".pileup"
writeListToTextFile(pileup,newfname)

            