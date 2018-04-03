# stats.py collets information related to the mapping process from the SAM and VCF files.

import os, csv, time, sys, os.path, collections, math, numpy
from Bio import SeqIO



def mean(data):
    return float(sum(data))/len(data)

def std(data):
    n = len(data)
    c = mean(data)
    ss = sum((x-c)**2 for x in data)	
    return math.sqrt(ss/(n-1))

def sortFlag(x):
	return int('{0:08b}'.format(x)[4])
 
def countReadsFastq(fname):
    num_lines = sum(1 for line in open(fname))
    return num_lines/4
    
def meanLentghReadsFastq(fname,rn):
    fastq_iter=SeqIO.parse(open(fname),"fastq")
    array=[]
    for rec in fastq_iter:
        array.append(len(rec.seq))
    return round(numpy.mean(array),rn),round(numpy.std(array),rn)
 
def countUniqueReadsInFastq(fname):
    outfile=open(fname[:-5]+"unique.fastq","w")
    #done_ids=set()
    done_seqs=[]
    fastq_iter=SeqIO.parse(open(fname),"fastq")
    cont=0
    contT=0
    for rec in fastq_iter:
        contT=contT+1
        #if str(rec.id) not in done_ids:
        if str(rec.seq) not in done_seqs: 
           done_seqs=done_seqs+[str(rec.seq)]
           SeqIO.write([rec],outfile,"fastq")
           cont=cont+1
    outfile.close()
    cmd = "gzip "+ fname[:-5]+"unique.fastq"
    os.system(cmd)           
    return contT,cont

def doSam2(fname,rn):
    print "Processing sam file: "+fname
    fileIn = open(fname, 'r')
    cont1=0
    for line in fileIn:
        if not line[0]=="@":
            cont1=cont1+1
    fileIn.close()
    unmap = 0
    mapqs =[0]*cont1
    readLen =[0]*cont1
    cont = 0
    fileIn = open(fname, 'r')
    for line in fileIn:
        if not line[0]=="@":
            #if cont%100000==0:
            #    print str(cont)+" of "+str(cont1)+" rows done"
            line = line.strip()
            line = line.split("\t")
            if int(line[1]) & 4:
                unmap=unmap+1
            else:
                mapqs[cont]=int(line[4])
                if len(line)>=9:
                    readLen[cont]=len(line[9])
            cont=cont+1
    fileIn.close()
    per_map = round(float(100*(cont-unmap))/cont,rn)
    per_unmap = round(float(100*(unmap))/cont,rn)
    mapqsM = round(mean(mapqs),rn)
    mapqsMSd = round(std(mapqs),rn)
    mapqsMedian=round(numpy.median(mapqs),rn)
    readLenM = round(mean(readLen),rn)
    readLenSd = round(std(readLen),rn)
    cfmappedsam1,cfmappedsam1u = countUniqueReadsInFastq(os.path.join(args[6],args[1]+".Unique.Sorted.Mapped.TrimIn.nodup.sam.fastq"))
   
    return [cont,cont-unmap,per_map,cfmappedsam1u,mapqsMedian,mapqsM,mapqsMSd,readLenM,readLenSd]

def readfqFile(fname,reflen):
    fileIn = open(fname, 'rb')
    dataOut = fileIn.readlines()
    pplus=dataOut.index("+\n")
    dataOut=dataOut[1:pplus]
    dataOut=[x.rstrip() for x in dataOut]
    dataOut="".join(dataOut)
    fileIn.close()
    dataOut=dataOut+"n"*(reflen-len(dataOut))
    print "fq "+fname+" loadded. "+str(len(dataOut))+" bases." 
    return dataOut

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."

def readRef(fname):
    fileIn = open(fname, 'rb')
    dataOut = fileIn.readlines()[1:]
    dataOut=[x.rstrip() for x in dataOut]
    dataOut="".join(dataOut)
    fileIn.close()
    print "Ref "+fname+" loadded. "+str(len(dataOut))+" bases." 
    return dataOut.upper()

def pileupStats(fname,refname,altname,rn):
	print "Procesing vcf file: "+fname
	ref = list(readRef(refname))
	reflen = len(ref)
	alt = list(readfqFile(altname,reflen))
	print "processing "+str(reflen)+" positions"
	header = ["pos","ref","alt_fq","alt_vcf","cov","qual","gcovRF","gcovRR","gcovAF","gcovAR"]
	stats =[header]+[[ind+1,ref[ind],alt[ind],"NA",0,0,0,0,0,0] for ind in range(reflen)]
	indels =0
	fileIn = open(fname, 'r')
	cont=0
	for line in fileIn:
		if line[0]!="#" and "INDEL" not in line:
			#if cont%1000000==0:
			#	print str(cont)+" of "+str(reflen)+" positions done" 
			line=line.split()
			pos=int(line[1])
			ref_vcf=line[3]
			if line[4]==".":
				alt_vcf=line[3]
			else:
				alt_vcf=line[4]
			qual=float(line[5])
			det=line[7].split(";")
			cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
			if "DP4=" in line[7]:
				gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
			stats[pos]=[pos,ref_vcf,alt[pos-1],alt_vcf,cov,qual]+gcov
			cont=cont+1
		else:
			if line[0]!="#" and "INDEL" in line:
				indels=indels+1
	writeCSV(fname.replace('.pileup.vcf','_alignment_stats.csv'),stats)
	stats = stats[1:]
	equal = [x for x in stats if x[1]==x[2].upper()]
	equaln = len(equal)
	per_equaln = round(float(100*equaln)/len(stats),rn)
	unequal = [x for x in stats if not x[1]==x[2].upper()]
	snps = [x for x in unequal if x[2] in ["a","t","c","g","A","G","T","C"]]
	snpsn = len(snps)
	gaps = [x for x in unequal if x[2] in ["n"]]
	gapsn = len(gaps)
	per_gapsn = round(float(100*gapsn)/len(stats),rn)
	other =  [x for x in unequal if x[2] not in ["n","a","t","c","g","A","G","T","C"]]
	othern = len(other)
	
	# meanCov is based in all the sites, even non-covered sites.
	aux = [x[4] for x in stats]# if x[3] not in ["n"] and not x[4]=="NA"] 
	meanCov = round(mean(aux),rn)
	stdCov = round(std(aux),rn)

	aux = [x[5] for x in stats]# if x[3] not in ["n"] and not x[4]=="NA"]
	meanQual = round(mean(aux),rn)
	stdQual = round(std(aux),rn)

	aux = [sum(x[6:]) for x in stats]# if x[3] not in ["n"] and not x[4]=="NA"]
	meangCov = round(mean(aux),rn)
	stdgCov = round(std(aux),rn)
	return [reflen,equaln,per_equaln,meanCov,meanQual]

def pileupDiversity(fname,callqth,errorth):
    fpileup=open(fname,"r")
    #thmapq=fname.split("_")[-1].split(".")[0]
    contT=0
    contV=0
    for line in fpileup:
        line=line.strip().split("\t")
        cov = int(line[3])
        ref=line[2]
        if cov>0 and ref.upper()!="N":
            bases=list(line[4])
            qual = map(ord,list(line[5]))
            qual = [x-33 for x in qual]
            basesFin = [bases[i] for i in range(len(bases)) if qual[i]>=callqth]            
            basesVin=[basesFin[i] for i in range(len(basesFin)) if not basesFin[i] in [".",","] ]
            if len(basesFin)>0:
                if float(len(basesVin))/len(basesFin)>=errorth:
                    contV=contV+len(basesVin)
            contT=contT+len(basesFin)
            
            
    if contT>0:
        return round(contV*1000000/float(contT),rn)
    else:
        return 0

##############################################
##############################################
########## Alignment stats
#arg[1] isolate new name $nf
#arg[2] run name $7
#arg[3] reference genome $2 reference genome
#arg[4] stats file name $4 stats file name
#arg[5] isolate original name $of
#arg[6] isolate directory $isoldir
#arg[7] results path $5
#Rscript $1/myStats.r $nf $3 $2 $4 $of $isoldir $5

rn=3
args=sys.argv


filename=args[1]
refFile=args[2]
pileupFile=args[3]
callqth=int(args[4])

stats_table=[["File_name","R1reads","R2reads","R1R2reads","FirstMappingReads","UniqueReads","TrimmedReads","SecondMappingReads","AverageLength","stdLength","ReferenceLength","EqualBases","Percentage","MeanCoverage","MeanQuality","DiversityIndexError0","DiversityIndexError0.01","DiversityIndexError0.05"]]

fastqR1=filename+"_R1.fastq"
fastqR2=filename+"_R2.fastq"
fastq=filename+".fastq"
firstMapped=filename+".mapped.sorted.fastq"
unique=filename+".mapped.sorted.unique.fastq"
trimmed=filename+".mapped.sorted.unique.trimmed.fastq"
secondMapped=filename+".mapped.sorted.unique.trimmed.mapped.sorted.fastq"
fqFile=filename+".mapped.sorted.unique.trimmed.mapped.sorted.fq"
vcfFile=filename+".mapped.sorted.unique.trimmed.mapped.sorted.vcf"




R1reads=countReadsFastq(fastqR1)
R2reads=countReadsFastq(fastqR2)
R1R2reads=countReadsFastq(fastq)
firstMappedreads=countReadsFastq(firstMapped)
uniquereads=countReadsFastq(unique)
trimmedreads=countReadsFastq(trimmed)
secondMappedreads=countReadsFastq(secondMapped)
meanlen,stdlen=meanLentghReadsFastq(secondMapped,rn)
[reflen,equaln,per_equaln,meanCov,meanQual]=pileupStats(vcfFile,refFile,fqFile,rn)
diverindex0=pileupDiversity(pileupFile,callqth,0)
diverindex01=pileupDiversity(pileupFile,callqth,0.01)
diverindex005=pileupDiversity(pileupFile,callqth,0.05)



results=[filename.split(os.sep)[-1],R1reads,R2reads,R1R2reads,firstMappedreads,uniquereads,trimmedreads,secondMappedreads,meanlen,stdlen,reflen,equaln,per_equaln,meanCov,meanQual,diverindex0,diverindex01,diverindex005]
stats_table.append([str(x) for x in results])

writeCSV(filename+'.summary.csv',stats_table)
