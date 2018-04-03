'''
Sep 2017
Javier Nunez
Animal and Plant Agency
The variant analyser software was created to measure and locate the population variability
of viruses. Using a reference sequence (viral genome) the software used SMALT and SAMTOOLS
to map Illumina sequencing data (in format FASTQ) onto the reference. The FASTQ files are
fisrt fitered to eliminate identical short reads, then trimmed and mapped onto the reference.
Unmapped reads (probably from the host or other contaminants) are then removed. The main outputs
are:
- a table ("sample_name.summary.csv") specifying the number of short read left at each stage, and a diversity index 
(number of calls for each and all mapped short reads that disagree with the consensus sequence
(not the reference) per million of calls). The index is provided for three levels of assumed sequenceing
error (0, 0.01, 0.05). 

- A table (sample_name_mapQth_$$_callQth_$$.csv") with one row per site in the reference sequence specifying coverage, number of As, Cs, Gs and Ts
and the percentages of the most common nucleotides.In this table the variant sites can be easily identified

- A graphical visualisation of the main parameters of the previous table using three sequencing error levels 0, 0.01 and 0.05.         

'''

import os,os.path,sys,numpy,csv
from Bio import SeqIO

args=sys.argv
if len(args)>1:
    parametersFile=args[1]
else:
    #just to test
    parametersFile="/media/javier/3TB/Work/WorkVLA/Software/VariantAnalyser/example.args"

print "Reading arguments from file: "+parametersFile

fil=open(parametersFile,"r")
par=""
for line in fil:
    exec(line.strip())
    print "Loading argument ",line.strip()

'''
List of arguments loaded:
ref_sequence............. path to and DNA sequence used as the reference template
fastq_R1..................path to and R1 FASTQ file
fastq_R2..................path to and R2 FASTQ file
call_quality_threshold....minimum quality considered for the sequencing calls (0-40)
map_quality_threshold.....minimum quality considered for the SAMTOOLS mapping short reads (0-60) 
software_path.............path to the directory where the variant analyser lives
results_path..............path to the directory where the results will be stored
trimmomatic...............path to and trimomatric software
vcfutils..................path to and vcfutils script
'''

def mapping_se_return_mapped(ref,fastq_in):
    '''
    This function creates:
    sample.sam
    sample.bam
    sample.mapped.bam
    sample.mapped.sorted.bam
    sample.mapped.sorted.sam
    sample.mapped.sorted.fastq
    '''
    sample=fastq_in.replace('.fastq','')
    os.system('smalt index -k 13 -s 6 '+ref+' '+ref)
    os.system('smalt map  -d -1 -y 0.95 -x -f samsoft -o '+sample+'.sam'+' '+ref+' '+fastq_in)
    os.system('samtools view -Shu '+sample+'.sam > '+sample+'.bam')
    os.system('samtools view -b -F 4 '+sample+'.bam > '+sample+'.mapped.bam')
    os.system('samtools sort '+sample+'.mapped.bam'+' '+sample+'.mapped.sorted')
    os.system('picard-tools SamToFastq I='+sample+'.mapped.sorted.bam F='+sample+'.mapped.sorted.fastq')
    os.system('samtools view -h -o '+sample+'.mapped.sorted.sam'+' '+sample+'.mapped.sorted.bam')

def bam_to_vcf(ref,bam_in):
    sample=bam_in.replace('.bam','')
    os.system('samtools index '+bam_in)
    os.system('samtools faidx '+ref)
    os.system('samtools mpileup -uf '+ref+' '+bam_in+' > '+sample+'.bcf')
    os.system('bcftools view -cg '+sample+'.bcf'+' > '+sample+'.vcf')
    os.system('perl '+vcfutils+' vcf2fq '+sample+'.vcf'+' > '+sample+'.fq')
    fq_to_fasta(sample+'.fq')

    
def isACGT(c):
    return c in ["a","A","c","C","g","G","t","T"]
    
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
    letters = list(s)
    new=[]
    for l in letters:
        if isACGT(l):
            new=new+[basecomplement[l]]
        else:
            new=new+[l]
    return ''.join(new)

def revcom(s):
    return complement(s[::-1])

def repeats_filter_se(fastq,fastq_out_unique_file,fastq_out_repeated_file):
    fastq_iter=SeqIO.parse(open(fastq),"fastq")
    contT=0
    contU=0
    contR=0
    readsDict=dict()
    fastq_out_repeated=open(fastq_out_repeated_file,"wb")
    for record in fastq_iter:
        contT=contT+1
        rec=record.format("fastq").split("\n")[:-1]  
        sequence=rec[1]
        qualities=record.letter_annotations["phred_quality"]
        if sequence not in readsDict and revcom(sequence) not in readsDict:
            readsDict[sequence]=[rec[0],rec[3],numpy.min(qualities),numpy.mean(qualities)]
            contU=contU+1
        else:
            contR=contR+1
            fastq_out_repeated.write(rec[0]+"\n")
            fastq_out_repeated.write(rec[1]+"\n")
            fastq_out_repeated.write(rec[2]+"\n")
            fastq_out_repeated.write(rec[3]+"\n")
            if rec[1] in readsDict:
                sequence=rec[1]
            else:
                sequence=revcom(rec[1])
            minQualInDict=readsDict[sequence][-2]
            minQualrec=numpy.min(qualities)
            meanQualInDict=readsDict[sequence][-1]
            meanQualrec=numpy.mean(qualities)
            if meanQualrec>meanQualInDict and minQualrec>minQualInDict:
                readsDict[sequence]=[rec[0],rec[3],minQualrec,meanQualrec]
            else:
                readsDict[sequence][-1]=readsDict[sequence][-1]+1
    fastq_out_repeated.close()
    fastq_out_unique=open(fastq_out_unique_file,"wb")
    for key in readsDict:
        fastq_out_unique.write(readsDict[key][0]+"\n")
        fastq_out_unique.write(key+"\n")
        fastq_out_unique.write("+"+"\n")
        fastq_out_unique.write(readsDict[key][1]+"\n")
    fastq_out_unique.close()
    print "Total number of reads: "+str(contT)
    print "Unique reads: "+str(contU)
    print "Repeated reads: "+str(contR)
    return 0

def read_fq(fname):
    fileIn = open(fname, 'rb')
    lines= fileIn.readlines()
    pluses = [i for i in range(0,len(lines)) if lines[i]=="+\n"]
    ini=0
    ids=[]
    seqs=[]
    for plus in pluses:
        ids.append(lines[ini][1:-1])
        seqs.append("".join([lines[i][:-1] for i in range(ini+1,plus)]))
        ini = plus+(plus-ini)
    return ids,seqs
    
def fq_to_fasta(fq_in,ncols=100):
    ids,seqs=read_fq(fq_in)
    fasta_out_file=fq_in.replace('.fq','.fasta')
    fasta_out=open(fasta_out_file,"w")
    for idi,seq in zip(ids,seqs):
        fasta_out.write(">"+idi+"\n")
        for i in range(0,len(seq),ncols):
            fasta_out.write(seq[i:i+ncols]+"\n")
    fasta_out.close()

def spn_filter(snp_in,thqual=150,thprop=0.2,thmincov=10):
    csv_out_file=snp_in.replace('.spn','.filtered.csv')
    csv_out=open(csv_out_file,"w")
    writer = csv.writer(csv_out, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
    with open(snp_in) as infile:
        for line in infile:
            if line[0]!="#" and "INDEL" not in line:
                line=line.split()
                refGenome=line[0]
                pos=int(line[1])
                ref=line[3]
                alt=line[4]
                qual=float(line[5])
                det=line[7].split(";")
                cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
                gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
                if len(alt)==1 and qual>=thqual and gcov[0]<=thprop*gcov[2] and gcov[1]<=thprop*gcov[3] and (gcov[2]>thmin or gcov[3]>thmin):
                    writer.writerow([refGenome,pos,ref,alt,qual,cov]+gcov)
    csv_out.close()

def variant_calling(bcf_in):
    '''
    This function creates a csv file containing the SNPs filtered by several QC measures.
    '''
    sample=bcf_in.replace('.bcf','')
    os.system('bcftools view -vcg '+bcf_in+' > '+sample+'.snp')
    #os.system('perl '+vcfutils+' varFilter '++sample+'.snp'+' > '+sample+'.flt.snp')
    spn_filter(sample+'.snp')



ref_name=ref_sequence.split(os.sep)[-1].split(".")[0]
sample_name=fastq_R1.split(os.sep)[-1].split("_")[0]
results_path=os.path.join(results_path,sample_name+"__Versus__"+ref_name+"__callQth_"+str(call_quality_threshold)+"__mapQth_"+str(map_quality_threshold))
stats_File=os.path.join(results_path,sample_name+"__Versus__"+ref_name+"__callQth_"+str(call_quality_threshold)+"__mapQth_"+str(map_quality_threshold)+".summary.txt")

print "Creating results directory "+results_path
os.system('mkdir -p '+results_path)

print "Unzipping and copying FASTQ files and reference sequence"
R1=os.path.join(results_path,sample_name+'_R1.fastq')
R2=os.path.join(results_path,sample_name+'_R2.fastq')
R1R2=os.path.join(results_path,sample_name+'.fastq')
sample=R1R2.replace('.fastq','')
ref=os.path.join(results_path,ref_sequence.split(os.sep)[-1])

os.system('gunzip -c '+fastq_R1+' > '+R1)
os.system('gunzip -c '+fastq_R2+' > '+R2)
os.system('cat '+R1+' '+R2+' > '+R1R2)
os.system('cp '+ref_sequence+' '+ref)

print "Fist mapping"
mapping_se_return_mapped(ref,R1R2)

print "Filtering exact repeats"
repeats_filter_se(R1R2,R1R2.replace('.fastq','.mapped.sorted.unique.fastq'),R1R2.replace('.fastq','.mapped.sorted.repeated.fastq'))

print "Trimming mapped reads" 
os.system('java -jar '+trimmomatric+' SE -phred33 -trimlog '+R1R2.replace('.fastq','.trim.log')+' '+R1R2.replace('.fastq','.mapped.sorted.unique.fastq')+' '+R1R2.replace('.fastq','.mapped.sorted.unique.trimmed.fastq')+' SLIDINGWINDOW:10:20 MINLEN:36')

print "Second mapping"
mapping_se_return_mapped(ref,R1R2.replace('.fastq','.mapped.sorted.unique.trimmed.fastq'))

print "Generating the vcf file"
bam_to_vcf(ref, R1R2.replace('.fastq','.mapped.sorted.unique.trimmed.mapped.sorted.bam'))

print "Variant calling"
variant_calling(R1R2.replace('.fastq','.mapped.sorted.unique.trimmed.mapped.sorted.bcf'))

print "Generating the old pileup version using the CIGAR format infomation"
os.system('python '+software_path+'/OldPileupGenerator.py'+' '+sample+'.mapped.sorted.unique.trimmed.mapped.sorted.sam'+' '+ref+' '+str(map_quality_threshold))

print "Generating diversity plots"
os.system('python '+software_path+'/pileupCounter.py '+sample+'_mapQth_'+str(map_quality_threshold)+'.pileup'+' '+str(call_quality_threshold)+' '+software_path+' 0_0.01_0.05 '+sample+'_mapQth_'+str(map_quality_threshold))

print "Collecting info to fill up stats table"
os.system('python '+software_path+'/stats_SE.py '+sample+' '+ref+' '+sample+'_mapQth_'+str(map_quality_threshold)+'.pileup'+' '+str(map_quality_threshold))

print "Deleting unwanted files"
os.system('rm '+sample+'*.sam')
os.system('rm '+sample+'*.fastq')
os.system('rm '+sample+'*.log')
os.system('rm '+sample+'*.pileup')
os.system('rm '+sample+'*.vcf')
os.system('rm '+sample+'*.bcf')
os.system('rm '+sample+'*.fq')
os.system('rm '+sample+'*.snp')
os.system('rm '+sample+'*.trimmed.mapped.bam')
os.system('rm '+sample+'*.mapped.bam')
os.system('rm '+sample+'*.unique.trimmed.bam')
os.system('rm '+sample+'.mapped.sorted.bam')
os.system('rm '+sample+'.bam')
os.system('rm '+ref+'.fai')
os.system('rm '+ref)
os.system('rm '+ref+'*.smi')
os.system('rm '+ref+'*.sma')

