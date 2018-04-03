# ViralVariantAnalysis
# Javier Nunez
# CSU, Animal and Plant Health Agency, Gov, UK

Provides the location-site where viral sub-populations divert from a comon consensus genome.

To run this software you will need to have installed:

- python 2.7.6
- Perl 5.18.2
- R (any version)
- SMALT (from the Sange Institute)
- samtools (commands index, map, view, faidx, mpileup) and the script vcfutils.pl
- bcftools (view)
- Picard-tools (SamToFastq)
- trimmomatric

The arguments used during the execution of variantAnalyser.py are set in the text file example.args. Edit this file and modify the right hand side for each row according to your settings and preferences.

    ref_sequence=path/to/reference/reference.fasta (only one sequence allowed)
    fastq_R1=path/to/R1_fastq_file/R1_fastq_file.fastq.gz
    fastq_R2=path/to/R2_fastq_file/R2_fastq_file.fastq.gz
    call_quality_threshold=integer between 0 and 40 with the illumina quality calls
    map_quality_threshold=integer between 0 and 60 with the samtools mapping quality
    software_path=path/to/the/folder/of/this/software
    results_path=path/to/where/you/want/the/results/to/show/up
    trimmomatric=path/to/trimmomatric/software/trimomatric.jar
    vcfutils=path/to/vcfutils/script/vcfutils.pl

Once the arguments are set, save the file and run:
  
    python variantAnalyser.py

If everything run smoothly you should get:

   - A folder called "SampleName_Versus_ReferenceName__callQth_"call_quality_threshold"__mapQth_"map_quality_threshold
that contains the results of the variant analysis.

   - A CSV document called SampleName.summary.csv. This file provides a summary of the different stages of the analysis; R1reads or number of short reads in the R1 fastq file, R2reads or naumber of short reads in the R2 fastq file, R1R2 both fastq file combined, FirstMappingreads or the number of reads that maps to reference on a first mapping, UniqueReads or number of unique short reads (repeated reads are filtered out except for one copy), TrimmedRedas or number of reads left after trimming (max length allowed is 35bp), SecondMappingReads or number of reads that maps onto the reference after carring out another mapping stage, AverageLength or average of the reads length, ReferenceLength in bp, EquaBases or number of bases that equal between the reference and the consensus sequences, Percentage of ReferenceLength, MeanCoverage or average of the number of reads that cover each site in the reference sequence, DiversityIndexErrorXX or the number of calls in all and each of the reads that mapped onto the reference that differ from the consensus sequence per million of calls (i.e. diversity of the sample per million sequenced nucleotides).

   - A table called "SampleName_mapQth_"call_quality_threshold"__mapQth_"map_quality_threshold.csv. This file contains a row for each position-site in the reference sequence with columns; the postion, coverage, coverage after filtering the short reads by the mapping threshold, number of A, C, G and Ts, proportion of the most common nucleotide, proportion of the second most common nucleotide, proportion of the third most common nucleotide and the most uncommon nuclotide.

   - Three plots showing the numbers from the previous table. The plots only shows those sites that contain variavility. 
