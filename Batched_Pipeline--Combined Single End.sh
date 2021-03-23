#!/bin/bash

#Make the list of files to analyze, based on whats in the folder containing this script
SUFFIX="_R1_001.fastq.gz"
ls | grep $SUFFIX | grep -v "Undetermined" | rev | cut -c17- | rev > filelist.txt


#Build the reference library
#bowtie2-build Guide_RNA_Reference_Genome.fasta grna_ref
bowtie2-build ../GuidesOnly.fasta guides
#bowtie2-build mergedPromoterGuides.fasta merged
bowtie2-build ../PromotersOnly.fasta promoters

while read each; do

  mkdir -p "$each"

  BASENAME1="${each}/Read1"
  BASENAME2="${each}/Read2"

  INDIR="${each}/${each}"

  READFILE1="${each}_R1_001.fastq.gz"
  READFILE2="${each}_R2_001.fastq.gz"


  #FIVEPRIMEADAPTER is the forward sequence of the 5' end of the predicted PCR product and only contains gDNA. This matches the beginning of read 2
  #THREEPRIMEADAPTER is the reverse complement of the 3' end of the predicted PCR product and only contains gDNA. This RC matches the beginning of read 1
  READ1THREEPRIME="ACTAGAATTATACGAGCCGGATGATTAATTGTCAACAGCTCATTTCAGAATATTTGCCAGAACCGGAATTCTTTCAGCTCAGTCGATAGGTAGTAGGCAAGAGTAC"
  READ1FIVEPRIME="CTGTTTCCAGCATAGCTCTTAAAC"

  READ2FIVEPRIME="TACTAGTTTTCTCCTCTTTAGATC"
  READ2THREEPRIME="GAACCGAGGTAACTGGCTTGGTGATAAGCTGTCAAACCAGATCAATTCGCACTAGTCTGCAGACGTAAAAAAAGCGGCGTGGTTAGCCGCT"


  #Dependencies
  #Install cutadapt from the command line by running "pip3 install cutadapt
  #Download FastQC from the following site and put in your Applications folder: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


  #Check quality scores and average read length for the two files.
  /Applications/FastQC.app/Contents/MacOS/fastqc "$READFILE1" "$READFILE2"

  #Remove adapters and genomic DNA sequences
  #Leave the promoter directly at the beginning, and the guide directly at the end of each read
  cutadapt -j 0 -a "$READ1FIVEPRIME"..."$READ1THREEPRIME" -A "$READ2FIVEPRIME"..."$READ2THREEPRIME" -o "$BASENAME1"_trimmed.fastq.gz -p "$BASENAME2"_trimmed.fastq.gz "$READFILE1" "$READFILE2" > "$INDIR"_cutadapt_output_fullset.txt

  #Align reads. -p flag is for number of processors.
  bowtie2 --score-min 'C,0,-1' -p 8 -x guides -U "$BASENAME1"_trimmed.fastq.gz -S "$BASENAME1"_Guides.sam --un "$BASENAME1"_UnalignedGuides.fastq &> "$BASENAME1"_bowtie_output.txt
  bowtie2 --score-min 'C,0,-1' -p 8 -x promoters -U "$BASENAME2"_trimmed.fastq.gz -S "$BASENAME2"_Promoters.sam --un "$BASENAME2"_UnalignedPromoters.fastq &> "$BASENAME2"_bowtie_output.txt
  
  #Compress reads that won't be used further
  gzip -f9 "$BASENAME1"_UnalignedGuides.fastq
  gzip -f9 "$BASENAME2"_UnalignedPromoters.fastq
  


  #Sort the alignments and convert to a binary file
  samtools view -bS "$BASENAME1"_Guides.sam | samtools sort -o "$BASENAME1"_Guides.bam
  samtools view -bS "$BASENAME2"_Promoters.sam | samtools sort -o "$BASENAME2"_Promoters.bam

  #Index and count the resulting alignments
  samtools index "$BASENAME1"_Guides.bam; samtools idxstats "$BASENAME1"_Guides.bam > "$INDIR"_countsR1.txt
  samtools index "$BASENAME2"_Promoters.bam; samtools idxstats "$BASENAME2"_Promoters.bam > "$INDIR"_countsR2.txt

  #Sort sam files by name
  samtools sort -n "$BASENAME1"_Guides.sam -o "$BASENAME1"_Guides--namesort.sam
  samtools sort -n "$BASENAME2"_Promoters.sam -o "$BASENAME2"_Promoters--namesort.sam

  #Pseudo Paired-end
  python ../matchMates.py "$BASENAME1"_Guides--namesort.sam "$BASENAME2"_Promoters--namesort.sam ../GuidesOnly.fasta "${INDIR} Pseudo Paired-end Counts.csv" > "${INDIR} matching.txt"

  #Remove the large files
  rm "$BASENAME1"_Guides--namesort.sam
  rm "$BASENAME2"_Promoters--namesort.sam
  rm "$BASENAME1"_Guides.sam
  rm "$BASENAME2"_Promoters.sam


done <filelist.txt
