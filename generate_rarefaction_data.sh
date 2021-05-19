#Iterate through the sample list
for sample in `cat sample_list.txt`; do
  for readcount in {2500,5000,12500,25000,50000,100000,250000,500000,750000,1000000,1250000,1500000}; do
    #reformat the paired read files to contain a random subset of reads at the specified count
    reformat.sh in1=${sample}_r1.fq in2=${sample}_r2.fq out=temp${readcount}.fq samplereadstarget=${readcount};
    #map reads to the target genome, assumes bowtie2-build has been run first
    bowtie2 -x ../target_genome --interleaved temp${readcount}.fq -S temp${readcount}.sam;
    samtools sort temp${readcount}.sam > temp${readcount}_sorted.sam;
    tempcount=`samtools view -c temp${readcount}_sorted.sam`;
    printf "${sample}\t${tempcount}\t" >> data.txt;
    #generates coverage for the mapping data and collates this information
    samtools coverage temp${readcount}_sorted.sam | cut -f5 | awk '{s+=$1} END {print s}' >> data.txt;
  done
  rm temp*
done
