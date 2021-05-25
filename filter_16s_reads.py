#read filtering

from Bio import SeqIO
import os

log = open('output_log.txt','w')

print('Beginning read filtering.', file=log)
print('-'*15, file=log)

meta_reads_file = '[community profile reads file]'
reads_dir = '[read files to filter]'
outfile = 'all_reads_filtered.txt'

meta_reads=[]
for record in SeqIO.parse(meta_reads_file,'fastq'): meta_reads.append(str(record.seq))

print('Community contains ' + str(len(meta_reads)) + ' total ASVs', file=log)

meta_reads = set(meta_reads)

print('Community contains ' + str(len(meta_reads)) + ' unique ASVs', file=log)
print('-'*15, file=log)

good_samples = []
for sample in os.listdir(reads_dir):
    sample_reads = {} #need to count these
    total = 0
    for record in SeqIO.parse(reads_dir + '/' + sample, 'fastq'): 
        sample_reads[str(record.seq)] = sample_reads.get(str(record.seq), 0) + 1
        total += 1
    
    sample_reads = sorted([(x,float(sample_reads[x])/total) for x in sample_reads], key=lambda asv: asv[1], reverse=True)
    print('Sample ' + sample + ' contains ' + str(len(sample_reads)) + ' unique ASVs ( ' + str(total) + 'total reads)' , file=log)
    print('Top ASV accounts for ' + str(sample_reads[0][1]) + ' percent of all ASVs', file=log)
    
    if sample_reads[0][1] < 0.7:
        print('Top ASV does not pass 70% threshold.', file=log)
        pass
    else:
        print('Top ASV passes 70% threshold, checking against community...', file=log)
        if sample_reads[0][0] in meta_reads:
            print('Top ASV is in community, sample is good.', file=log)
            good_samples.append(sample)
        else:
            print('Top ASV is not found in community, sample is bad.', file=log)
            pass
    
    print('-'*15, file=log)
    
print('Of ' + str(len(os.listdir(reads_dir))) + ' samples, ' + str(len(good_samples)) + ' passed filtering.', file=log)

log.close()

with open(outfile,'w') as f: 
    for entry in good_samples: f.write(entry + '\n')
