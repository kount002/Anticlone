#bash script to make a test sample of the sequenced runs for testing un annotatied search

mkdir -p test_sample/
for i in Processed_data/*
do 
sample=$(echo $i | cut -d "/" -f2)
mkdir -p test_sample/$sample
samtools sort -on Processed_data/$sample/tophat/accepted_hits.bam - | samtools view -h - | head -n 100000 |\
samtools view -S -h -b - > test_sample/$sample/accepted_hits.bam 
done

