#!/bin/bash
# takes a list of genes and creates a file for IGV (integrated genomic viewer) 
# synthasis: igvprep.sh [gene_list] [output_directory]
# gene list should contain one gene per line
# output dir is optional default is igv_out_files


#check that input list present
#check that output directory present or use default


if [ -z "$1" ];
then
echo "No path to files suppled, exiting"
exit
fi


#directory with sam files
sams='clone_count'
outdir='igv_out_files'


if [ -n "$2" ];
then
outdir="$2"
fi



mkdir -p "$outdir"

#sort files
for i in "$sams"/*.sam ; do
    (
    si=$(basename $i)
    echo "Sample in progress: $i"
    #head -n1 "$i" > "$outdir"/"$si"
    while IFS='' read -r line ; do
        echo 'Parsing for gene:' "$line"
        cat "$i" | grep $line >> "$outdir"/"$si"
    done <"$1"
    
    /home/kount002/binaries/IGVTools/igvtools sort "$outdir"/"$si" "$outdir"/"$si"_temp.sort && \
    
    samtools view -Sb -t /nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai \
             "$outdir"/"$si"_temp.sort > "$outdir"/"$si".bam
 
    rm "$outdir"/"$si"_temp.sort "$outdir"/"$si"
    #mv "$outdir"/"$si"_temp.sort "$outdir"/"$si"
    /home/kount002/binaries/IGVTools/igvtools count "$outdir"/"$si".bam "$outdir"/"$si".tdf $HOME/binaries/IGVTools/genomes/sizes/hg19.chrom.sizes
    
    /home/kount002/binaries/IGVTools/igvtools index "$outdir"/"$si".bam
    )&
done
wait
echo 'Im done'

