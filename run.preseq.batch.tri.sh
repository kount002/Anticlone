#based on run.preseq.sh and modified to run a batch of the samples at once
#scrip runs preseq lc-extrapolate and C-curve


echo "$0" "sript was started" "$(date)"

mkdir -p preseq_results

for i in Processed_data/*/tophat/accepted_hits.bam;
do
(
file=$i
lfile=$(basename $i)
dir=$(dirname $i)
sample=$(echo $i | cut -d '/' -f 2)

##############
# Test to run bed format conversion

samtools view -F 0x0004 $file |	awk '{OFS="\t"; if (and($2, 16)) \
		print $3,$4,$4+length($10),$1,$5,"-"; \
		else print $3,$4,$4+length($10),$1,$5,"+" }' > "$dir"/out_sam.bed
echo "Finished BED conversion of ..." "$file"


sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 "$dir"/out_sam.bed > "$dir"/out.sort_sam.bed
rm "$dir"/out_sam.bed
echo "Finished sorting of ..." "$file"

##############


~/binaries/preseq-1.0.2.Linux_x86_64/preseq lc_extrap -P -o preseq_results/lc.extrap."$sample".txt "$dir"/out.sort_sam.bed

~/binaries/preseq-1.0.2.Linux_x86_64/preseq c_curve -P -v "$dir"/out.sort_sam.bed &> preseq_results/c_curve."$sample".txt

rm "$dir"/out.sort_sam.bed
)&
done
wait

rm preseq_results/lc.summary.txt
for i in preseq_results/lc*
do
echo $i >> preseq_results/lc.summary.txt
tail -1 $i | cut -f2 >> preseq_results/lc.summary.txt
done

echo "The script has finished at" "$(date)"
