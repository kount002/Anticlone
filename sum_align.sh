#compiles result of tophat assambly in one log file for easy review

rm -f sum_align.log
for i in Processed_data/*/tophat/align_summary.txt
do
name=$(echo $i | cut -d "/" -f 2) 
echo $name "##################" >> sum_align.log
cat $i >> sum_align.log
echo -e "\n" >> sum_align.log
done

#aggregate log for cutadapt

for i in Processed_data/*/*cutadapt.log
do
name=$(echo $i | cut -d "/" -f 2)
echo $name "######################" >> sum_align.log
cat $i >> sum_align.log
echo -e "\n" >> sum_align.log
done
