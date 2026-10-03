# download
prefetch --option-file SRR_Acc_List.txt -O ./sra --max-size u
# decompress
while read acc; do
  fasterq-dump ./sra/$acc --split-files --include-technical -e 8 -O ./fastq
done < SRR_Acc_List.txt
# rename
lane=1 
for srr in SRR18028{398..417}; do
L=$(printf "L%03d" $lane);
ln -sf /data/fastq2/${srr}_1.fastq.gz pool1_S1_${L}_I1_001.fastq.gz;
ln -sf /data/fastq2/${srr}_2.fastq.gz pool1_S1_${L}_R1_001.fastq.gz;
ln -sf /data/fastq2/${srr}_3.fastq.gz pool1_S1_${L}_R2_001.fastq.gz;
lane=$((lane+1));
done
# cell ranger
cellranger count \
  --id=pool1 \
  --create-bam=true \
  --transcriptome=/home/mw/refdata-cellranger-GRCh38-3.0.0 \
  --fastqs=/data/fastq2/ \
  --sample=pool1