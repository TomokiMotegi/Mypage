### Run RSEM (calculation of expression counts)
ls ./chr22_Fastq/*_alignable_1.fq |xargs -I {} basename {} .fq |perl -pe s/_alignable_1//g|sort|uniq|xargs -I {} mkdir -p {}
#　--hisat2-hcaは --no-softclipなのでfastqのadapter trimmingが必要になる
# 16 threads 256 GB in 8 min
./RSEM-1.3.3/rsem-calculate-expression --paired-end --hisat2-hca \
 --hisat2-path ./hisat2-2.2.1 \
 --estimate-rspd --calc-ci --ci-memory 4096 --append-names --output-genome-bam \
 ./chr22_Fastq/ERR204882_alignable_1.fq ./chr22_Fastq/ERR204882_alignable_2.fq \
 ./hg38_hisat2/hg38_chr22 ./ERR204882/ERR204882;

# if you have enough memory and threads try as below
ls ./chr22_Fastq/*_alignable_1.fq |xargs -I {} basename {} .fq |perl -pe s/_alignable_1//g|sort|uniq|xargs -I {} mkdir -p {}
ls ./chr22_Fastq/*_alignable_1.fq |xargs -I {} basename {} .fq |perl -pe s/_alignable_1//g|xargs -I {} bash -c "nohup ./RSEM-1.3.3/rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /restricted/projectnb/camplab/home/tmotegi/Statistics/RNA_seq_Training/hisat2-2.2.1 --estimate-rspd --calc-ci --ci-memory 4096 --append-names --output-genome-bam ./chr22_Fastq/{}_alignable_1.fq ./chr22_Fastq/{}_alignable_2.fq ./hg38_hisat2/hg38_chr22 ./{}/{} &"


### add sample name and data collection ###
ls ./ERR*/ERR*.genes.results|xargs -I {} basename {} .genes.results|xargs -I {} sed -i -e "s/expected_count\t/expected_count_{}\t/g" ./{}/{}.genes.results
ls ./ERR*/ERR*.genes.results|xargs -I {} basename {} .genes.results|xargs -I {} sed -i -e "s/TPM\t/TPM_{}\t/g" ./{}/{}.genes.results
paste ./ERR*/ERR*.genes.results > ./all.genes.results;
awk '{print$1"\t"$2"\t"$3"\t"$5"\t"$22"\t"$39"\t"$56"\t"$73"\t"$90"\t"$107}' all.genes.results > expected_count_all.genes.results;
awk '{print$1"\t"$2"\t"$3"\t"$6"\t"$23"\t"$40"\t"$57"\t"$74"\t"$91"\t"$108}' all.genes.results > TPM_all.genes.results;

