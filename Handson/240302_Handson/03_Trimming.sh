# DL softwares
wget http://opengene.org/fastp/fastp.0.23.4
chmod a+x ./fastp.0.23.4
# pip install multiqc


#Trim_QC
mkdir -p ./Trimmed_Fastq
mkdir -p ./Trimmed_Report
fastp.0.23.4 -i ./Fastq/ERR204882_1.fastq.gz -I ./Fastq/ERR204882_2.fastq.gz -3 -o ./Trimmed_Fastq/ERR204882_Trimmed_R1.fastq.gz -O ./Trimmed_Fastq/ERR204882_Trimmed_R2.fastq.gz -h ./Trimmed_Report/ERR204882_trimmed_report_fastp.html -j ./Trimmed_Report/ERR204882_trimmed_report_fastp.json -q 30 -n 10 -t 1 -T 1 -l 20 -w 4;


# if you have enough memory
ls ./Fastq/*_1.fastq.gz|xargs -I {} basename {} .fastq.gz|perl -pe s/_1//g| xargs -I {} fastp.0.23.4 -i ./Fastq/{}_1.fastq.gz -I ./Fastq/{}_2.fastq.gz -3 -o ./Trimmed_Fastq/{}_Trimmed_R1.fastq.gz -O ./Trimmed_Fastq/{}_Trimmed_R2.fastq.gz -h ./Trimmed_Report/{}_trimmed_report_fastp.html -j ./Trimmed_Report/{}_trimmed_report_fastp.json -q 30 -n 10 -t 1 -T 1 -l 20 -w 4

# Generalized QC reports
multiqc ./Trimmed_Report/ -o ./Trimmed_Report/;


