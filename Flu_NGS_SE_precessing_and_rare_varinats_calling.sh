#!/bin/bash


WorkDir=/media/nvaulin/Winux/BI/Projects/Prj2_Flu/Flu_project_BI_2022

DataDir=${WorkDir}/data
AlDir=${WorkDir}/alignment
FastQC=${WorkDir}/FastQC_reports
VarDir=${WorkDir}/variant_calling

reference=${DataDir}/Influenza_A_virus_H3N2_segmenе_4_hemagglutinin_gene.fasta

n=0

for f in ${DataDir}/*.fastq.gz
do  
    echo $f
    sample_code=${f##${DataDir}/}
    

    
    if [[ $n -gt 0 ]]
    then
        sample=control_${n}
    else
        sample=person
    fi	
    echo $sample
    
    fastqc -q -t 8 $f -o ${FastQC}
    
    nreads=$(cat ${FastQC}/${sample_code%.fastq.gz}_fastqc.html | grep -Poi "Total Sequences.*?(\d)\D" | grep -o1 "[0-9]*")
    cov=$((${nreads}*151/1690))
    
    bwa mem -t 14 ${reference} $f 2>> ${WorkDir}/alignment_flu_messages.log | samtools view -Sb -@ 13 - | samtools sort -@ 13 - > ${AlDir}/alignment_${sample}_sorted.bam
    
    samtools flagstat ${AlDir}/alignment_${sample}_sorted.bam >> ${WorkDir}/alignment_flu_statistics.log

    
    samtools index ${AlDir}/alignment_${sample}_sorted.bam 
    samtools mpileup -d ${cov} -f ${reference} ${AlDir}/alignment_${sample}_sorted.bam > ${VarDir}/influenza_variants_${sample}.mpileup
 
   java -jar /opt/VarScan/VarScan.v2.4.4.jar   mpileup2snp ${VarDir}/influenza_variants_${sample}.mpileup --min-var-freq 0.001 --variants --output-vcf 1 \
   > ${VarDir}/influenza_${sample}_rare_variants_flu_minfreq1.vcf 2>> ${WorkDir}/varscan_flu_messages.log
   
   awk 'BEGIN{FS=OFS="\t"}{gsub(":","\t",$10)}1' ${VarDir}/influenza_${sample}_rare_variants_flu_minfreq1.vcf | awk 'NR>24 {print $2, $4, $5, $16}' > ${VarDir}/variants_list_${sample}.txt
   
    echo >> ${WorkDir}/alignment_flu_messages.log
    echo >> ${WorkDir}/alignment_flu_statistics.log
    echo >> ${WorkDir}/varscan_flu_messages.log
  
    n=$(($n+1))  
    
done
