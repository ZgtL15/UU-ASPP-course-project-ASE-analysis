python  generate_reads_info_from_sam_4.0.py -s SNP.frq -b ${sample}.sort.bam -l 150 -o ${sample}.reads.info

python distingguish_species_sam.py -b   ${sample}.sort.bam  -i ${sample}.reads.info -r ${sample}.sort.ref.sam -a ${sample}.sort.alt.sam 

samtools view -bS ${sample}.sort.ref.sam > ${sample}.sort.ref.bam
samtools view -bS ${sample}.sort.alt.sam > ${sample}.sort.alt.bam
rm -f ${sample}.sort.ref.sam
rm -f ${sample}.sort.alt.sam
