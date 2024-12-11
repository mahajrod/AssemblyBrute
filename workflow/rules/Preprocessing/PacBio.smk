

#extracthifi - extract HiFi reads (>= Q20) from full CCS reads.bam output

extracthifi -j 10 <input.bam> <output.bam>

#extract fastq from hifi bam
samtools fastq -0 m64241e_230526_215620.hifi.fastq  -@ 4 /maps/projects/codon_0000/data/EBP_DanOME/PB014_Fringillaria_tahapisi/m64241e_230526_215620.hifi_reads.bam &