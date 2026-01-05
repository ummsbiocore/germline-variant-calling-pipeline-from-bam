$HOSTNAME = ""
params.outdir = 'results'  


if (!params.firstRound){params.firstRound = ""} 
if (!params.secondRound){params.secondRound = ""} 
if (!params.genome){params.genome = ""} 
if (!params.bams){params.bams = ""} 

Channel.value(params.firstRound).set{g_7_2_g_56}
Channel.value(params.secondRound).set{g_21_2_g_57}
g_48_0_g_55 = file(params.genome, type: 'any')
if (params.bams){
   Channel.fromPath(params.bams, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_58_0_g_35}
} else {
	g_58_0_g_35 = Channel.empty()
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 64
}
//* platform
//* platform
//* autofill

process markDuplicatesSpark {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_dedup_metrics.txt$/) "metrics/$filename"}
input:
 tuple val(name), file(reads)

output:
 tuple val(name), file("${name}_sorted_dedup.bam")  ,emit:g_35_mapped_reads00_g_4 
 path "${name}_dedup_metrics.txt"  ,emit:g_35_txtFile11 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
mkdir -p tmp/${name}
gatk --java-options '-Djava.io.tmpdir=tmp/${name}' \
 MarkDuplicates \
-I ${reads} \
-M ${name}_dedup_metrics.txt \
-O ${name}_sorted_dedup.bam 
rm -r tmp/${name}
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process build_gatk4_genome_dictionary {

input:
 path genome

output:
 path "genomeDict"  ,emit:g_55_genomeDict01_g_4 

container 'quay.io/viascientific/gatk:1.0.0'
stageInMode 'copy'

script:

basename = genome.baseName

"""
mkdir -p genomeDict && mv ${genome} genomeDict/ && cd genomeDict
samtools faidx ${genome}
gatk CreateSequenceDictionary -R ${genome} 
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 2
}
//* platform
//* platform
//* autofill

process getMetrics {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_alignment_metrics.txt$/) "metrics/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_insert_metrics.txt$/) "metrics/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_insert_size_histogram.pdf$/) "report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_depth_out.txt$/) "metrics/$filename"}
input:
 tuple val(name),file(sorted_dedup_reads)
 path ref

output:
 path "${name}_alignment_metrics.txt"  ,emit:g_4_txtFile00 
 path "${name}_insert_metrics.txt" ,optional:true  ,emit:g_4_txtFile11 
 path "${name}_insert_size_histogram.pdf" ,optional:true  ,emit:g_4_outputFilePdf22 
 path "${name}_depth_out.txt"  ,emit:g_4_txtFile33 

container 'quay.io/viascientific/variant_calling:1.0'

errorStrategy 'retry'
maxRetries 1

when:
(params.run_Metrics && (params.run_Metrics == "yes")) || !params.run_Metrics

script:
"""
indexname=\$(ls ${ref}/*.fa )
    picard \
       CollectAlignmentSummaryMetrics \
	   R=\$indexname \
       I=${sorted_dedup_reads} \
	   O=${name}_alignment_metrics.txt
    picard \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	    OUTPUT=${name}_insert_metrics.txt \
        HISTOGRAM_FILE=${name}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${name}_depth_out.txt
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 8
}
//* platform
//* platform
//* autofill

process HaplotypeCaller {

input:
 tuple val(name), file(input_bam)
 path ref
 val round

output:
 tuple val(name), val(round), file("${input_bam}"), file("${name}_raw_variants_${round}.vcf")  ,emit:g_56_VCFset00_g_8 

container 'quay.io/viascientific/gatk:1.0.0'

script:

"""
samtools index $input_bam
indexname=\$(ls ${ref}/*.fa )
 gatk HaplotypeCaller \
	-R \${indexname} \
	-I $input_bam \
	-O ${name}_raw_variants_${round}.vcf
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process selectVariants {

input:
 tuple val(name), val(round),  file(inputbam), file(variants)
 path ref

output:
 tuple val(name), val(round),  file("${inputbam}"), file("${name}_raw_*_${round}.vcf*")  ,emit:g_8_VCFset00_g_9 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
indexname=\$(ls ${ref}/*.fa )
gatk SelectVariants \
	-R \${indexname} \
	-V ${variants} \
	-select-type SNP \
	-O ${name}_raw_snps_${round}.vcf
gatk SelectVariants \
	-R \${indexname} \
	-V ${variants} \
	-select-type INDEL \
	-O ${name}_raw_indels_${round}.vcf
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process VariantFiltration {

input:
 tuple val(name), val(round), file(inputbam), file(variants)
 path ref

output:
 tuple val(name), val(round), file("${inputbam}"),  file("${name}_filtered_*_${round}.vcf*")  ,emit:g_9_VCFset00_g_14 
 tuple val(name), file("${name}_filtered_*_${round}.vcf*")  ,emit:g_9_VCFset11 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
#FOR SNPS
indexname=\$(ls ${ref}/*.fa )
gatk VariantFiltration \
	-R \${indexname} \
	-V ${name}_raw_snps_${round}.vcf \
	-O ${name}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#FOR INDELS
gatk VariantFiltration \
        -R \${indexname} \
        -V ${name}_raw_indels_${round}.vcf \
        -O ${name}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 2
}
//* platform
//* platform
//* autofill

process BaseRecalibrator {

input:
 tuple val(name), val(round), file(input_bam), file(filtered_variants)
 path ref

output:
 tuple val(name), file("${input_bam}"), file("${name}_bqsr_*.vcf*"), file("${name}_recal_data.table")  ,emit:g_14_VCFset00_g_16 
 val basecal_parameters  ,emit:g_14_run_parameters12_g_16 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
basecal_parameters = params.BaseRecalibrator.basecal_parameters


"""
indexname=\$(ls ${ref}/*.fa )
gatk SelectVariants \
	--exclude-filtered \
	-V ${name}_filtered_snps_${round}.vcf \
	-O ${name}_bqsr_snps.vcf

gatk SelectVariants \
	--exclude-filtered \
	-V ${name}_filtered_indels_${round}.vcf \
	-O ${name}_bqsr_indels.vcf
	
gatk BaseRecalibrator \
	-R \${indexname}  ${basecal_parameters} \
	-I ${input_bam} \
	--known-sites ${name}_bqsr_snps.vcf \
    --known-sites ${name}_bqsr_indels.vcf \
	-O ${name}_recal_data.table

"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 5
}
//* platform
//* platform
//* autofill

process applyBSQRS {

input:
 tuple val(name), file(input_bam), file(variants), file(recal_table)
 path ref
 val basecal_parameters

output:
 tuple val(name),  file("${name}_recal_data.table"),  file("${name}_post_recal_data.table")  ,emit:g_16_outputFileTab00_g_17 
 tuple val(name), file("${name}_recal.bam")  ,emit:g_16_mapped_reads10_g_57 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
"""
indexname=\$(ls ${ref}/*.fa )
	gatk ApplyBQSR \
        -R \$indexname \
        -I $input_bam \
        -bqsr ${name}_recal_data.table \
        -O ${name}_recal.bam
    gatk BaseRecalibrator \
	    -R \$indexname ${basecal_parameters} \
		-I ${name}_recal.bam \
	    --known-sites ${name}_bqsr_snps.vcf\
		--known-sites ${name}_bqsr_indels.vcf \
		-O ${name}_post_recal_data.table
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process AnalyzeCovariates {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_recalibration_plots.pdf$/) "AnalyzeCovarage/$filename"}
input:
 tuple val(name),file(recal_table), file(post_recal_table)

output:
 path "${name}_recalibration_plots.pdf"  ,emit:g_17_outputFilePdf00 

container 'quay.io/viascientific/variant_calling:1.0'

errorStrategy 'retry'
maxRetries 1

when:
(params.run_AnalyzeCoverage && (params.run_AnalyzeCoverage == "yes")) || !params.run_AnalyzeCoverage
script:

"""
gatk AnalyzeCovariates \
	-before $recal_table \
	-after $post_recal_table \
	-plots ${name}_recalibration_plots.pdf
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 8
}
//* platform
//* platform
//* autofill

process CalibHaplotypeCaller {

input:
 tuple val(name), file(input_bam)
 path ref
 val round

output:
 tuple val(name), val(round), file("${input_bam}"), file("${name}_raw_variants_${round}.vcf")  ,emit:g_57_VCFset00_g_22 

container 'quay.io/viascientific/gatk:1.0.0'

script:

"""
samtools index $input_bam
indexname=\$(ls ${ref}/*.fa )
 gatk HaplotypeCaller \
	-R \${indexname} \
	-I $input_bam \
	-O ${name}_raw_variants_${round}.vcf
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process selectPostSNPs {

input:
 tuple val(name), val(round),  file(inputbam), file(variants)
 path ref

output:
 tuple val(name), val(round),  file("${inputbam}"), file("${name}_raw_*_${round}.vcf*")  ,emit:g_22_VCFset00_g_24 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
indexname=\$(ls ${ref}/*.fa )
gatk SelectVariants \
	-R \${indexname} \
	-V ${variants} \
	-select-type SNP \
	-O ${name}_raw_snps_${round}.vcf
gatk SelectVariants \
	-R \${indexname} \
	-V ${variants} \
	-select-type INDEL \
	-O ${name}_raw_indels_${round}.vcf
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 2
    $MEMORY = 4
}
//* platform
//* platform
//* autofill

process SNPPostFIltration {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_filtered_.*_${round}.vcf.*$/) "filtered_snps/$filename"}
input:
 tuple val(name), val(round), file(inputbam), file(variants)
 path ref

output:
 tuple val(name), val(round), file("${inputbam}"),  file("${name}_filtered_*_${round}.vcf*")  ,emit:g_24_VCFset00 
 tuple val(name), file("${name}_filtered_*_${round}.vcf*")  ,emit:g_24_VCFset11 

container 'quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0'

script:
	
"""
#FOR SNPS
indexname=\$(ls ${ref}/*.fa )
gatk VariantFiltration \
	-R \${indexname} \
	-V ${name}_raw_snps_${round}.vcf \
	-O ${name}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
#FOR INDELS
gatk VariantFiltration \
        -R \${indexname} \
        -V ${name}_raw_indels_${round}.vcf \
        -O ${name}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"

"""
}


workflow {


markDuplicatesSpark(g_58_0_g_35)
g_35_mapped_reads00_g_4 = markDuplicatesSpark.out.g_35_mapped_reads00_g_4
(g_35_mapped_reads00_g_56) = [g_35_mapped_reads00_g_4]
g_35_txtFile11 = markDuplicatesSpark.out.g_35_txtFile11


build_gatk4_genome_dictionary(g_48_0_g_55)
g_55_genomeDict01_g_4 = build_gatk4_genome_dictionary.out.g_55_genomeDict01_g_4
(g_55_genomeDict01_g_8,g_55_genomeDict01_g_9,g_55_genomeDict01_g_14,g_55_genomeDict01_g_16,g_55_genomeDict01_g_56,g_55_genomeDict01_g_57,g_55_genomeDict01_g_24,g_55_genomeDict01_g_22) = [g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4,g_55_genomeDict01_g_4]


getMetrics(g_35_mapped_reads00_g_4,g_55_genomeDict01_g_4)
g_4_txtFile00 = getMetrics.out.g_4_txtFile00
g_4_txtFile11 = getMetrics.out.g_4_txtFile11
g_4_outputFilePdf22 = getMetrics.out.g_4_outputFilePdf22
g_4_txtFile33 = getMetrics.out.g_4_txtFile33


HaplotypeCaller(g_35_mapped_reads00_g_56,g_55_genomeDict01_g_56,g_7_2_g_56)
g_56_VCFset00_g_8 = HaplotypeCaller.out.g_56_VCFset00_g_8


selectVariants(g_56_VCFset00_g_8,g_55_genomeDict01_g_8)
g_8_VCFset00_g_9 = selectVariants.out.g_8_VCFset00_g_9


VariantFiltration(g_8_VCFset00_g_9,g_55_genomeDict01_g_9)
g_9_VCFset00_g_14 = VariantFiltration.out.g_9_VCFset00_g_14
g_9_VCFset11 = VariantFiltration.out.g_9_VCFset11


BaseRecalibrator(g_9_VCFset00_g_14,g_55_genomeDict01_g_14)
g_14_VCFset00_g_16 = BaseRecalibrator.out.g_14_VCFset00_g_16
g_14_run_parameters12_g_16 = BaseRecalibrator.out.g_14_run_parameters12_g_16


applyBSQRS(g_14_VCFset00_g_16,g_55_genomeDict01_g_16,g_14_run_parameters12_g_16)
g_16_outputFileTab00_g_17 = applyBSQRS.out.g_16_outputFileTab00_g_17
g_16_mapped_reads10_g_57 = applyBSQRS.out.g_16_mapped_reads10_g_57


AnalyzeCovariates(g_16_outputFileTab00_g_17)
g_17_outputFilePdf00 = AnalyzeCovariates.out.g_17_outputFilePdf00


CalibHaplotypeCaller(g_16_mapped_reads10_g_57,g_55_genomeDict01_g_57,g_21_2_g_57)
g_57_VCFset00_g_22 = CalibHaplotypeCaller.out.g_57_VCFset00_g_22


selectPostSNPs(g_57_VCFset00_g_22,g_55_genomeDict01_g_22)
g_22_VCFset00_g_24 = selectPostSNPs.out.g_22_VCFset00_g_24


SNPPostFIltration(g_22_VCFset00_g_24,g_55_genomeDict01_g_24)
g_24_VCFset00 = SNPPostFIltration.out.g_24_VCFset00
g_24_VCFset11 = SNPPostFIltration.out.g_24_VCFset11


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
