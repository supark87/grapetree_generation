#!/usr/bin/env nextflow
params.reads = params.input.fastq_path+'*{_1,_2}*.fastq'
params.ref = "$baseDir/ref/TruSeq3-PE.fa"
params.ref2 = "$baseDir/ref/GCF_000006945.2_ASM694v2_genomic.fna"
//params.ref3 = "$baseDir/ref/midas_exports.sh"
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
out_path = Channel.fromPath(params.output.folder)
pyscripts="$baseDir/pyscripts"

process combineFastq {
    publishDir "$params.output.folder/trimFastq/${pair_id}", pattern: "*.fastq", mode : "copy"
    input:
        //set val(pair_id), path(fastq_group) from fastq_path
		set val(pair_id), path(fastq_group) from fastq_path
    output:
        //tuple val(pair_id), path("${pair_id}_R1.fastq.gz"), path("${pair_id}_R2.fastq.gz") into comb_out
		tuple val(pair_id), path("${pair_id}_R1.fastq"), path("${pair_id}_R2.fastq") into comb_out

    script:
        """

		cat *_1*.fastq > ${pair_id}_R1.fastq
		cat *_2*.fastq >  ${pair_id}_R2.fastq

        """
}
comb_out.into{preqc_path;trim_path}
// process midas {
//     publishDir "$params.output.folder/midas", mode : "copy"
//     input:
//         set val(sample), path(read_one), path(read_two) from midas_path
// 		path ref from params.ref3
    
//     output:
//         tuple val(sample), path("*") into midas_out
    
//     script:
//         """
//         run_midas.py species ${sample} -1 ${read_one} -2 ${read_two}
//         """
// }


process preFastQC {
	// conda "$baseDir/envs/fastqc.yml"
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
		module load fastqc/0.11.5
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
		module purge
        """
}

// process shovill{
// 	publishDir "$params.output.folder/shovill/${sample}", mode : "copy"
// 	errorStrategy 'ignore'
//     input:
//         set val(sample), path(read_one), path(read_two) from trim_path
//         path ref from params.ref
    
//     output:
//         tuple val(sample),  path("${sample}_contigs.fa") into shovill_out
// 		set val(sample), path(${sample}_R1_trim.fq.gz), path(${sample}_R2_trim.fq.gz) into trim_seq
// 		//tuple val(sample), path(${sample}_R2_trim.fq.gz) into trim_seq

//     script:
//         """          
//     	shovill --outdir ${sample} -R1 ${read_one} -R2 ${read_two} --ram 8 --force
// 		mv ${sample}/contigs.fasta ${sample}_contigs.fasta
// 		mv ${sample}/R1.cor.fq.gz ${sample}_R1_trim.fq.gz
// 		mv ${sample}/R2.cor.fq.gz ${sample}_R2_trim.fq.gz
//         """
// }

// shovill_out.into{shovill_out1; shovill_out2; shovill_out3}
// trim_seq.into{trim_seq1,trim_seq2}


// process postFastQC {
//     publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
//     input:
//         set val(sample), path(read_one), path(read_two) from trim_seq1
//     output:
//         tuple val(sample), path("*") into postqc_out
    
//     script:
//         """
// 		module load fastqc/0.11.5
//         fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
// 		module purge
//         """
// }

// // process quast{
// // 	// conda "$baseDir/envs/quast.yml"
// // 	publishDir "$params.output.folder/quast/${sample}", mode : "copy"
// // 	input:
// // 		tuple val(sample), path("${sample}_contigs.fa") from shovill_out2
// // 		set val(sample), path(read_one), path(read_two) from assembly_path from trim_seq2
// // 		path ref from params.ref2

		
// // 	output:
// // 		tuple val(sample), path("*") into quast_out

// // 	script:
// // 	"""
// // 	module load quast/5.0.2
// // 	quast.py ${sample}_contigs.fasta -1 ${read_one} -2 ${read_two} -o ${sample} -r ${ref}
// // 	mv ${sample}/report.html ${sample}_report.html
// // 	module purge
// // 	"""
	
// // }

// // process mashtree{
// // 	publishDir "$params.output.folder/mashtree/${sample}", mode : "copy"
// // 	input:
// // 		tuple val(sample), path("${sample}_contigs.fa") from shovill_out2
// // 	output:
// // 		tuple val(sample), path("${sample}_tree.dnd") into quast_out

// // 	script:
// // 	"""
// // 	module load mashtree/1.0.4
// // 	mashtree "${sample}_contigs.fa" > ${sample}_tree.dnd
// // 	"""

// // }


// // process resfinder{
// // 	// conda "$baseDir/envs/resfinder.yml"
// // 	publishDir "$params.output.folder/resfinder/${sample}", mode : "copy"
// // 	input:
// // 		tuple val(sample), path("${sample}_scaffolds.fasta"), path("${sample}_contigs.fasta") from spades_out2

// // 	output:
// // 		tuple val(sample), path("${sample}")

// // 	script:
// // 	"""
// // 	module load abricate/1.0.1
// // 	abricate --db resfinder ${sample}_scaffolds.fasta > ${sample}
// // 	module purge
// // 	"""
// // }

// // process seqsero{

// // 	publishDir "$params.output.folder/seqsero/${sample}", mode : "copy"
// // 	input:
// // 		tuple val(sample), path("${sample}_scaffolds.fasta"), path("${sample}_contigs.fasta") from spades_out3

// // 	output:
// // 		path("${sample}")
		
// // 	script:
// // 	"""
// // 	SeqSero2_package.py -t 4 -i ${sample}_scaffolds.fasta -d ${sample} -p 10 -m k
// // 	"""

// // }
	
	
// // process mummer{
// // 	publishDir "$params.output.folder/mummer/${sample}", mode : "copy"
// // 	input:
// // 		tuple val(sample), path("${sample}_scaffolds.fasta"), path("${sample}_contigs.fasta") from spades_out4
// // 		path ref from params.ref2

// // 	output:
// // 		path("${sample}.delta")

// // 	script:
// // 	"""
// // 	module load MUMmer/4
// // 	nucmer -p ${sample} ${sample}_scaffolds.fasta ${ref}
// // 	module purge
// // 	"""
// // }	
	
