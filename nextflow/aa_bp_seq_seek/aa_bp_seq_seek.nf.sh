/*
=================================================================
=================================================================

workflow

=================================================================
=================================================================
*/

workflow aa_bp_seq_seek {

	take:
        aligned_bams

	main:

        println "aa_gain = ${params.aa_gain}"
        println "aa_cnsize_min = ${params.aa_cnsize_min}"
        println "aa_downsample = ${params.aa_downsample}"

        R_ch = Channel.value(file("${params.nf_home}/R/nextflow"))
        
        ch_genome_fasta = Channel.value(file(params.genome_fasta))
        ch_genome_index = Channel.value(file(params.genome_index))
        ch_genome_dict = Channel.value(file(params.genome_dict))
        ch_genome_sa = Channel.value(file(params.genome_sa))
        ch_genome_bwt = Channel.value(file(params.genome_bwt))
        ch_genome_ann = Channel.value(file(params.genome_ann))
        ch_genome_amb = Channel.value(file(params.genome_amb))
        ch_genome_pac = Channel.value(file(params.genome_pac))
        ch_genome_alt = Channel.value(file(params.genome_alt))
        ch_svaba_dbsnp_vcf = Channel.value(file(params.svaba_dbsnp_vcf))


        ch_samples = make_samples ( aligned_bams, params.dna_legacy_file_sheet )
		//OUTPUT: aliquot, bam, bai
		
		ch_input = ch_samples
							.map{aliquot, bam, bai ->
						  return(aliquot)}
		//OUTPUT: aliquot

		
		// process order: 1 & 7 -> 2-1 or 2-2 -> 3 -> 4-1 or 4-2 -> 5 -> 6

		// process 1
        count_aa_amp_num(ch_input)


		// process 2-1 or 2-2
		amp = count_aa_amp_num.out.map{ aliquot, aa_sum_file_path_txt, amp_num_txt, amp_num, cmds ->
							    return[ aliquot,             					    amp_num ] }
        no_aa_amp(amp)
		make_input_for_breakpoints_to_bed(amp, R_ch)
    	

		// process 3
		input_breakpoints_to_bed = make_input_for_breakpoints_to_bed.out
																	.map{ aliquot, bpinput_txt, aainterval_txt, env, cmds -> 
									 							  return[ aliquot, bpinput_txt, aainterval_txt ]}
        get_aa_bp(input_breakpoints_to_bed)


        // process 4-1 or 4-2
        output_breakpoints_to_bed = get_aa_bp.out
    										 .map{ aliquot, bp_bed, env, cmds -> 
    									  return [ aliquot, bp_bed ]}

        convert_aa_bp_to_svaba_target_bed(output_breakpoints_to_bed, R_ch)
		no_aa_bp(output_breakpoints_to_bed)


        // process 5
        input_svaba_local_assembly = ch_samples
    		    							   .combine(convert_aa_bp_to_svaba_target_bed.out, by:0)
											   .map{ aliquot, bam, bai, aa_1kb_bp_svaba_input_bed, env, cmds->
			 							     return[ aliquot, bam, bai, aa_1kb_bp_svaba_input_bed ]}
	    
	    run_svaba_ss_targeted_local_assembly(input_svaba_local_assembly, ch_genome_fasta, ch_genome_index, ch_genome_dict,
             								 ch_genome_sa, ch_genome_bwt, ch_genome_ann, ch_genome_amb, ch_genome_pac,
             								 ch_genome_alt, ch_svaba_dbsnp_vcf)


	    // process 6 (preliminary process)
	    input_align_aa_bp_to_svaba = output_breakpoints_to_bed
        														 .combine(run_svaba_ss_targeted_local_assembly.out, by:0)
        														 .map{ aliquot, bp_bed, svaba_sv_vcf, svaba_unfiltered_sv_vcf, svaba_output, env ->
              											      return [ aliquot, bp_bed, svaba_sv_vcf, svaba_unfiltered_sv_vcf ]}
        
        align_aa_bp_to_svaba(input_align_aa_bp_to_svaba, R_ch)

        // process 7
		copy_script(ch_input, R_ch)
		
}



/*
=================================================================
=================================================================

Define processes

=================================================================
=================================================================
*/

// --------------------------------------------------------------
// 
//
// Tumor-only analysis 
//
//
// --------------------------------------------------------------
process count_aa_amp_num {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/aa_amp_num", pattern: "*", mode: 'copy'

	input: 
		val(aliquot_barcode)

	output:
		tuple val(aliquot_barcode), path("aa_summary_file_path.txt"), path("amplicon_num.txt"), env(amp_num), path("*{command,exitcode}*",hidden:true)

	shell:

		'''

		find  !{params.scratch_dir}/results/!{params.aa_workflow}/!{params.aasuite_ver}/!{params.genome}/minCN!{params.aa_gain}/cnsizeMin!{params.aa_cnsize_min}/!{params.aa_downsample}X/calls/!{aliquot_barcode}/ -name "!{aliquot_barcode}_summary.txt" > ./aa_summary_file_path.txt
		cat ./aa_summary_file_path.txt | xargs -I{} sh -c 'head $0 -n 1 | cut -d " " -f3 > ./amplicon_num.txt' {} ;

		amp_num=`cat ./amplicon_num.txt`

		'''

}

process no_aa_amp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}", pattern: "*", mode: 'copy'

	input: 
		tuple val(aliquot_barcode), val(amp_num)

	output:
		tuple val(aliquot_barcode), path(no_amplicon_flag), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num == "0"
	
	script:

		"""
		#!/bin/bash

		touch no_amplicon_flag

		"""
		

}


process make_input_for_breakpoints_to_bed {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/get_aa_bp", pattern: "*", mode: 'copy'

	label "aa_bp_seq_seek"

	input: 
		tuple val(aliquot_barcode), val(amp_num)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("${aliquot_barcode}_bp.input.txt"), path("${aliquot_barcode}_aa_summary_intervals.txt"), path("*env.txt"), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num != "0"

	script:

		"""
		#!/bin/bash

			Rscript ${R_dir}/input_for_breakpoints_to_bed_script.R \
		    --aa_output_path ${params.scratch_dir}/results/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/calls/${aliquot_barcode} \
			--aliquot_barcode ${aliquot_barcode} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> make_input_for_breakpoints_to_bed_env.txt

		"""

}

process get_aa_bp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/get_aa_bp", pattern: "*", mode: 'copy'

	label "get_aa_bp"

	input:
		tuple val(aliquot_barcode), path(bpinput_txt), path(aainterval_txt)

	output:
		tuple val(aliquot_barcode), path("${aliquot_barcode}_bp.input_breakpoints.bed"), path("*env.txt"), path("*{command,exitcode}*",hidden:true)
	
	
	script:

		"""
		#!/bin/bash

		aainterval=`cat ${aainterval_txt}`

		python3 /home/breakpoints_to_bed.py \
		-i ${bpinput_txt} \
        -r \$aainterval --add_chr_tag && \
        ls -al -R ./ >> get_aa_bp_env.txt

		"""

}

process convert_aa_bp_to_svaba_target_bed {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/get_aa_bp", pattern: "*", mode: 'copy'

	label "aa_bp_seq_seek"

	input:
		tuple val(aliquot_barcode), path(bp_bed)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("${aliquot_barcode}_aa_1kb_bp_svaba_input.bed"), path("*env.txt"), path("*{command,exitcode}*",hidden:true)
	
	when:
		bp_bed.size() > 0
	
	script:

		"""
		#!/bin/bash

			Rscript ${R_dir}/input_for_svaba.R \
		    --breakpoints_bed ${bp_bed} \
			--aliquot_barcode ${aliquot_barcode} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> convert_aa_bp_to_svaba_target_bed_env.txt

		"""

}

process no_aa_bp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}", pattern: "*", mode: 'copy'

	input:
		tuple val(aliquot_barcode), path(bp_bed)

	output:
		tuple val(aliquot_barcode), path(no_aa_bp_flag), path("*{command,exitcode}*",hidden:true)
	
	when:
		bp_bed.size() == 0
	
	script:

		"""
		#!/bin/bash

		touch no_aa_bp_flag

		"""

}

process run_svaba_ss_targeted_local_assembly {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/svaba_ss/local_assembly", pattern: "*", mode: 'copy'

	label "svaba"

	input:
		tuple val(aliquot_barcode), path(bam), path(bai), path(aa_1kb_bp_svaba_input_bed)
		file(genome_fasta)
        file(genome_index)
        file(genome_fasta_dict)
        file(genome_sa) 
        file(genome_bwt) 
        file(genome_ann) 
        file(genome_amb) 
        file(genome_pac) 
        file(genome_alt)
        file(svaba_dbsnp_vcf)

	output:
		tuple val(aliquot_barcode), path("*.svaba.sv.vcf"), path("*unfiltered.sv.vcf"), path("*{.svaba.indel.vcf,txt.gz,contigs.bam,log}"), path("env.txt") 

	
	script:
		"""
    	#!/bin/bash
    	
    	svaba run \
    	-t ${bam} \
    	-p ${params.svaba_threads} \
    	-D ${svaba_dbsnp_vcf} \
    	-a ${aliquot_barcode} \
    	-G ${genome_fasta} \
    	-k ${aa_1kb_bp_svaba_input_bed} && \
    	ls -al -R ./ >> env.txt
    	"""

}

// preliminary step
process align_aa_bp_to_svaba {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/align_aa_bp_to_svaba", pattern: "*", mode: 'copy'

	label "aa_bp_seq_seek"

	input:
		tuple val(aliquot_barcode), path(bp_bed), path(svaba_sv_vcf), path(svaba_unfiltered_sv_vcf)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("*.tsv"), path("env.txt"), path("*{command,exitcode}*",hidden:true)
	
	
	script:

		"""
    	#!/bin/bash
    	
    		Rscript ${R_dir}/matching_aa_svaba.R \
		    --breakpoints_bed ${bp_bed} \
		    --aa_output_path ${params.scratch_dir}/results/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/calls/${aliquot_barcode} \
		    --svaba_sv_vcf ${svaba_sv_vcf} \
		    --svaba_unfiltered_sv_vcf ${svaba_unfiltered_sv_vcf} \
			--aliquot_barcode ${aliquot_barcode} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> env.txt

    	"""

}

process copy_script {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/${aliquot_barcode}/scripts", pattern: "*", mode: 'copy'

	input:
		val(aliquot_barcode)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("*"), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num != "0"
	
	script:

		"""
    	#!/bin/bash
    	
    		cp ${R_dir}/input_for_breakpoints_to_bed_script.R ./`date +%Y.%m.%d`_input_for_breakpoints_to_bed_script.R && \
    		cp ${R_dir}/input_for_svaba.R ./`date +%Y.%m.%d`_input_for_svaba.R && \
    		cp ${R_dir}/matching_aa_svaba.R ./`date +%Y.%m.%d`_matching_aa_svaba.R && \
    		cp ${params.nf_home}/nfs/aa_bp_seq_seek/${params.step}.nf.sh ./`date +%Y.%m.%d`_${params.step}.nf.sh

    	"""

}


/* 
=================================================================
=================================================================

Define functions

=================================================================
=================================================================
*/

// Prepare sample map
def make_samples ( ch_aligned_bams, dna_legacy_file_sheet  ) {
    
    ch_samples = Channel
        .fromPath(dna_legacy_file_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aliquot_barcode,
                            row.patient_barcode,
                            row.source_barcode,
                            row.sequence_type,
                            row.gender,
                            row.tumor_or_normal,
                            row.action ) }
        .filter{      aliquot, patient, source, seqtype, gender, tmnm, action ->
        				action == "run" }
        .map{         aliquot, patient, source, seqtype, gender, tmnm, action ->
        return       [aliquot] }
        .unique()

    ch_aligned_bams
        .combine( ch_samples, by: 0 )
        .map{           aliquot, bam, bai, md5 ->
                return [aliquot, bam, bai ]}

        
}
