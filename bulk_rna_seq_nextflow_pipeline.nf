#!/usr/bin/env nextflow

// Bulk RNA-seq pipeline is a sequential pipeline that generates counts from the provided fastq data
// Here are the steps involved in the Bulk RNA-seq Analysis:
// 1. Quality Control - Generate FastQC and MultiQC reports
// 2. Genome Alignment - Map reads to the reference genome using STAR
// 2_1. Mapping Metrics - Generate mapping statistics and quality reports
// 3. Post-Alignment Processing - Filter, deduplicate, and index aligned reads
// 3_1. Quality Assessment - Generate Qualimap reports for aligned reads
// 4. Transcript Assembly and Quantification - Generate counts using StringTie
// 5. Raw Count Generation - Generate raw counts using HTSeq && Feature Counts - Generate counts using Rsubread's featureCounts


println "Bulk RNA-seq Pipeline"
println " "
println "Steps involved in the Bulk RNA-seq Analysis:"
println "1. Quality Control - Generate FastQC and MultiQC reports"
println "2. Genome Alignment - Map reads to the reference genome using STAR"
println "2.1. Mapping Metrics - Generate mapping statistics and quality reports"
println "3. Post-Alignment Processing - Filter, deduplicate, and index aligned reads"
println "3.1. Quality Assessment - Generate Qualimap reports for aligned reads"
println "4. Transcript Assembly and Quantification - Generate counts using StringTie"
println "5. Raw Count Generation - Generate raw counts using HTSeq && Feature Counts"


process generate_fastqc_multiqc_reports {
    // Process to generate FastQC and MultiQC reports
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val fastq_files
    val fastqc_cores

    output:
    val "${config_directory}/fastqc_and_multiqc_reports", emit: fastqc_output

    script:
    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    echo "[INFO] Starting FastQC and MultiQC at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log

    if [ ! -d "${fastq_files}" ]; then
        echo "[ERROR] FASTQ files directory not found: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
        exit 1
    fi
    
    echo '${config_directory}/generate_fastqc_reports.sh ${fastq_files} ${fastqc_cores} ${config_directory}'
    bash ${config_directory}/generate_fastqc_reports.sh ${fastq_files} ${fastqc_cores} ${config_directory}
    
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] FastQC failed on files in ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
        exit 1
    fi
    echo "[INFO] FastQC completed successfully for all FASTQ files" | tee -a ${config_directory}/0_nextflow_logs/fastqc_multiqc.log
    """
}


process star_mapping {
    // Process to perform STAR alignment for RNA-seq data with dynamic file format support
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val samples_file
    val fastq_files
    val is_paired_end
    val read1_suffix
    val read2_suffix
    val num_threads 
    val clip5_num
    val clip3_num
    val index_path
    val file_extension  // New parameter for file extension (e.g., "fastq.gz", "fastq", "fq", "fa", "fasta")
    
    output:
    val "${config_directory}/2_star_mapping_output", emit: mapped_files
    
    script:
    def map_output_fpath = "${config_directory}/2_star_mapping_output"
    def map_log_fpath = "${map_output_fpath}/log"
    
    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    echo "[INFO] Starting STAR mapping at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
    echo "[INFO] Working with file extension: ${file_extension}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
    
    # Validate input files and directories
    if [ ! -f "${samples_file}" ]; then
        echo "[ERROR] Samples file not found: ${samples_file}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        exit 1
    fi
    if [ ! -d "${fastq_files}" ]; then
        echo "[ERROR] FASTQ files directory not found: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        exit 1
    fi
    if [ ! -d "${index_path}" ]; then
        echo "[ERROR] STAR index directory not found: ${index_path}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        exit 1
    fi
    
    mkdir -p ${map_output_fpath} ${map_log_fpath}
    
    # Determine readFilesCommand based on file extension
    case "${file_extension}" in
        *.gz)
            READ_FILES_CMD="zcat"
            echo "[INFO] Using zcat for compressed files (.gz)" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            ;;
        *.bz2)
            READ_FILES_CMD="bzcat"
            echo "[INFO] Using bzcat for bzip2 compressed files (.bz2)" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            ;;
        *.xz)
            READ_FILES_CMD="xzcat"
            echo "[INFO] Using xzcat for xz compressed files (.xz)" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            ;;
        *)
            READ_FILES_CMD="cat"
            echo "[INFO] Using cat for uncompressed files" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            ;;
    esac
    
    # Process each sample
    for sample_name in \$(cat ${samples_file})
    do
        echo "[INFO] Running STAR for sample: \$sample_name" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        
        # Check if sample files exist before processing
        if [ "${is_paired_end}" = true ]; then
            read1_file="${fastq_files}/\${sample_name}${read1_suffix}.${file_extension}"
            read2_file="${fastq_files}/\${sample_name}${read2_suffix}.${file_extension}"
            
            if [ ! -f "\$read1_file" ]; then
                echo "[ERROR] Read 1 file not found: \$read1_file" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
                exit 1
            fi
            if [ ! -f "\$read2_file" ]; then
                echo "[ERROR] Read 2 file not found: \$read2_file" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
                exit 1
            fi
            
            # Paired-end STAR command
            star_cmd="STAR --runThreadN ${num_threads} \\
                      --runMode alignReads \\
                      --genomeDir ${index_path} \\
                      --readFilesIn \$read1_file \$read2_file \\
                      --readFilesCommand \$READ_FILES_CMD \\
                      --clip5pNbases ${clip5_num} ${clip5_num} \\
                      --clip3pNbases ${clip3_num} ${clip3_num} \\
                      --outFileNamePrefix ${map_output_fpath}/\${sample_name} \\
                      --outSAMtype BAM SortedByCoordinate \\
                      --quantMode GeneCounts"
        else
            read1_file="${fastq_files}/\${sample_name}${read1_suffix}.${file_extension}"
            
            if [ ! -f "\$read1_file" ]; then
                echo "[ERROR] Read file not found: \$read1_file" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
                exit 1
            fi
            
            # Single-end STAR command  
            star_cmd="STAR --runThreadN ${num_threads} \\
                      --runMode alignReads \\
                      --genomeDir ${index_path} \\
                      --readFilesIn \$read1_file \\
                      --readFilesCommand \$READ_FILES_CMD \\
                      --clip5pNbases ${clip5_num} \\
                      --clip3pNbases ${clip3_num} \\
                      --outFileNamePrefix ${map_output_fpath}/\${sample_name} \\
                      --outSAMtype BAM SortedByCoordinate \\
                      --quantMode GeneCounts"
        fi
        
        echo "[INFO] Executing STAR command: \${star_cmd}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        eval \${star_cmd} 2> ${map_log_fpath}/\${sample_name}_mapping_step.log
        
        if [ "\$?" -ne 0 ]; then
            echo "[ERROR] STAR mapping failed for sample: \${sample_name}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            exit 1
        fi
        
        # Index the BAM file
        echo "[INFO] Indexing BAM file for sample: \${sample_name}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        index_cmd="samtools index ${map_output_fpath}/\${sample_name}Aligned.sortedByCoord.out.bam"
        echo "[INFO] Executing samtools index command: \${index_cmd}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
        eval \${index_cmd}
        
        if [ "\$?" -ne 0 ]; then
            echo "[ERROR] BAM indexing failed for sample: \${sample_name}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
            exit 1
        fi
        
        echo "[INFO] Successfully processed sample: \${sample_name}" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
    done
    
    echo "[INFO] STAR mapping process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/star_mapping.log
    """
}


process filter_samples {
    // Process to filter and deduplicate aligned samples
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val num_threads
    val blist_bed_file
    val exclude_bed_file
    val filter_samples

    output:
    val "${config_directory}/3_filter_output", emit: filtered_files

    script:
    """
    echo "[INFO] Starting filtering and deduplication at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log

    if [ ! -f "${config_directory}/dedup_and_filtering.sh" ]; then
        echo "[ERROR] Deduplication and filtering script not found: ${config_directory}/dedup_and_filtering.sh" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    fi

    if [ ! -f "${blist_bed_file}" ] && [ ! -f "${exclude_bed_file}" ]; then
        echo "[ERROR] Both blacklist and exclude BED files are missing." | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    fi

    echo "[INFO] Running deduplication and filtering with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Number of threads: ${num_threads}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Blacklist BED file: ${blist_bed_file}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    echo "  Exclude BED file: ${exclude_bed_file}" | tee -a ${config_directory}/0_nextflow_logs/filtering.log

    echo '${config_directory}/dedup_and_filtering.sh $config_directory $num_threads $blist_bed_file $exclude_bed_file' | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    bash ${config_directory}/dedup_and_filtering.sh $config_directory $num_threads $blist_bed_file $exclude_bed_file | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Filtering and deduplication failed" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
        exit 1
    else
        echo "[INFO] Filtering and deduplication completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    fi

    echo "[INFO] Filtering process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/filtering.log
    """
}


process generate_stringtie_counts {
    // Process to generate counts using StringTie
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val num_threads
    val reference_gtf
    val strand_st
    val species
    val filtered_files

    output:
    val "${config_directory}/4_stringtie_counts_output", emit: counts_files

    script:
    """
    echo "[INFO] Starting StringTie counts generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log

    if [ ! -f "${config_directory}/generate_stringtie_counts.sh" ]; then
        echo "[ERROR] StringTie counts script not found: ${config_directory}/generate_stringtie_counts.sh" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    fi

    echo "[INFO] Running StringTie with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Number of threads: ${num_threads}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Strand-specific: ${strand_st}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    echo "  Species: ${species}" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log

    echo 'bash ${config_directory}/generate_stringtie_counts.sh $config_directory $num_threads $reference_gtf $strand_st $species' | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    bash ${config_directory}/generate_stringtie_counts.sh $config_directory $num_threads $reference_gtf $strand_st $species | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] StringTie counts generation failed" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
        exit 1
    else
        echo "[INFO] StringTie counts generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    fi

    echo "[INFO] StringTie counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/stringtie_counts.log
    """
}


process generate_raw_counts {
    // Process to generate raw counts using HTSeq-count
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val reference_gtf
    val is_paired_end
    val strand_hts
    val paired_hts
    val species
    val feature_counts_files

    output:
    val "${config_directory}/5_raw_counts_output", emit: raw_counts_files

    script:
    """
    echo "[INFO] Starting raw counts generation with HTSeq at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log

    if [ ! -f "${config_directory}/generate_raw_counts.sh" ]; then
        echo "[ERROR] HTSeq-count script not found: ${config_directory}/generate_raw_counts.sh" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    fi

    echo "[INFO] Running HTSeq-count with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Strand: ${strand_hts}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    echo "  Paired HTSeq option: ${paired_hts}" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log

    echo 'bash ${config_directory}/generate_raw_counts.sh $config_directory $reference_gtf $is_paired_end $strand_hts $paired_hts $species' | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    bash ${config_directory}/generate_raw_counts.sh $config_directory $reference_gtf $is_paired_end $strand_hts $paired_hts $species | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] HTSeq-count generation failed" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
        exit 1
    else
        echo "[INFO] HTSeq-count generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    fi

    echo "[INFO] Raw counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/htseq_counts.log
    """
}


process generate_feature_counts {
    // Process to generate feature counts using Rsubread's featureCounts
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val num_threads
    val reference_gtf
    val is_paired_end
    val filtered_files

    output:
    val "${config_directory}/5_raw_counts_output", emit: feature_counts_files

    script:
    """
    echo "[INFO] Starting feature counts generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    
    if [ ! -f "${config_directory}/rsubread_featurecount_script.R" ]; then
        echo "[ERROR] Rsubread feature counts script not found: ${config_directory}/rsubread_featurecount_script.R" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] Reference GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    fi

    cd "${config_directory}"
    mkdir -p "${config_directory}/temp_feature_counts_dir"
    cp -r "${config_directory}/3_filter_output"/*.filt.bam "${config_directory}/temp_feature_counts_dir"
    
    echo "[INFO] Running featureCounts with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Reference GTF: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    echo "  Filtered files directory: ${config_directory}/temp_feature_counts_dir" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log

    echo 'Rscript ${config_directory}/rsubread_featurecount_script.R "${config_directory}/temp_feature_counts_dir" $reference_gtf $is_paired_end $num_threads "${config_directory}/5_raw_counts_output/raw_feature_counts.csv"' | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    Rscript ${config_directory}/rsubread_featurecount_script.R "${config_directory}/temp_feature_counts_dir" $reference_gtf $is_paired_end $num_threads "${config_directory}/5_raw_counts_output/raw_feature_counts.csv" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    rm -r "${config_directory}/temp_feature_counts_dir"

    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Feature counts generation failed" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
        exit 1
    else
        echo "[INFO] Feature counts generated successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    fi

    echo "[INFO] Feature counts process completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/feature_counts.log
    """
}


process generate_mapping_metrics {
    // Process to generate mapping metrics and MultiQC report on the mapped files
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val mapped_files

    output:
    val "${config_directory}/5_raw_counts_output", emit: map_metric_files

    script:
    """
    echo "[INFO] Starting mapping metrics generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log

    if [ ! -f "${config_directory}/generate_map_metrics.py" ]; then
        echo "[ERROR] Mapping metrics script not found: ${config_directory}/generate_map_metrics.py" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
        exit 1
    fi

    if [ ! -d "${config_directory}/2_star_mapping_output" ]; then
        echo "[ERROR] STAR mapping output directory not found: ${config_directory}/2_star_mapping_output" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
        exit 1
    fi

    echo "[INFO] Running mapping metrics with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    echo "  Mapped files directory: ${config_directory}/2_star_mapping_output" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log

    echo 'python "${config_directory}/generate_map_metrics.py" "${config_directory}" "${config_directory}/2_star_mapping_output"' | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    python "${config_directory}/generate_map_metrics.py" "${config_directory}" "${config_directory}/2_star_mapping_output" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    
    cd "${config_directory}/2_star_mapping_output"
    multiqc . | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    mkdir -p "${config_directory}/2_1_map_metrics_output_qc"
    scp -pr multiqc_report.html multiqc_data "${config_directory}/2_1_map_metrics_output_qc" | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    rm -r multiqc_report.html multiqc_data | tee -a ${config_directory}/0_nextflow_logs/mapping_metrics.log
    """
}


process generate_qualimap_reports {
    // Process to generate Qualimap reports for quality control of BAM files
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val reference_gtf
    val is_paired_end
    val num_threads
    val filtered_files

    output:
    val "${config_directory}/3_1_qualimap_filter_output_qc", emit: qualimap_files

    script:
    """
    echo "[INFO] Starting Qualimap report generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log

    if [ ! -f "${config_directory}/generate_qualimap_report.sh" ]; then
        echo "[ERROR] Qualimap script not found: ${config_directory}/generate_qualimap_report.sh" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    fi

    if [ ! -f "${reference_gtf}" ]; then
        echo "[ERROR] GTF file not found: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    fi

    echo "[INFO] Running Qualimap with the following parameters:" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Config directory: ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  GTF file: ${reference_gtf}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Paired-end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    echo "  Filtered BAM files: ${filtered_files}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log

    echo "bash "${config_directory}/generate_qualimap_report.sh" "${config_directory}" "${num_threads}" "${reference_gtf}" "${is_paired_end}"" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    bash "${config_directory}/generate_qualimap_report.sh" "${config_directory}" "${num_threads}" "${reference_gtf}" "${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Qualimap report generation failed" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
        exit 1
    else
        echo "[INFO] Qualimap report generation completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/qualimap_report.log
    fi
    """
}


process generate_stats {
    // Process to generate stats about the pipeline
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val fastq_files
    val is_paired_end
    val clip5_num
    val clip3_num
    val feature_counts_files

    script:
    """
    echo "[INFO] Starting stats generation at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log

    if [ ! -f "${config_directory}/generate_stats.sh" ]; then
        echo "[ERROR] generate_stats.sh script not found in ${config_directory}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
        exit 1
    fi

    echo "[INFO] Running stats with the following inputs:" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    echo "  Fastq files: ${fastq_files}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    echo "  Paired end: ${is_paired_end}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    echo "  5' clipping: ${clip5_num}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    echo "  3' clipping: ${clip3_num}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log

    echo "bash ${config_directory}/generate_stats.sh ${config_directory} ${fastq_files} ${is_paired_end} ${clip5_num} ${clip3_num}" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    bash ${config_directory}/generate_stats.sh ${config_directory} ${fastq_files} ${is_paired_end} ${clip5_num} ${clip3_num} | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    
    if [ "\$?" -ne 0 ]; then
        echo "[ERROR] Stats generation failed" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
        exit 1
    else
        echo "[INFO] Stats generation completed successfully at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/6_pipeline_stats.log
    fi
    """
}


process infer_strandedness {
    debug true
    errorStrategy 'terminate'

    input:
    val config_directory
    val threshold

    output:
    tuple val("sample"), val("strandedness_call"), emit: strand_calls

    script:
    """
    mkdir -p "${config_directory}/0_nextflow_logs"
    mkdir -p "${config_directory}/stranded_calls"
    echo "[INFO] Starting strandedness inference at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/infer_strandedness.log
    
    # Process all ReadsPerGene.out.tab files
    for rpg_file in ${config_directory}/2_star_mapping_output/*ReadsPerGene.out.tab
    do
        if [ -f "\$rpg_file" ]; then
            sample_name=\$(basename "\$rpg_file" | sed 's/ReadsPerGene.out.tab//')
            echo "[INFO] Processing sample: \$sample_name" | tee -a ${config_directory}/0_nextflow_logs/infer_strandedness.log
            
            call=\$(python ${config_directory}/infer_strandedness.py "\$rpg_file" -t ${threshold} \\
                    | tail -n1 \\
                    | cut -f2)
            
            echo -e "\${sample_name}\\t\${call}" > ${config_directory}/stranded_calls/\${sample_name}_calls.tsv
            echo "[INFO] Sample \$sample_name: \$call" | tee -a ${config_directory}/0_nextflow_logs/infer_strandedness.log
        fi
    done
    
    # Combine all calls into a single file
    cat ${config_directory}/stranded_calls/*_calls.tsv > ${config_directory}/stranded_calls/all_calls.tsv
    echo "[INFO] Strandedness inference completed at \$(date)" | tee -a ${config_directory}/0_nextflow_logs/infer_strandedness.log
    """
}


workflow {

    // input parameters
    config_directory = params.config_directory
    fastq_files = params.fastq_files
    samples_file = params.samples_file
    fastqc_cores = params.fastqc_cores
    is_paired_end = params.paired_end
    read1_suffix = params.read1_suffix
    read2_suffix = params.read2_suffix
    file_extension = params.file_extension ?: "fastq.gz"
    // star variables
    clip5_num = params.clip5_num
    clip3_num = params.clip3_num ?: 0
    index_path = params.star_index
    // filtering variables
    blist_bed_file = params.blacklist_bed_file_path
    exclude_bed_file = params.exclude_bed_file_path
    // reference file and variables used for stringtie / raw-counts 
    reference_gtf = params.reference_gtf
    // strand_st = params.strand_st
    // strand_hts = params.strand_hts
    paired_hts = params.paired_hts
    threshold = params.threshold ?: 0.8
    species = params.species
    
    // logging
    log.info " "
    log.info "Config directory: ${config_directory}"
    log.info "Fastq files directory: ${fastq_files}"
    log.info "File extension: ${file_extension}"
    log.info " "

    // Assertions to ensure that critical parameters are set
    assert config_directory != null : "ERROR: 'config_directory' must be defined!"
    assert fastq_files != null : "ERROR: 'fastq_files' must be defined!"
    assert samples_file != null : "ERROR: 'samples_file' must be defined!"
    assert fastqc_cores != null : "ERROR: 'fastqc_cores' must be defined!"
    assert is_paired_end != null : "ERROR: 'is_paired_end' must be defined!"
    assert read1_suffix != null : "ERROR: 'read1_suffix' must be defined!"
    assert file_extension != null : "ERROR: 'file_extension' must be defined!"

    // Running FastQC and MultiQC Reports
    if (params.run_fastqc) {
        log.info "Running FastQC and MultiQC"
        println "FastQC and MultiQC output directory: ${config_directory}/fastqc_and_multiqc_reports"
        println "No. of cores to be used for generating FastQC reports: ${fastqc_cores}"
        fastqc_output = generate_fastqc_multiqc_reports(config_directory, fastq_files, fastqc_cores)
        // fastqc_output.view()
    }

    if (params.run_rna_pipeline) {
        assert clip5_num != null : "ERROR: 'clip5_num' must be defined!"
        assert clip3_num != null : "ERROR: 'clip3_num' must be defined!"
        assert index_path != null : "ERROR: 'index_path' must be defined!"
        assert reference_gtf != null : "ERROR: 'reference_gtf' must be defined!"
        assert paired_hts != null : "ERROR: 'paired_hts' must be defined!"
        assert file_extension != null : "ERROR: 'file_extension' must be defined!"
        assert threshold != null : "ERROR: 'threshold' must be defined!"
    }

    // Running STAR Mapping
    if (params.run_rna_pipeline) {
        log.info "Running STAR Mapping"
        println "STAR mapping output directory: ${config_directory}/2_star_mapping_output"
        mapped_files = star_mapping(
            config_directory, samples_file, fastq_files, is_paired_end, read1_suffix, read2_suffix, 
            fastqc_cores, clip5_num, clip3_num, index_path, file_extension
        )

        // generating mapping metrics
        log.info "Generating Mapping Metrics"
        map_metric_files = generate_mapping_metrics(config_directory, mapped_files.mapped_files)
        
        // inferring strandedness
        log.info "Inferring strandedness"
        mapped_files.mapped_files.map { mapping_dir ->
            // 1) pick one ReadsPerGene file
            def rpg = file("${config_directory}/2_star_mapping_output/*ReadsPerGene.out.tab")[0]
            // 2) run the Python inference
            def call = "python ${config_directory}/infer_strandedness.py ${rpg} -t ${threshold}".execute().text.trim().split('\n')[-1].split('\t')[1]
            // 3) map to StringTie / HTSeq flags
            def strand_st_dynamic = (call=='second' ? '--fr' :
                              call=='first'  ? '--rf' : '')
            def strand_hts_dynamic = (call=='second' ? 'yes'  :
                              call=='first'  ? 'reverse' : 'no')
            
            log.info "Detected strandedness: ${call}"
            log.info "StringTie flag: ${strand_st_dynamic}"
            log.info "HTSeq flag: ${strand_hts_dynamic}"
            
            return tuple(mapping_dir, strand_st_dynamic, strand_hts_dynamic)
        }.set { mapping_with_strand }

        // Filtering, deduplication and indexing
        log.info "Running Filtering and Deduplication"
        println "Filtering, deduplication and sorting output directory: ${config_directory}/3_filter_output"
        filtered_files = filter_samples(config_directory, fastqc_cores, blist_bed_file, exclude_bed_file, mapped_files.mapped_files)

        // Use the dynamically determined strand parameters
        mapping_with_strand.map { mapping_dir, strand_st_val, strand_hts_val ->
            // Generating stringtie and raw-counts with dynamic strand parameters
            log.info "Running StringTie and Raw Counts Generation with dynamic strand parameters"
            log.info "StringTie strand parameter: ${strand_st_val}"
            log.info "HTSeq strand parameter: ${strand_hts_val}"
            
            println "Generating StringTie and raw counts in: ${config_directory}/4_stringtie_counts_output and ${config_directory}/5_raw_counts_output"
            stringtie_files = generate_stringtie_counts(config_directory, fastqc_cores, reference_gtf, strand_st_val, species, filtered_files.filtered_files)
            feature_count_files = generate_feature_counts(config_directory, fastqc_cores, reference_gtf, is_paired_end, filtered_files.filtered_files)
            raw_count_files = generate_raw_counts(config_directory, reference_gtf, is_paired_end, strand_hts_val, paired_hts, species, feature_count_files.feature_counts_files)  

            // Generate Stats
            log.info "Generating Pipeline Stats"
            println "Generating Stats for the current run: ${config_directory}/6_pipeline_stats_<TIMESTAMP>.log"
            generate_stats(config_directory, fastq_files, is_paired_end, clip5_num, clip3_num, feature_count_files.feature_counts_files)
            
            return tuple(stringtie_files, feature_count_files, raw_count_files)
        }

        // backup
        // // Filtering, deduplication and indexing
        // log.info "Running Filtering and Deduplication"
        // println "Filtering, deduplication and sorting output directory: ${config_directory}/3_filter_output"
        // filtered_files = filter_samples(config_directory, fastqc_cores, blist_bed_file, exclude_bed_file, mapped_files.mapped_files)

        // // executes rest of the steps right after filtering
        // // filtered_files.branch {
        // //     stringtie: true
        // //     rawcounts: true
        // //     featurecounts: true
        // //     qualimap: true
        // // }.set { branched_filtered_files }
        
        // // Generate Qualimap reports
        // // log.info "Running Qualimap for quality control"
        // // println "Generating qualimap reports for filtered BAM files: ${config_directory}/3_1_qualimap_filter_output_qc"
        // // qualimap_files = generate_qualimap_reports(config_directory, reference_gtf, is_paired_end, fastqc_cores, filtered_files.filtered_files)

        // // Generating stringtie and raw-counts 
        // log.info "Running StringTie and Raw Counts Generation"
        // println "Generating StringTie and raw counts in: ${config_directory}/4_stringtie_counts_output and ${config_directory}/5_raw_counts_output"
        // stringtie_files = generate_stringtie_counts(config_directory, fastqc_cores, reference_gtf, strand_st, species, filtered_files.filtered_files)
        // feature_count_files = generate_feature_counts(config_directory, fastqc_cores, reference_gtf, is_paired_end, filtered_files.filtered_files)
        // raw_count_files = generate_raw_counts(config_directory, reference_gtf, is_paired_end, strand_hts, paired_hts, species, feature_count_files.feature_counts_files)  

        // // Generate Stats
        // log.info "Generating Pipeline Stats"
        // println "Generating Stats for the current run: ${config_directory}/6_pipeline_stats_<TIMESTAMP>.log"
        // generate_stats(config_directory, fastq_files, is_paired_end, clip5_num, clip3_num, feature_count_files.feature_counts_files)

    }
}