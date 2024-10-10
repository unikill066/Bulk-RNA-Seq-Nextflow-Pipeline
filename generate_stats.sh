
#!/bin/bash

set -e
set -o pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

trap 'echo -e "${RED}Error occurred on line $LINENO. Exiting...${NC}"' ERR

# logging function
log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

log_to_file() {
    log -e "[$(date +'%Y-%m-%d %H:%M:%S')] $1" >> "$log_file"
}

check_version() {
    local tool=$1
    local version_cmd=$2
    log -e "\n$tool Version:" >> "$log_file"
    if command -v $version_cmd &> /dev/null; then
        $version_cmd >> "$log_file"
    else
        echo -e "${RED}$tool is not installed.${NC}" >> "$log_file"
    fi
}

usage() {
    log "Usage: $0 <CONFIG_DIRECTORY> <FASTQ_FILES> <PAIRED_END> <CLIP5_NUM> <CLIP3_NUM>"
    exit 1
}


if [ "$#" -ne 5 ]; then
    usage
fi

config_directory=$1
fastq_files=$2
paired_end=$3
clip5_num=$4
clip3_num=$5

if ! [[ "$clip5_num" =~ ^[0-9]+$ ]] || ! [[ "$clip3_num" =~ ^[0-9]+$ ]]; then
    echo -e "${RED}Error: CLIP5_NUM and CLIP3_NUM must be non-negative integers.${NC}"
    exit 1
fi

mkdir -p "${config_directory}" || { echo -e "${RED}Failed to create ${config_directory}${NC}"; exit 1; }

fastqc_multiqc_output_dir="${config_directory}/1_fastqc_and_multiqc_reports"
map_output_fpath="${config_directory}/2_star_mapping_output"
map_metrics_output_fpath="${config_directory}/2_1_map_metrics_output_qc"
filter_output_fpath="${config_directory}/3_filter_output"
qualimap_output_fpath="${config_directory}/3_1_qualimap_filter_output_qc"
counts_output_fpath="${config_directory}/4_stringtie_counts_output"
raw_counts_output_fpath="${config_directory}/5_raw_counts_output"

timestamp=$(date +"%Y%m%d_%H%M%S")
log_file="${config_directory}/6_pipeline_stats_${timestamp}.log"

{
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - Bulk RNA-seq Pipeline Snapshot - x - x - x - x - x - x - x - x - x - x -x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log " "
  log "Configuration Directory: ${config_directory}"
  log "Fastq files: ${fastq_files}"
  log "Paired End: ${paired_end}"
  log "5' Clipping: ${clip5_num}"
  log "3' Clipping: ${clip3_num}"
  log " "
  
  log "Pipeline Stages:"
  log "1. FastQC and MultiQC reports stored at: ${fastqc_multiqc_output_dir}"
  log "2. STAR genome alignment stored at: ${map_output_fpath}"
  log "3. Post-alignment processing (Filtering, Deduplication) stored at: ${filter_output_fpath}"
  log "4. StringTie counts stored at: ${counts_output_fpath}"
  log "5. Raw counts generated using HTSeq & Feature Counts stored at: ${raw_counts_output_fpath}"
  log " "

  log "Checking software versions..."
  
} > "$log_file"


check_version "Nextflow" "nextflow -v"
check_version "Conda" "conda --version"
check_version "FastQC" "fastqc --version"
check_version "MultiQC" "multiqc --version"
check_version "Samtools" "samtools --version"
check_version "Sambamba" "sambamba --version"
check_version "Bedtools" "bedtools --version"
check_version "StringTie" "stringtie --version"
check_version "Qualimap" "qualimap --version"
check_version "HTSeq" "pip show HTSeq"
check_version "R" "R --version"

{
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x"
  log "Pipeline stats saved to: $log_file"
} >> "$log_file"

log "Pipeline stats saved to: $log_file"

exit 0

