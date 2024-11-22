#!/bin/bash

set -e
set -o pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

trap 'echo -e "${RED}Error occurred on line $LINENO. Exiting...${NC}"' ERR

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

usage() {
    log "Usage: $0 <CONFIG_DIRECTORY> <NUM_THREADS> <REFERENCE_GTF> <PAIRED_END (true/false)>"
    exit 1
}

if [ "$#" -ne 4 ]; then
    usage
fi

config_directory=$1
num_threads=$2
ref_genome_file=$3
paired_end=$3

memory_size="10G"
available_mem=$(free -g | awk '/^Mem:/{print $2}')

# if [[ "$available_mem" -gt 10 ]]; then
#     memory_size="${available_mem}G"
#     log "Memory size dynamically set to: $memory_size"
# else
#     log "Using default memory size: $memory_size"
# fi

# INPUT & OUTPUT PATHS
filter_output_fpath="${config_directory}/3_filter_output"
qualimap_output_fpath="${config_directory}/3_1_qualimap_filter_output_qc"

mkdir -p "${qualimap_output_fpath}"
cd "${qualimap_output_fpath}"

# Check if reference GTF file exists
if [ ! -f "$ref_genome_file" ]; then
    echo -e "${RED}Error: Reference GTF file not found: $ref_genome_file${NC}"
    exit 1
fi

if [ ! -d "${filter_output_fpath}" ]; then
    echo -e "${RED}Error: Filtered files directory not found: $filter_output_fpath${NC}"
    exit 1
fi

if [ -d "${filter_output_fpath}" ]; then
    for file in "$filter_output_fpath"/*.filt.bam; do
        log " "
        log "Running Qualimap on $file"
        
        if [ "$paired_end" = true ]; then
            log "Running Qualimap for paired-end data"
            log "qualimap rnaseq -pe -bam $file -outdir $qualimap_output_fpath -nt $num_threads -gtf $ref_genome_file --java-mem-size=$memory_size"
            qualimap rnaseq -pe -bam "$file" -outdir "$qualimap_output_fpath" -nt "$num_threads" -gtf "$ref_genome_file" --java-mem-size="$memory_size"
            # stackoverflow.com/questions/45964751/python3-opencv-cv2-error-215-unable-to-show-captured-image (-nt threads (8) documentation)
        else
            log "Running Qualimap for single-end data"
            log "qualimap rnaseq -bam $file -outdir $qualimap_output_fpath -nt $num_threads -gtf $ref_genome_file --java-mem-size=$memory_size"
            qualimap rnaseq -bam "$file" -outdir "$qualimap_output_fpath" -nt "$num_threads" -gtf "$ref_genome_file" --java-mem-size="$memory_size"
        fi
        log " "
    done
else
    log "Filtered files directory not found: $filter_output_fpath"
    exit 1
fi

log "Qualimap output directory: $qualimap_output_fpath"
log "Qualimap pipeline completed successfully."

exit 0