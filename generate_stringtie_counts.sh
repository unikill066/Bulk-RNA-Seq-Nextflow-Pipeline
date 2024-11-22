#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No color

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

usage() {
    echo -e "${RED}Usage:${NC} $0 <CONFIG_DIRECTORY> <NUM_THREADS> <REFERENCE_GTF> <STRAND_ST> <SPECIES>"
    exit 1
}

if [ "$#" -ne 5 ]; then
    usage
fi

# if [ "$#" -ne 5 ]; then
#     log "Usage: $0 <CONFIG_DIRECTORY> <NUM_THREADS> <REFERENCE_GTF> <STRAND_ST> <SPECIES>"
#     exit 1
# fi

config_directory=$1
num_threads=$2
reference_gtf=$3
strand_st=$4
species=$5

cd $config_directory

# INPUT & OUTPUT PATHS
filter_output_fpath="${config_directory}/3_filter_output"
counts_output_fpath="${config_directory}/4_stringtie_counts_output"
counts_log_fpath="${counts_output_fpath}/log"

mkdir -p ${counts_output_fpath} ${counts_log_fpath}

log "Processing stringtie counts for files in ${filter_output_fpath}"
log " " 
for filtered_file in "${filter_output_fpath}"/*.filt.bam; do
    if [[ -f "$filtered_file" ]]; then
        temp_file="${filtered_file##*/}" 
        libr_name="${temp_file%.filt.bam}"
        log "${libr_name}"
        log "Running stringtie: ${filtered_file}"
        log "( stringtie -p ${num_threads} -e -v -G ${reference_gtf} ${strand_st} ${filtered_file} -o ${counts_output_fpath}/${libr_name}.assemb.gtf -A ${counts_output_fpath}/${libr_name}.tab ) 2> ${counts_log_fpath}/${libr_name}_stringtie_step.log"
        ( stringtie -p ${num_threads} -e -v -G ${reference_gtf} ${strand_st} ${filtered_file} -o ${counts_output_fpath}/${libr_name}.assemb.gtf -A ${counts_output_fpath}/${libr_name}.tab ) 2> ${counts_log_fpath}/${libr_name}_stringtie_step.log
        # https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
        # -p <int>	Specify the number of processing threads (CPUs) to use for transcript assembly. The default is 1.
    else
        log "No ${filtered_file}.filt.bam file found in $filter_output_fpath"
    fi

done


# merging the stringtie outputs
log "Merging stringtie outputs..."
mkdir -p "${config_directory}/stringtie"
scp -pr "${counts_output_fpath}"/*.tab "${config_directory}/stringtie"

ls ./stringtie/ > sample_list
FILE_NO=`cat sample_list | wc -l`
COUNTER=1
STATUS=0
START_TIME=$(date +%s)

x=`cat sample_list | head -n $COUNTER | tail -n 1`

while read x
do
    cat stringtie/$x | awk '{ print $1"#"$2"#"$3"#"$4"#"$5"#"$6,$9;}' | sort -k1,1 > ./stringtie/$x.txt
done < sample_list

counter=0

while read x
do
    counter=$(($counter + 1))
    if [ $counter -eq 1 ]
    then
        cat ./stringtie/$x.txt > ./stringtie/tmp
    else
        join ./stringtie/tmp ./stringtie/$x.txt > ./stringtie/tmp2
        mv ./stringtie/tmp2 ./stringtie/tmp
    fi
done < sample_list

echo ensmbl_gene_id gene_name chr strand start stop `cat sample_list | tr '\n' " "` > ./stringtie/genes.tpm.txt
cat ./stringtie/tmp | tr "#" " " > ./stringtie/tmp2
cat ./stringtie/tmp2 | sort -k1,1 > ./stringtie/tmp
cat ./stringtie/tmp >> ./stringtie/genes.tpm.txt

log "Done merging output!"

FINISH_TIME=$(($(date +%s) - $START_TIME))
log "$FILE_NO stringtie jobs finished in $FINISH_TIME seconds!"

scp -pr "${config_directory}/stringtie/genes.tpm.txt" "${counts_output_fpath}"
rm -rf  "${config_directory}/stringtie" "${config_directory}/sample_list" "${counts_output_fpath}/tmp_*"

cd $counts_output_fpath

if [[ "$species" =~ ^[Hh][Uu][Mm][Aa][Nn]$ ]]; then
    ref_gene_id_name_human_file="${config_directory}/ref_human_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_gene_id_name_human_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
elif [[ "$species" =~ ^[Mm][Oo][Uu][Ss][Ee]$ ]]; then
    ref_gene_id_name_mouse_file="${config_directory}/ref_mouse_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_gene_id_name_mouse_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
elif [[ "$species" =~ ^[Rr][Aa][Tt][Tt][Uu][Ss]$ ]]; then
    ref_gene_id_name_rattus_file="${config_directory}/ref_rattus_geneid_genename_genebiotype.tsv"
    awk 'NR==FNR {gene_biotype[$1]=$3; next} FNR==1 {print $0, "gene_biotype"; next} ($1 in gene_biotype) {print $0, gene_biotype[$1]}' $ref_gene_id_name_rattus_file "${counts_output_fpath}/genes.tpm.txt" > "${counts_output_fpath}/genes_tpm.txt"
fi

log "Adding gene biotype for $species using $ref_file"
log "Processing completed."

exit