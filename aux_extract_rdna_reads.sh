#!/bin/bash

store_selected_rdna_reads_to_fasta() {
    local main_script_dir="$1"
    local selected_rdna_reads_id_path="$2"
    local pacbio_reads_path="$3"
    local output_fasta_path="$4"
python - <<EOF
import os
os.chdir("${main_script_dir}")
import rdna_reads_utils as utils
utils.extract_rdna_reads("${pacbio_reads_path}", "${selected_rdna_reads_id_path}", "${output_fasta_path}")
EOF
}


tmp_prefix="$1"
reads_path="$2"
main_script_dir="$3"
selected_rdna_reads_id_path="$4"
out_fasta_path="${tmp_prefix}`basename ${reads_path} | cut -f1 -d.`_selected_rdna_reads.fasta"

if [ -f "${out_fasta_path}" ]
then
    continue
fi
out_fasta_tmp_path="${tmp_prefix}`basename ${reads_path} | cut -f1 -d.`_selected_rdna_reads.tmp"
store_selected_rdna_reads_to_fasta "${main_script_dir}" "${selected_rdna_reads_id_path}" "${reads_path}" "${out_fasta_tmp_path}" && mv "${out_fasta_tmp_path}" "${out_fasta_path}"