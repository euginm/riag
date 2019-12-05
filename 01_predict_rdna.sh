#!/bin/bash
#============================================================================================================
# FUNCTIONS
# Define a timestamp function
echo_with_timestamp() {
    date +"[%F %T] $1"
}

check_var() {
    local var_name="${1}"
    if ! [ -v "${var_name}" ]
    then
        echo_with_timestamp "ERROR: ${var_name} not defined in config file. Exiting..."
        exit 1
    fi
}

# Define executable check function
check_executable() {
    local command_name="$1"
    local install_link="$2"
    if [ -v "${command_name}" ]
    then
        local command_bin="${!command_name}"
        if ! [ -f "$command_bin" ]
        then
            echo_with_timestamp "ERROR: ${command_bin} not found, please correct the path for ${command_name} in config file or install ( ${install_link} )'. Exiting..."
            exit 1
        elif ! [ -x "$command_bin" ]
        then
            echo_with_timestamp "ERROR: ${command_bin} is not excecutable, please use 'chmod +x ${command_bin}'. Exiting..."
            exit 1
        else
            local command_dir=`dirname "${!command_name}"`
            export PATH="${command_dir}":$PATH || { echo_with_timestamp "ERROR: cannot add ${command_dir} to PATH. Exiting..." ; exit 1; }
        fi
    elif ! [ -x "$(command -v $command_name)" ]
    then
        echo_with_timestamp "ERROR: ${command_name} not found, please install ( ${install_link} ) or provide path in config file. Exiting..."
        exit 1
    fi
}

check_dir() {
    local dir="$1"
    if ! [ -d "$dir" ]
    then
        echo_with_timestamp "INFO: ${dir} not found, creating..."
        mkdir -p "$dir" || { echo_with_timestamp "ERROR: Can't create ${dir}. Exiting..." ; exit 1; }
    fi
}

check_path() {
    local path="$1"
    if ! [ -e "$path" ]
    then
        echo_with_timestamp "ERROR: ${path} not found. Exiting..."
        exit 1
    fi
}

check_if_readable() {
    local path="$1"
    if ! [ -r "$path" ]
    then
        echo_with_timestamp "ERROR: ${path} requires read premission, please use 'chmod +r ${path}'. Exiting..."
        exit 1
    fi
}

check_if_writable() {
    local path="$1"
    if ! [ -w "$path" ]
    then
        echo_with_timestamp "ERROR: ${path} requires write premission, please use 'chmod +w ${path}'. Exiting..."
        exit 1
    fi
}

fastq_to_fasta() {
    paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n'
}

predict_rdna() {
    barrnap --kingdom 'euk' --quiet --threads ${threads} --evalue 1e-04 --reject 0.1 -
}

select_rdna_reads_for_assembly() {
    local main_script_dir="$1"
    local rdna_reads_gff_path="$2"
    local out_prefix="$3"
    local evalue="$4"
python - <<EOF
import os
os.chdir("${main_script_dir}")
import rdna_reads_utils as utils
gff_dict = utils.store_gff_as_dict("${rdna_reads_gff_path}", "${evalue}")
utils.pick_rdna_reads_to_assembly(gff_dict, "${out_prefix}")
EOF
}

assemble_rrna_gene_consensus_flye() {
    local selected_rdna_reads_path="$1"
    local output_dir="$2"
    local rdna_repeat_size="$3"
    local n_iterations=$4
    flye --${lr_type} "${selected_rdna_reads_path}" --genome-size "${rdna_repeat_size}" --out-dir "${output_dir}" \
      --threads "${threads}" --iterations ${n_iterations} --asm-coverage 100
}

#============================================================================================================
# PIPELINE START
main_script_dir=`dirname "$0"`
python_utils_path="${main_script_dir}/rdna_reads_utils.py"
if ! [ -f "${python_utils_path}" ]
then
    echo_with_timestamp "ERROR: ${python_utils_path} not found, please place rdna_reads_utils.py to the directory with 01_predict_rdna.sh script (${main_script_dir}). Exiting..."
    exit 1
fi

aux_extract_rdna_reads_script="${main_script_dir}/aux_extract_rdna_reads.sh"
if ! [ -f "${aux_extract_rdna_reads_script}" ]
then
    echo_with_timestamp "ERROR: ${aux_extract_rdna_reads_script} not found, please place aux_extract_rdna_reads.sh to the directory with 01_predict_rdna.sh script (${main_script_dir}). Exiting..."
    exit 1
fi

# Export variables and check them
config_file="$1"
if ! [ -f "$config_file" ]
then
    echo_with_timestamp "ERROR: Config file not found or not a file. Please provide config file path as an argument: bash 01_predict_rdna.sh path/to/config.sh. Exiting..."
    exit 1
elif ! [ -s "$config_file" ]
then
    echo_with_timestamp "ERROR: Config has zero size. Please correct config file or the path. Exiting..."
    exit 1
else
    check_if_readable "$config_file"
fi
echo_with_timestamp "INFO: Reading ${config_file}..."
source "${config_file}" || { echo_with_timestamp "ERROR: Can't source ${config_file}. Exiting..." ; exit 1; }

echo_with_timestamp "INFO: Preflight check..."

## Check executables and set them to PATH if needed
check_executable "barrnap" "https://github.com/tseemann/barrnap"
check_executable "samtools" "https://github.com/samtools/samtools"
check_executable "bedtools" "https://github.com/arq5x/bedtools2/releases"
check_executable "nhmmer" "http://hmmer.org"
check_executable "flye" "https://github.com/fenderglass/Flye"


## Check if required python modules installed
python -c "from Bio import SeqIO, Seq" || { echo_with_timestamp "ERROR: biopython library not available, please install (e.g., with 'pip install biopython --user'). Exiting..." ; exit 1; }

## Check general variables
### Check output directory
check_var "out_dir"
check_dir "${out_dir}"
check_if_readable "${out_dir}"
check_if_writable "${out_dir}"

### Check if tmp_dir variable is defined, use tmp dir from PATH otherwise
if [ -v "tmp_dir" ]
then
    echo_with_timestamp "INFO: Using ${tmp_dir} as tmp directory"
    check_dir "${tmp_dir}"
    check_if_readable "${tmp_dir}"
    check_if_writable "${tmp_dir}"
    export TMPDIR="${tmp_dir}"
else
    echo_with_timestamp "INFO: tmp_dir not defined, using ${TMPDIR} as tmp directory"
    tmp_dir=${TMPDIR}
fi

### Check long reads directory
check_var "lr_dir"
check_path "${lr_dir}"
check_if_readable "${lr_dir}"

if ! find "${lr_dir}" -mindepth 1 -print -quit 2>/dev/null | grep -q .;
then
    echo_with_timestamp "ERROR: ${lr_dir} is empty. Exiting..."
    exit 1
fi

### Check threads number
if ! [ -v "threads" ]
then
    threads=`nproc --all`
elif ! [ "${threads}" -ge 0 ]
then
    echo_with_timestamp "ERROR: threads variable should be a positive integer. Exiting..."
    exit 1
elif [ "${threads}" -eq 0 ]
then
    threads=`nproc --all`
    echo_with_timestamp "WARNING: threads value from config equals 0, using threads=${threads} instead"
fi

### Check e-value variable
if ! [ -v "evalue" ]
then
    evalue="1e-06"
fi

### Check rdna gene size variable
if ! [ -v "rdna_repeat_size" ]
then
    rdna_repeat_size="30k"
fi

### Check number of polishing iterations variable
if ! [ -v "n_polish_iterations" ]
then
    n_polish_iterations=3
elif ! [ "${n_polish_iterations}" -ge 1 ]
then
    echo_with_timestamp "ERROR: n_polish_iterations variable should be a positive integer. Exiting..."
    exit 1
fi

### Check reads type variable
if ! [ -v "lr_type" ]
then
    lr_type='pacbio-raw'
elif [ ! "${lr_type}" = 'pacbio-raw' ] && [ ! "${lr_type}" = 'nano-raw' ]
then
    echo_with_timestamp "ERROR: lr_type variable should be 'pacbio-raw' or 'nano-raw', current value: '${lr_type}'. Exiting..."
    exit 1
fi

### Check if prefix is set and concatenate directory path and prefix
if [ -v "prefix" ]
then
    output_prefix="${out_dir}/${prefix}_"
    tmp_prefix="${tmp_dir}/${prefix}_"
else
    output_prefix="${out_dir}/rdna_pipeline_"
    tmp_prefix="${tmp_dir}/rdna_pipeline"
fi

### Check long reads coverage-related variables
if [ -v "lr_coverage" ]
then
    if ! [ "${lr_coverage}" -gt 0 ]
    then
        echo_with_timestamp "ERROR: lr_coverage variable should be a positive integer. Exiting..."
        exit 1
    else
        echo_with_timestamp "INFO: using lr_coverage=${lr_coverage}"
    fi
elif ! [ -v "lr_to_assembly_bam_path" ]
then
    echo_with_timestamp "INFO: lr_coverage and lr_to_assembly_bam_path variables not provided, will map long reads to assembly to compute median long reads coverage"
    check_executable "minimap2" "https://github.com/lh3/minimap2"
    lr_to_assembly_bam_path="${output_prefix}lr_to_${assembly_basename}.bam"
    #### Check assembly
    check_var "assembly_path"
    check_path "${assembly_path}"
    check_if_readable "${assembly_path}"
    if ! [ -s "${assembly_path}" ]
    then
        echo_with_timestamp "ERROR: ${assembly_path} has zero size. Exiting..."
        exit 1
    fi
    assembly_basename=`basename "${assembly_path}" | cut -d. -f1`
else
    check_path "${lr_to_assembly_bam_path}"
    check_if_readable "${lr_to_assembly_bam_path}"
    echo_with_timestamp "INFO: lr_coverage variable not provided, using ${lr_to_assembly_bam_path} to compute median long reads coverage"
    if ! [ -v "assembly_path" ]
    then
        assembly_basename=`basename "${lr_to_assembly_bam_path}" | cut -d. -f1`
    else
        assembly_basename=`basename "${assembly_path}" | cut -d. -f1`
    fi
fi

echo_with_timestamp "INFO: Preflight check success, starting rDNA prediction in long reads and rDNA copy number estimation"

# Predict rDNA in long reads
rdna_reads_gff_path="${output_prefix}rdna_reads.gff"
if [ -s "${rdna_reads_gff_path}" ]
then
    echo_with_timestamp "INFO: Found ${rdna_reads_gff_path}, starting rDNA reads assessment and selection"
else
    for reads_path in "${lr_dir}"/*.*
    do
        reads_basename=`basename "${reads_path}" | cut -d. -f1`
        tmp_gff_path="${tmp_prefix}${reads_basename}_rdna_prediction.gff"
        if [ -s "${tmp_gff_path}" ]
        then
            echo_with_timestamp "INFO: Found ${tmp_gff_path} with `grep 'barrnap' ${tmp_gff_path} | cut -f1 | sort | uniq | wc -l` rDNA reads"
            continue
        fi
        reads_extention=`echo "${reads_path}" | rev | cut -d. -f1 | rev`
        case ${reads_extention} in
            "gz"|"GZ")
                reads_opener="zcat"
                reads_extention=`echo "${reads_path}" | rev | cut -d. -f2 | rev`;;
            *)
                reads_opener="cat";;
        esac
        case ${reads_extention} in
            "fa"|"fasta"|"FA"|"FASTA"|"fna"|"FNA")
            echo_with_timestamp "INFO: Predicting rDNA in ${reads_path}..."
            "${reads_opener}" "${reads_path}" | predict_rdna > "${tmp_gff_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while predicting rDNA in ${reads_path}. Exiting..." ; rm "${tmp_gff_path}"; exit 1; }
            echo_with_timestamp "INFO: Identified `grep 'barrnap' ${tmp_gff_path} | cut -f1 | sort | uniq | wc -l` reads with rDNA features.";;
        "fq"|"fastq"|"FQ"|"FASTQ")
            echo_with_timestamp "INFO: Predicting rDNA in ${reads_path}..."
            "${reads_opener}" "${reads_path}" | fastq_to_fasta | predict_rdna > "${tmp_gff_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while predicting rDNA in ${reads_path}. Exiting..." ; rm "${tmp_gff_path}"; exit 1; }
            echo_with_timestamp "INFO: Identified `grep 'barrnap' ${tmp_gff_path} | cut -f1 | sort | uniq | wc -l` reads with rDNA features";;
        *)
            echo_with_timestamp "WARNING: long reads file ${reads_path} has unknown format. Allowed formats: fasta, fa, fna, fastq, fq (not case-sensitive, can be compressed with gzip). Skipping file...";;
        esac
    done

    ## Concatenate to one file
    echo "##gff-version 3" > "${rdna_reads_gff_path}"
    for tmp_gff_path in "${tmp_prefix}"*_rdna_prediction.gff
    do
        sed 1d "${tmp_gff_path}" >> "${rdna_reads_gff_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while concatenating rDNA reads annotations to ${rdna_reads_gff_path}. Exiting..." ; rm "${rdna_reads_gff_path}"; exit 1; }
    done
    for tmp_gff_path in "${tmp_prefix}"*_rdna_prediction.gff
    do
        rm "${tmp_gff_path}"
    done
    echo_with_timestamp "INFO: rDNA reads annotation stored to ${rdna_reads_gff_path}, starting rDNA reads assessment and selection"
fi

# Run rDNA reads qualitative assessment and selection
rdna_reads_stats_path="${output_prefix}rdna_reads.stats"
selected_rdna_reads_id_path="${output_prefix}rdna_reads_for_assembly.list"
selected_rdna_reads_path="${output_prefix}_rdna_reads_for_assembly.fasta"
if [ ! -s "${rdna_reads_stats_path}" ] || [ ! -s "${selected_rdna_reads_id_path}" ]
then
    rdna_out_prefix="${output_prefix}rdna_reads"
    echo_with_timestamp "INFO: Running rDNA reads qualitative assessment and selecting rDNA reads for assembly"
    select_rdna_reads_for_assembly "${main_script_dir}" "${rdna_reads_gff_path}" "${rdna_out_prefix}" "${evalue}"
    cat "${rdna_reads_stats_path}" || { echo_with_timestamp "ERROR: can't open ${rdna_reads_stats_path}. Exiting..." ; rm "${rdna_reads_gff_path}"; exit 1; }
    echo_with_timestamp "INFO: `cat ${selected_rdna_reads_id_path} | sort | uniq | wc -l` long reads selected for rRNA gene assembly" || { echo_with_timestamp "ERROR: can't open ${selected_rdna_reads_id_path}. Exiting..." ; rm "${rdna_reads_gff_path}"; exit 1; }
    echo_with_timestamp "INFO: rDNA reads statistics saved to ${rdna_reads_stats_path}"
    echo_with_timestamp "INFO: Selected long rDNA read ids saved to ${selected_rdna_reads_id_path}"
    echo_with_timestamp "INFO: All long rDNA read ids (including choosed for assemly) grouped by rdna composition saved to ${rdna_out_prefix}.json"
fi

# Store selected rDNA reads to file
n_selected_rdna_reads_in_list=`cat "${selected_rdna_reads_id_path}" | sort | uniq | wc -l`
if ! [ -s "${selected_rdna_reads_path}" ]
then
    echo_with_timestamp "INFO: Extracting ${n_selected_rdna_reads_in_list} selected long reads to ${tmp_dir} ..."
    parallel -j ${threads} --halt now,fail=1 --linebuffer "{1} {2} {3} {4} {5}" ::: "${aux_extract_rdna_reads_script}" ::: "${tmp_prefix}" ::: "${lr_dir}"/*.* ::: "${main_script_dir}" ::: "${selected_rdna_reads_id_path}" || exit 1
    echo_with_timestamp "INFO: Storing ${n_selected_rdna_reads_in_list} selected long reads to ${selected_rdna_reads_path}..."
    touch "${selected_rdna_reads_path}"
    for tmp_rdna_reads in "${tmp_prefix}"*_selected_rdna_reads.fasta
    do
        cat "${tmp_rdna_reads}" >> "${selected_rdna_reads_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while concatenating selected rDNA reads to ${selected_rdna_reads_path}. Exiting..." ; rm "${selected_rdna_reads_path}"; exit 1; }
    done
    for tmp_rdna_reads in "${tmp_prefix}"*_selected_rdna_reads.fasta
    do
        rm "${tmp_rdna_reads}"
    done
fi

n_selected_rdna_reads_in_fasta=`grep '>' "${selected_rdna_reads_path}" | wc -l`
if ! [ ${n_selected_rdna_reads_in_fasta} -eq ${n_selected_rdna_reads_in_list} ]
then
    echo_with_timestamp "WARNING: Found ${n_selected_rdna_reads_in_fasta} reads in ${selected_rdna_reads_path}, while there are ${n_selected_rdna_reads_in_list} read ids in ${selected_rdna_reads_id_path}. It may indicate pipeline errors."
else
    echo_with_timestamp "INFO: Found ${n_selected_rdna_reads_in_fasta} reads in ${selected_rdna_reads_path}, starting rRNA gene assembly."
fi

# Assemble rRNA gene consensus with flye
flye_out_dir="${output_prefix}flye_assembly"
flye_assembly_link="${output_prefix}rdna_flye_assembly.fasta"
if [ -L "${flye_assembly_link}" ]
then
    echo_with_timestamp "INFO: Found flye assembly ${flye_assembly_link}."
else
    assemble_rrna_gene_consensus_flye "${selected_rdna_reads_path}" "${flye_out_dir}" "${rdna_repeat_size}" ${n_polish_iterations} || { echo_with_timestamp "ERROR: Nonzero exit status while rRNA gene consensus assembly with flye. Exiting..." ; rm -rf "${flye_out_dir}"; exit 1; }
    ln -s "${output_prefix}flye_assembly/assembly.fasta" "${flye_assembly_link}"
    samtools faidx "${flye_assembly_link}"
    echo_with_timestamp "INFO: Predicting rDNA in ${flye_assembly_link}..."
    cat "${flye_assembly_link}" | predict_rdna > "${flye_assembly_gff_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while predicting rDNA in ${flye_assembly_gff_path}. Exiting..." ; rm "${flye_assembly_gff_path}"; exit 1; }
    echo_with_timestamp "INFO: Predicted rDNA subunits in flye assembly stored to ${flye_assembly_gff_path}."
fi

# Copy number estimation
rdna_lr_coverage=`grep 'Total NOR long reads:' "${rdna_reads_stats_path}" | rev | cut -d' ' -f1 | rev`
if ! [ -v "lr_coverage" ]
then
    if ! [ -s "${lr_to_assembly_bam_path}" ]
    then
        echo_with_timestamp "INFO: Long read alignments to assembly not found, creating..."
        for reads_path in "${lr_dir}"/*.*
        do
            reads_basename=`basename "${reads_path}" | cut -d. -f1`
            alignments_tmp_prefix="${tmp_prefix}${reads_basename}_to_${assembly_basename}"
            if [ -s "${alignments_tmp_prefix}.bam" ]
            then
                echo_with_timestamp "INFO: Found ${alignments_tmp_prefix}.bam"
                continue
            fi
            reads_extention=`echo "${reads_path}" | rev | cut -d. -f1 | rev`
            case ${reads_extention} in
                "gz"|"GZ")
                    reads_extention=`echo "${reads_path}" | rev | cut -d. -f2 | rev`;;
            esac
            case ${reads_extention} in
                "fa"|"fasta"|"FA"|"FASTA"|"fna"|"FNA"|"fq"|"fastq"|"FQ"|"FASTQ")
                    echo_with_timestamp "INFO: Mapping ${reads_path} to ${assembly_path}...";;
                *)
                    echo_with_timestamp "WARNING: long reads file ${reads_path} has unknown format. Allowed formats: fasta, fa, fna, fastq, fq (not case-sensitive, can be compressed with gzip). Skipping file..."
                    continue;;
            esac
            if ! [ -s "${alignments_tmp_prefix}".sam ]
            then
                minimap2 -x map-pb -a -t ${threads} -2 "${assembly_path}" "${reads_path}" > "${alignments_tmp_prefix}".sam.tmp || { echo_with_timestamp "ERROR: Nonzero exit status while mapping ${reads_path} to ${assembly_path}. Exiting..." ; rm "${alignments_tmp_prefix}".sam.tmp; exit 1; }
                mv "${alignments_tmp_prefix}".sam.tmp "${alignments_tmp_prefix}".sam
            fi
            samtools view -F 256 -@ ${threads} -b "${alignments_tmp_prefix}".sam > "${alignments_tmp_prefix}".bam || { echo_with_timestamp "ERROR: Nonzero exit status while converting ${alignments_tmp_prefix}.sam to ${alignments_tmp_prefix}.bam. Exiting..." ; rm "${alignments_tmp_prefix}".bam; exit 1; }
            rm "${alignments_tmp_prefix}".sam
        done
        echo_with_timestamp "INFO: Merging multiple bam files to one..."
        samtools merge -@ ${threads} "${lr_to_assembly_bam_path}" "${tmp_prefix}"*_to_${assembly_basename}.bam || { echo_with_timestamp "ERROR: Nonzero exit status while merging temporary bam files to ${lr_to_assembly_bam_path}. Exiting..." ; rm "${lr_to_assembly_bam_path}"; exit 1; }
        echo_with_timestamp "INFO: Sorting final bam..."
        samtools sort -@ ${threads} -T "${tmp_prefix}all_lr_to_${assembly_basename}" -o "${lr_to_assembly_bam_path}".sorted "${lr_to_assembly_bam_path}"
        mv "${lr_to_assembly_bam_path}".sorted "${lr_to_assembly_bam_path}"
        echo_with_timestamp "INFO: Indexing final bam..."
        samtools index "${lr_to_assembly_bam_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while indexing ${lr_to_assembly_bam_path}. Exiting..." ; exit 1; }
        for tmp_bam_file in "${tmp_prefix}"*_to_${assembly_basename}.bam
        do
            rm "${tmp_bam_file}"
        done
        echo_with_timestamp "INFO: Done! long read alignments to assembly stored in ${lr_to_assembly_bam_path}"
    fi
    longest_scaffold=`samtools view -H "${lr_to_assembly_bam_path}" | grep 'SN:' | sort -k 3 -t':' -g -r | head -1 | cut -f2 | cut -d':' -f2` || { echo_with_timestamp "ERROR: Nonzero exit status while trying to infer the longest scaffold in ${lr_to_assembly_bam_path}. Exiting..." ; exit 1; }
    longest_scaffold_lr_median_coverage_file="${output_prefix}${assembly_basename}_${longest_scaffold}_median_coverage.txt"
    if ! [ -s "${longest_scaffold_lr_median_coverage_file}" ]
    then
        echo_with_timestamp "INFO: Computing ${longest_scaffold} median long reads coverage..."
        samtools depth -r "${longest_scaffold}" "${lr_to_assembly_bam_path}" | cut -f3 | sort -n | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' >> "${longest_scaffold_lr_median_coverage_file}" || { echo_with_timestamp "ERROR: Nonzero exit status while trying to compute ${longest_scaffold} median long reads coverage in ${lr_to_assembly_bam_path}. Exiting..." ; exit 1; }
        echo_with_timestamp "INFO: Stored median long reads coverage value to ${longest_scaffold_lr_median_coverage_file}"
    fi
    lr_coverage=`cat "${longest_scaffold_lr_median_coverage_file}"`
    echo_with_timestamp "INFO: Median long reads coverage of ${longest_scaffold} is ${lr_coverage}"
fi
rdna_copy_number=`echo ${rdna_lr_coverage} ${lr_coverage} | awk '{print int($1/$2 + 0.5)}'`
echo_with_timestamp "INFO: Approximate number of rDNA copies in genome is ${rdna_copy_number}. Note that this value is an underestimation and therefore a lower bound"
echo -e "Approximate* number of rDNA copies in genome: ${rdna_copy_number} \n\t*based on conservative estimation of rDNA long reads number (${rdna_lr_coverage}) and long reads coverage (${lr_coverage})" >> "${rdna_reads_stats_path}"
echo_with_timestamp "INFO: ${rdna_reads_stats_path} updated with approximate rDNA copy number value"


# End
echo_with_timestamp "END OF RNDA PREDICTION PIPELINE: Please review the predicted rDNA reads properties, the assembly and the rDNA prediction in the assembly. If the results are not satisfying, try other methods: (A) Pick other rDNA reads for assembly, store their id's to ${selected_rdna_reads_id_path}, remove ${selected_rdna_reads_path}, ${flye_assembly_link}, and run the script again; (B) Use self-correction methods on selected reads (e.g., canu, FLAS, CONSENT, LoRMA), pick a corrected read with highest overlap and cut it so it represents an rDNA operon (optional, wit intergenic spacer). After that, update config file with rDNA operon sequence and run next script (02_find_rdna_in_genome.sh)"