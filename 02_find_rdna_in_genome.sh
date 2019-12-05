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

predict_rdna() {
    barrnap --kingdom 'euk' --quiet --threads ${threads} --evalue 1e-06 --reject 0.2 -
}

#============================================================================================================
# PIPELINE START

# Export variables and check them
config_file="$1"
if ! [ -f "$config_file" ]
then
    echo_with_timestamp "ERROR: Config file not found or not a file. Please provide config file path as an argument: bash 02_find_rdna_in_genome.sh path/to/config.sh. Exiting..."
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
check_executable "blastn" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
check_executable "makeblastdb" "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
check_executable "barrnap" "https://github.com/tseemann/barrnap"
check_executable "bedtools" "https://github.com/arq5x/bedtools2/releases"
check_executable "nhmmer" "http://hmmer.org"

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

### Check if prefix is set and concatenate directory path and prefix
if [ -v "prefix" ]
then
    output_prefix="${out_dir}/${prefix}_"
    tmp_prefix="${tmp_dir}/${prefix}_"
else
    output_prefix="${out_dir}/rdna_pipeline_"
    tmp_prefix="${tmp_dir}/rdna_pipeline"
fi

### Check genome assembly
check_var "assembly_path"
check_path "${assembly_path}"
check_if_readable "${assembly_path}"
if ! [ -s "${assembly_path}" ]
then
    echo_with_timestamp "ERROR: ${assembly_path} has zero size. Exiting..."
    exit 1
fi
assembly_basename=`basename "${assembly_path}" | rev | cut -d. -f1 --complement | rev | tr '.' '_'`

### Check rDNA operon fasta
check_var "rdna_operon_path"
check_path "${rdna_operon_path}"
check_if_readable "${rdna_operon_path}"
if ! [ -s "${rdna_operon_path}" ]
then
    echo_with_timestamp "ERROR: ${rdna_operon_path} has zero size. Exiting..."
    exit 1
fi

echo_with_timestamp "INFO: Preflight check success, starting rDNA prediction in genome assembly"

# Predict rDNA in assembly
assembly_gff_path="${output_prefix}${assembly_basename}_rdna_prediction.gff"
if [ -s "${assembly_gff_path}" ]
then
    echo_with_timestamp "INFO: Found ${assembly_gff_path}"
else
    assembly_extention=`echo "${assembly_path}" | rev | cut -d. -f1 | rev`
    case ${assembly_extention} in
        "gz"|"GZ")
            assembly_opener="zcat"
            assembly_extention=`echo "${assembly_path}" | rev | cut -d. -f2 | rev`;;
        *)
            assembly_opener="cat";;
    esac
    case ${assembly_extention} in
        "fa"|"fasta"|"FA"|"FASTA"|"fna"|"FNA")
            echo_with_timestamp "INFO: Predicting rDNA in ${assembly_path}..."
            "${assembly_opener}" "${assembly_path}" | predict_rdna > "${assembly_gff_path}" || { echo_with_timestamp "ERROR: Nonzero exit status while predicting rDNA in ${assembly_path}. Exiting..." ; rm "${assembly_gff_path}"; exit 1; }
            echo_with_timestamp "INFO: rDNA predictions in genome stored to ${assembly_gff_path}";;
        *)
            echo_with_timestamp "ERROR: Unknown genome assembly format (${assembly_path}). Allowed formats: fasta, fa, fna (not case-sensitive, can be compressed with gzip). Exiting...";;
    esac
fi

# Map rDNA operon to gneome assembly
echo_with_timestamp "INFO: Mapping rDNA operon to genome assembly"
blast_db_prefix="${out_dir}/blast_db/${assembly_basename}_blast_db"
if [ -s "${blast_db_prefix}.nsq" ]
then
    echo_with_timestamp "INFO: Using BLAST database ${blast_db_prefix}"
else
    assembly_extention=`echo "${assembly_path}" | rev | cut -d. -f1 | rev`
    case ${assembly_extention} in
        "gz"|"GZ")
            zcat "${assembly_path}" | makeblastdb -dbtype nucl -max_file_sz '4GB' -out "${blast_db_prefix}" -title "${assembly_basename}" || { echo_with_timestamp "ERROR: Nonzero exit status while building BLAST database. Exiting..." ; exit 1 ; } ;;
        *)
            makeblastdb -in "${assembly_path}" -dbtype nucl -max_file_sz '4GB' -out "${blast_db_prefix}" -title "${assembly_basename}" || { echo_with_timestamp "ERROR: Nonzero exit status while building BLAST database. Exiting..." ; exit 1 ; } ;;
    esac
fi

blast_out_path="${output_prefix}rdna_operon_to_${assembly_basename}.tab"
blastn -task megablast -query "${rdna_operon_path}" -db "${blast_db_prefix}" -outfmt "7 qseqid sseqid length qstart qend sstart send pident qcovhsp" -out "${blast_out_path}" -num_threads ${threads}
cat "${blast_out_path}"
echo_with_timestamp "INFO: Mapping results stored in ${blast_out_path}"
echo_with_timestamp "END"