#!/usr/bin/env bash
################################################################################
##
## THE MIT LICENSE
##
## Copyright 2023 Jacob Morrison <jacob.morrison@vai.org>
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the “Software”), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
## DESCRIPTION
##     Script to process STORM data through STARsolo
##
## NOTES
##     - synthbar and STAR (version 2.7.9+) must be in PATH
##
## CREATED BY
##     Jacob Morrison
##
## CREATION DATE
##     October 2022
##
## UPDATE NOTES
##     - Oct 2022
##         - Initial creation
##     - May 2023
##         - MIT License added
##
################################################################################

set -euo pipefail

# Check for synthbar and STAR in PATH
function check_path {
    # Check for synthbar
    if [[ `which synthbar 2>&1 > /dev/null` ]]; then
        >&2 echo "synthbar does not exit in PATH"
        exit 1
    else
        >&2 echo "Using synthbar found at: `which synthbar`"
    fi

    # Check for STAR and it's version
    if [[ `which STAR 2>&1 > /dev/null` ]]; then
        >&2 echo "STAR does not exit in PATH"
        exit 1
    else
        CHECK=`STAR --version | \
            awk '{
                split($0,a,".")
                if (a[2] > 7) {
                    print "yes"
                } else if (a[2] == 7) {
                    match(a[3], /([[:digit:]]+)([[:alpha:]]+)/, b)
                    if (b[1] >= 9) {
                        print "yes"
                    } else {
                        print "no"
                    }
                } else {
                    print "no"
                }
            }'`
        if [[ "${CHECK}" == "yes" ]]; then
            >&2 echo -e "Using STAR found at: `which STAR`\n"
        else
            >&2 echo "Version of STAR (`STAR --version`) too old. Must be >= 2.7.9a"
        fi
    fi
}

function pipeline {
    # Cursory check of read 1 file name to see if FASTQ is gzipped
    # If not, then gzip it so we're all set to go for STAR running
    if [[ ! ${read1} =~ ".gz" ]]; then
        >&2 echo "${read1} is not gzipped. Compressing now."
        gzip ${read}

        bc_read1=${read1}.gz
    else 
        bc_read1=${read1}
    fi

    # Check if read 2 is gzipped (based on file name) and if it's .fastq or .fq
    if [[ "${read2}" =~ ".gz" ]]; then
        is_gz=".gz"
    else
        is_gz=
    fi

    if [[ "${read2}" =~ ".fastq" ]]; then
        fastq_name=".fastq"
    elif [[ "${read2}" =~ ".fq" ]]; then
        fastq_name=".fq"
    else
        >&2 echo "Unknown FASTQ extension: ${read2}"
        exit 1
    fi

    bc_read2=${read2/${fastq_name}${is_gz}/_barcoded${fastq_name}.gz}

    if [[ -f ${bc_read2} ]]; then
        >&2 echo "${bc_read2} exists. Skipping synthbar run since it seems like barcodes have already been added."
    else
        synthbar ${read2} | gzip > ${bc_read2}
    fi

    STAR \
        --runThreadN "${n_threads}" \
        --genomeDir "${index}" \
        --outFileNamePrefix "${o_path}" \
        --readFilesIn ${bc_read1} ${bc_read2} \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
        --readFilesCommand zcat \
        --soloStrand Reverse \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist None \
        --soloBarcodeMate 2 \
        --clip5pNbases 0 21 \
        --soloCBstart 1 \
        --soloCBlen 7 \
        --soloUMIstart 8 \
        --soloUMIlen 8 \
        --soloBarcodeReadLength 0 \
        --soloUMIdedup Exact \
	--soloMultiMappers EM \
        --soloFeatures GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx

    if [[ ! "${keep_bc}" == true ]]; then
        rm -f ${bc_read2}
    fi
}

function version {
    >&2 echo -e "\nProgram: STORMsolo.sh"
    >&2 echo -e "Creator: Jacob Morrison <jacob.morrison@vai.org>"
    >&2 echo -e "Version: 1.0.0\n"
}

function usage {
    version

    >&2 echo -e "\nUsage: STORMsolo.sh [options] <STAR index> <Read 1 FASTQ> <Read 2 FASTQ>\n"
    >&2 echo -e "Required Inputs:"
    >&2 echo -e "\tSTAR index          : Path to directory with STAR index (input to --genomeDir)"
    >&2 echo -e "\tRead 1 FASTQ        : FASTQ for read 1"
    >&2 echo -e "\tRead 2 FASTQ        : FASTQ for read 2 (will be run through synthbar to add synthetic barcode)"
    >&2 echo -e "Optional Inputs:"
    >&2 echo -e "\t-t, --threads       : Number of threads to use when running STAR [default: 1]"
    >&2 echo -e "\t-k, --keep-barcoded : Keep read 2 FASTQ file with added barcodes [default: delete FASTQ]"
    >&2 echo -e "\t-o, --output-path   : Path to location for output files [default: ./]"
    >&2 echo -e "\t-h, --help          : usage and version info"
    >&2 echo -e "\t-v, --version       : version info"
    >&2 echo -e ""
}

# Actual processing
################################################################################

# Check for required programs before going on
check_path

# Initialize default values
n_threads=1
keep_bc=false
o_path="./"

# Process command line arguments
OPTS=$(getopt \
    --options t:ko:hv \
    --long threads:,keep-barcoded,output-path:,help,version \
    --name "$(basename "$0")" \
    -- "$@"
)
eval set -- ${OPTS}

while true; do
    case "$1" in
        -h|--help )
            usage
            exit 0
            ;;
        -v|--version )
            version
            exit 0
            ;;
        -t|--threads )
            n_threads="$2"
            shift 2
            ;;
        -k|--keep-barcoded )
            keep_bc=true
            shift
            ;;
        -o|--output-path )
            o_path="$2"
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure all inputs are there
if [[ $# -ne 3 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
index=$1
read1=$2
read2=$3

# Perform checks on inputs
if [[ ! -d "${index}" ]]; then
    >&2 echo "Index directory missing: ${index}"
    exit 1
fi

if [[ ! -f "${read1}" ]]; then
    >&2 echo "Read 1 FASTQ missing: ${index}"
    exit 1
fi

if [[ ! -f "${read2}" ]]; then
    >&2 echo "Read 2 FASTQ missing: ${index}"
    exit 1
fi

# Configuration
>&2 echo "## Running STORMsolo with the following configuration ##"
>&2 echo "##----------------------------------------------------##"
>&2 echo "INDEX                      : ${index}"
>&2 echo "READ 1 FASTQ               : ${read1}"
>&2 echo "READ 2 FASTQ               : ${read2}"
>&2 echo "Number of Threads          : ${n_threads}"
>&2 echo "Keep Barcoded Read 2 FASTQ : ${keep_bc}"
>&2 echo "Output Path                : ${o_path}"
>&2 echo "##----------------------------------------------------##"

pipeline
