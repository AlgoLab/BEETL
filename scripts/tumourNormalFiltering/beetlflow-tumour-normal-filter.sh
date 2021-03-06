#!/bin/bash
################################################################################
##
## Copyright (c) 2013 Illumina, Inc.
##
## This software is covered by the "Illumina Non-Commercial Use Software
## and Source Code License Agreement" and any user of this software or
## source file is bound by the terms therein (see accompanying file
## Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
##
################################################################################

abscommand="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
abspath=`dirname "$abscommand"`


# Command line parsing

if [ $# != 5 ]
then
    echo "Usage:
  $0 <tumour_paired_reads1.fastq.gz> <tumour_paired_reads2.fastq.gz> <normal_paired_reads1.fastq.gz> <normal_paired_reads2.fastq.gz> <output directory>

EXAMPLE:
  $0 \\
    Tumour/lane1_NoIndex_L001_R1_001.fastq.gz \\
    Tumour/lane1_NoIndex_L001_R2_001.fastq.gz \\
    Normal/lane1_NoIndex_L001_R1_001.fastq.gz \\
    Normal/lane1_NoIndex_L001_R2_001.fastq.gz \\
    outDir
"
    exit 1
fi

TUMOUR1=`readlink -m "$1"`
TUMOUR2=`readlink -m "$2"`
NORMAL1=`readlink -m "$3"`
NORMAL2=`readlink -m "$4"`
OUTPUT_DIR=$5

echo "Tumour1: ${TUMOUR1}"
echo "Tumour2: ${TUMOUR2}"
echo "Normal1: ${NORMAL1}"
echo "Normal2: ${NORMAL2}"
echo "OutputDir: ${OUTPUT_DIR}"


BEETL_CONVERT=${abspath}/beetl-convert
BEETL_BWT=${abspath}/beetl-bwt
BEETL_INDEX=${abspath}/beetl-index
BEETL_COMPARE=${abspath}/beetl-compare
BEETL_EXTEND=${abspath}/beetl-extend

DATA_DIR=${abspath}/data
INPUT1_FASTA=${DATA_DIR}/testBeetlCompare.dataset1.fasta
INPUT2_FASTA=${DATA_DIR}/testBeetlCompare.dataset2.fasta
COMPARE_OUT=compare.out
SORTED_BKPT_A=${COMPARE_OUT}.sortedBKPT_setA
SORTED_BKPT_B=${COMPARE_OUT}.sortedBKPT_setB
EXTEND_A_OUT=extend_setA.out
EXTEND_B_OUT=extend_setB.out


function printAndRunCOMMAND {
    echo "${COMMAND}"
    echo "${COMMAND}" >> commands.log
    eval "${COMMAND}"
    if [ $? != 0 ]
    then
        echo "Error detected."
        exit 1
    fi
}

# Output directory initialisation
if [ -e "${OUTPUT_DIR}" ]
then
    echo "Error: Output directory \"${OUTPUT_DIR}\" already exists. Aborting to prevent overwriting."
    exit 2
fi

mkdir "${OUTPUT_DIR}"
if [ ! -e "${OUTPUT_DIR}" ]
then
    echo "Error: Cannot create output directory \"${OUTPUT_DIR}\". Aborting."
    exit 3
fi

cd ${OUTPUT_DIR}


# BWT1 creation
COMMAND="( zcat ${TUMOUR1} ${TUMOUR2} || cat ${TUMOUR1} ${TUMOUR2} ) | ${BEETL_BWT} -i - --input-format=fastq -o bwt1 --generate-end-pos-file --add-rev-comp --paired-reads-input=all1all2"
printAndRunCOMMAND


# BWT2 creation
COMMAND="( zcat ${NORMAL1} ${NORMAL2} || cat ${NORMAL1} ${NORMAL2} ) | ${BEETL_BWT} -i - --input-format=fastq -o bwt2 --generate-end-pos-file --add-rev-comp --paired-reads-input=all1all2"
printAndRunCOMMAND


# BWT1 indexing
COMMAND="${BEETL_INDEX} -i bwt1"
printAndRunCOMMAND


# BWT2 indexing
COMMAND="${BEETL_INDEX} -i bwt2"
printAndRunCOMMAND


# Comparison
COMMAND="${BEETL_COMPARE} -a bwt1 -b bwt2 -m tumour-normal -n 4 > ${COMPARE_OUT}"
printAndRunCOMMAND


# BKPT extraction for bwt1
COMMAND="grep BKPT ${COMPARE_OUT} | awk 'BEGIN { OFS=\"\t\"; FS=\" \" } {x+=\$7}   {if (\$7 < x/NR*3 && length(\$2)>12) print \$2, \$5, \$7}' > ${SORTED_BKPT_A}"
printAndRunCOMMAND


# BKPT extraction for bwt2
COMMAND="grep BKPT ${COMPARE_OUT} | awk 'BEGIN { OFS=\"\t\"; FS=\" \" } {x+=\$8}   {if (\$8 < x/NR*3 && length(\$2)>12) print \$2, \$6, \$8}' > ${SORTED_BKPT_B}"
printAndRunCOMMAND


# Extension from BWT1 position to read numbers
COMMAND="${BEETL_EXTEND} -i ${SORTED_BKPT_A} -b bwt1 -o ${EXTEND_A_OUT} --use-indexing"
printAndRunCOMMAND


# Extension from BWT2 position to read numbers
COMMAND="${BEETL_EXTEND} -i ${SORTED_BKPT_B} -b bwt2 -o ${EXTEND_B_OUT} --use-indexing"
printAndRunCOMMAND


# Extraction of reads for: tumour read 1
COMMAND="( zcat ${TUMOUR1} || cat ${TUMOUR1} ) | ${BEETL_CONVERT} -i - --input-format=fastq -o tumour_read1.fastq --extract-sequences=${EXTEND_A_OUT}"
printAndRunCOMMAND


# Extraction of reads for: tumour read 2
COMMAND="( zcat ${TUMOUR2} || cat ${TUMOUR2} ) | ${BEETL_CONVERT} -i - --input-format=fastq -o tumour_read2.fastq --extract-sequences=${EXTEND_A_OUT}"
printAndRunCOMMAND


# Extraction of reads for: normal read 1
COMMAND="( zcat ${NORMAL1} || cat ${NORMAL1} ) | ${BEETL_CONVERT} -i - --input-format=fastq -o normal_read1.fastq --extract-sequences=${EXTEND_A_OUT}"
printAndRunCOMMAND

# Extraction of reads for: normal read 2
COMMAND="( zcat ${NORMAL2} || cat ${NORMAL2} ) | ${BEETL_CONVERT} -i - --input-format=fastq -o normal_read2.fastq --extract-sequences=${EXTEND_A_OUT}"
printAndRunCOMMAND


# Cleanup
#rm -f bwt[12]-* compare.out*

echo ""
echo "Filtering completed"
