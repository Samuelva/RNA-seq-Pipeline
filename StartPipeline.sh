#!/bin/bash

# Prevent script from overwriting files.
set -o noclobber

# General parameters
READS="" # Directory with the reads.
REFERENCE="" # Reference genome file.
GFF="" # Genome annotation file (should be gff 3.0 from NCBI).
CONDITIONS="" # File with specified conditions for every read file.
METHOD=""
SPECIES="NULL"
# Directory to save pipeline output to. If no output directory is given as
# argument, the pipeline will make its own with the current date and time as name.
OUTPUT="Output_$(date +%F-%T)"
EXTENSION="fq" # Extension 
THREADS=1 # Threads to be used for programs.

# Trimmomatic parameters
ILLUMINACLIP="TruSeq3-SE.fa:2:30:10"
LEADING="3"
TRAILING="3"
SLIDINGWINDOW="4:15"
MINLEN="36"

# Locations to where files will be saved

if [ ${#1} == 0 ]; then
    # Checks if any argument is given. If not
    # the pipeline is aborted.
    echo "Please provide some input..."
    exit 1
fi

# Doesn't work without the double brackets.
# Checks if required parameters are given.
# -------------------------------------------------------------------
# Can currently be easily fooled if a given file or directory contains something
# like -f, -r, -g or -c.
# -------------------------------------------------------------------
if [[ ! $@ == *"-f"* || ! $@ == *"-r"* || ! $@ == *"-g"* || ! $@ == *"-c"* || ! $@ == *"-m"* ]]; then
    echo "Parameters -f, -r, -g, -c, or -m missing..."
    echo "Please see the userguide."
    exit 1
fi

checkArg () {
    # Checks whether an argument is given with a parameter.
    # Aborts the pipeline if not.
    if [ $# -eq 1 ]; then
        echo "No argument with parameter $1 given..."
        exit 1
    fi
}

checkFile () {
    # Checks whether the given file exists.
    # Aborts the pipeline if not.
    if [ ! -f $1 ]; then
        echo "File $1 does not exist..."
        exit 1
    fi
}

checkDirectory () {
    # Checks whether the given directory exists and if it is empty.
    # If it doesnt exist, or the directory is empty, the pipeline
    # will be aborted.
    if [ ! -d $1 ]; then
        echo "Directory $1 does not exist..."
        exit 1
    fi
    if [ -z "$(ls -A $1)" ]; then
        echo "Reads directory is empty..."
        exit 1
    fi
}

# checkMethod () {
#     if [[ "${1,,}" != "glmlrt" || "${1,,}" != "glmqlf" || "${1,,}" != "et" || "${1,,}" != "exacttest"]]; then
#         echo "No valid method specified..."
#         exit 1
#     fi
# }

# Parameter assigning
while [[ $# -gt 0 ]]; do
    arg="$1"
    # Check whether an argument is given with a parameter
    if [[ $(echo $2 | grep ^-) ]]; then
        echo "No argument with parameter $1 given..."
        exit 1
    fi
    # checkArg $arg $2
    case $arg in
        # General settings
        -f|--reads) checkDirectory $2; READS="$2"; shift;;
        -r|--reference) checkFile $2; REFERENCE="$2"; shift;;
        -g|--gff) checkFile $2; GFF="$2"; shift;;
        -c|--conditions) checkFile $2; CONDITIONS="$2"; shift;;
        -m|--method) METHOD="$2"; shift;;
        -s|--species) SPECIES="$2"; shift;;
        -o|--output) OUTPUT="$2"; shift;;
        -e|--extension) EXTENSION="$2"; shift;;
        -t|--threads) THREADS="$2"; shift;;
        -h|--help) help; shift;;
        # Trimmomatic settings
        -ti|--illuminaclip) ILLUMINACLIP="$2"; shift;;
        -tl|--leading) LEADING="$2"; shift;;
        -tt|--trailing) TRAILING="$2"; shift;;
        -ts|--slidingwindow) SLIDINGWINDOW="$2"; shift;;
        -tm|--minlen) MINLEN="$2"; shift;;
        *) echo $arg \(with value: $2\) is not a valid parameter...; exit 1;;
    esac
    shift
done

# Replace with function which displays the readme.
function help {
    echo ""
    echo "===================================================================="
    echo ""
    echo " RNA-seq pipeline"
    echo " ------------------------------------------------------------------"
    echo " Usage        : StartPipeline.sh -f [fastq] -r [reference] "
    echo "                -g [gff] -o [output] -t [threads]"
    echo " ------------------------------------------------------------------"
    echo " [fastq]      : folder with fastq files, make sure the extension is .fq"
    echo "                or specify with -e,--extension fa"
    echo " [reference]  : reference genome in fasta format"
    echo " [gff]        : gff file"
    echo " [output]     : output folder for generated data"
    echo " [conditions] : file with conditions specified for every read file"
    echo " [threads]    : number of threads to be used"
    echo " ------------------------------------------------------------------"
    echo " Example :"
    echo " ------------------------------------------------------------------"
    echo " $ bash StartPipeline.sh -f /home/samuel/Reads -r /home/samuel/genome/"
    echo "   Streptococcus_mutans.fa -g /home/samuel/gff/Streptococcus_mutans.gff"
    echo "   -o /home/samuel/output -t 8"
    echo ""
    echo "===================================================================="
    echo ""
    exit 1
}

# Output locations
TRIMMEDREADS=$OUTPUT/TrimmedReads
BOWTIEINDEX=$OUTPUT/BowtieIndex
BOWTIEOUTPUT=$OUTPUT/BowtieOutput
COUNTS=$OUTPUT/Counts
RESULTS=$OUTPUT/Results

function makeDirs {
    # Creates all the necessarily directories for the files created
    # by the pipeline.

    mkdir $OUTPUT

    # In case the directory cannot be made, the pipeline will be terminated.
    if [[ $? -ne 0 ]]; then
        echo "Cannot create directory..."
        exit 5
    fi

    mkdir $TRIMMEDREADS
    mkdir $BOWTIEINDEX
    mkdir $BOWTIEOUTPUT
    mkdir $COUNTS
    mkdir $RESULTS
    mkdir $RESULTS/FastQC
}

function syscall {
    # Error handling.
    # Terminates the pipeline in case of an error.
    echo In SYSCALL Running \""$@"\" >&2
    "$@"
    status=$?
    if [ $status -ne 0 ]; then
        echo "Error with $@" >&2
        echo "Status: $status" >&2
        echo "The pipeline has come to an unexpected end with exit status $status"
        echo "Please check the callback..."
        exit $status
    fi
    return $status
}

function trimming {
    # Runs trimmomatic on every file with the specified extension in the
    # specified location.
    for read in $READS*.$EXTENSION; do
        # Remove path from the fastq files
        TRIMOUT=${read##*/}
        echo "Trimming reads from $read..."

        ./FastQC/fastqc $read -outdir $RESULTS/FastQC

        syscall java -jar trimmomatic-0.36.jar SE $read \
        $TRIMMEDREADS/trimmed_$TRIMOUT ILLUMINACLIP:adapters/$ILLUMINACLIP \
        LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$SLIDINGWINDOW \
        MINLEN:$MINLEN -threads $THREADS

        ./FastQC/fastqc $TRIMMEDREADS/trimmed_$TRIMOUT -outdir $RESULTS/FastQC

    done
}

function mapping {
    # This functions builds a index file for the reference genome and
    # maps the trimmed reads against this.

    echo "Indexing reference genome..."

    # Bowtie2 indexing
    syscall bowtie2-build $REFERENCE $BOWTIEINDEX/index

    # Bowtie2 mapping
    # Loop loops over every file in the directory with the trimmed reads
    # and maps it against the reference genome. 
    for tReads in $TRIMMEDREADS/trimmed_*; do
        BASENAME=$(basename $tReads .$EXTENSION) # Strips path and extension from file
        echo "Mapping reads from $tReads..."

        syscall bowtie2 -p $THREADS -x $BOWTIEINDEX/index -q \
        $tReads -S $BOWTIEOUTPUT/$BASENAME.sam
    done
}

function quantification {
    # Quantification of bowtie allignments.
    echo "Quantifying..."
    
    # HTSeq
    # Loops over all the allignment files from bowtie in the output folder,
    # and performs htseq count on these. Counting is done on all
    # allignment files simultaneously.
    for alFile in $BOWTIEOUTPUT/*; do
        BASENAME=$(basename $alFile .sam) # Strip path and extension from file
        syscall htseq-count -t CDS -s no -i Name $alFile $GFF | grep -v ^__\
        > $COUNTS/$BASENAME.txt &
    done
    wait
}

function differentialExpression {
    echo "Differentializing..."

    for line in $(cat $CONDITIONS); do
        CONDITION=${line%,*}
        SAMPLES1=${line#*,}
        SAMPLES2=${SAMPLES1%.*}
        SAMPLESBASENAME=$(basename $SAMPLES2)

        echo $CONDITION,$(realpath $OUTPUT/Counts/trimmed_$SAMPLESBASENAME.txt)
    done
}

function differentialExpression {
    # Performs differential gene expression analysis

    echo "Differentializing..."
    
    # Loops over every line in the given conditions file. Every line consists
    # of the name of a condition, the reads associated with it, separated with
    # a ",". The loop separates the condition from the reads and puts them in
    # a new file with the location to the count files, belonging to the
    # corresponding reads.
    for line in $(cat $CONDITIONS); do
        # Condition, to which the reads belong to:
        CONDITION=${line%,*}

        # Remove path and extension from the sample files, so the basename can be used
        # to specify the location of the count tables.
        SAMPLES1=${line#*,}
        SAMPLES2=${SAMPLES1%.*}
        SAMPLESBASENAME=$(basename $SAMPLES2) # Remove path and extension from files

        # Combines the conditions, and sample basenames with the location and extension of the
        # count tables. Outputted file will be used by differential gene expression analysis script to link
        # the conditions to the count tables.
        echo $CONDITION,$(realpath $OUTPUT/Counts/trimmed_$SAMPLESBASENAME.txt) \
        >> $OUTPUT/conditions.txt
    done

    python3 gffRead.py $GFF > $OUTPUT/gffInfo.csv
    # Differential expression analysis
    syscall Rscript DGEanalysis.r -f $OUTPUT/conditions.txt \
    -m $METHOD -o $RESULTS/ -s $SPECIES -g $OUTPUT/gffInfo.csv
}

makeDirs
trimming
mapping
quantification
differentialExpression

echo "Done"
