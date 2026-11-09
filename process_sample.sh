#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -pe shared 8
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# Download, align and extract junctions for one row of the new-datasets sample sheet.
# Same steps and parameters as scripts/GSE*/execute_single_bam_junc_HN1.sh
# (HISAT2 2.0.5 --add-chrname, samtools, regtools junctions extract -s 0 -a 8 -m 50 -M 500000),
# except that paired-end data are aligned as pairs and FASTQs come from ENA / ENCODE.
#
# One task per sample sheet row; submit_samples.sh submits them, at most 5 running at once:
#   submit_samples.sh 1 33   (or a single row: qsub -t 3 process_sample.sh)
set -euo pipefail

SHEET=/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/metadata/sample_sheet_new_datasets.tsv
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
REGTOOLS=/gpfs0/tals/projects/software/leafcutter_0.2.9/bin/regtools
FASTERQ_DUMP=/gpfs0/tals/projects/software/anaconda3/envs/SRA/bin/fasterq-dump
CUTADAPT=/gpfs0/tals/projects/software/anaconda3/envs/cutadaptenv/bin/cutadapt
declare -A INDEX=([human]=/gpfs0/tals/projects/data/Genomes/hg38/grch38/genome
                  [mouse]=/gpfs0/tals/projects/data/Genomes/mouse/GRCm39/genome)

row=$(awk -F'\t' -v n=$((SGE_TASK_ID + 1)) 'NR==n' $SHEET)
IFS=$'\t' read -r dataset study species cell sample runs layout read_count fastq_1 fastq_2 <<< "$row"
echo "$SGE_TASK_ID $study $sample $layout"

bam_folder=$BASE/$study/bam_files
junc_folder=$BASE/$study/leafcutter_0.2.9
fq_folder=$bam_folder/fastq_$sample
mkdir -p $bam_folder $junc_folder $fq_folder
if [[ -s $junc_folder/$sample.sort.bam.junc ]]; then
    echo "already done"; exit 0
fi

download() {  # ';'-separated URLs -> ','-separated local paths
    local paths=()
    IFS=';' read -ra urls <<< "$1"
    for u in "${urls[@]}"; do
        local f=$fq_folder/$(basename "${u%%\?*}")
        wget -q -c --tries=10 -O "$f" "$u"
        gzip -t "$f"
        paths+=("$f")
    done
    (IFS=,; echo "${paths[*]}")
}

fetch_ranges() {  # url out: download as NPARTS parallel byte ranges (the bandwidth limit here is per
                  # connection), continuing from a partial file left by an earlier attempt
    local url=$1 out=$2 nparts=8
    local total=$(curl -sI --retry 5 "$url" | awk 'tolower($1)=="content-length:"{print $2+0}' | tail -n 1)
    local have=0
    [[ -f $out ]] && have=$(stat -c %s $out)
    (( have < total )) || return 0
    local chunk=$(( (total - have + nparts - 1) / nparts ))
    for i in $(seq 0 $((nparts - 1))); do
        local s=$((have + i * chunk)) e=$((have + (i + 1) * chunk - 1))
        (( e >= total )) && e=$((total - 1))
        (( s > e )) && continue
        (
            for try in 1 2 3 4 5 6 7 8 9 10; do
                curl -s --retry 10 --retry-delay 10 -r $s-$e -o $out.part$i "$url" || true
                [[ -f $out.part$i && $(stat -c %s $out.part$i) -eq $((e - s + 1)) ]] && exit 0
                sleep 30
            done
            exit 1
        ) &
    done
    wait
    local got=$have
    for i in $(seq 0 $((nparts - 1))); do
        [[ -f $out.part$i ]] && got=$((got + $(stat -c %s $out.part$i)))
    done
    (( got == total )) || { echo "incomplete download $url ($got of $total bytes)"; exit 1; }
    for i in $(seq 0 $((nparts - 1))); do
        [[ -f $out.part$i ]] && cat $out.part$i >> $out && rm -f $out.part$i
    done
}

sra_runs() {  # SRR runs: .sra from the NCBI AWS mirror (much faster here than ENA), then fasterq-dump
    local r1=() r2=()
    IFS=',' read -ra runs <<< "$1"
    for run in "${runs[@]}"; do
        if [[ ! -s $fq_folder/${run}_1.fastq && ! -s $fq_folder/$run.fastq ]]; then  # not dumped by an earlier attempt
            fetch_ranges https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run $fq_folder/$run.sra
            $FASTERQ_DUMP --split-files --threads 8 --temp $fq_folder -O $fq_folder $fq_folder/$run.sra
            rm -f $fq_folder/$run.sra
        fi
        # some submissions store index reads as reads 1-2, so the biological pair comes out as _3/_4
        if [[ ! -s $fq_folder/${run}_1.fastq && -s $fq_folder/${run}_3.fastq && -s $fq_folder/${run}_4.fastq ]]; then
            mv $fq_folder/${run}_3.fastq $fq_folder/${run}_1.fastq
            mv $fq_folder/${run}_4.fastq $fq_folder/${run}_2.fastq
        fi
        if [[ "$layout" == "PAIRED" ]]; then
            r1+=($fq_folder/${run}_1.fastq); r2+=($fq_folder/${run}_2.fastq)
        else
            r1+=($(ls $fq_folder/${run}*.fastq))
        fi
    done
    (IFS=,; echo "${r1[*]}") > $fq_folder/r1.txt
    (IFS=,; echo "${r2[*]}") > $fq_folder/r2.txt
}

if [[ "$runs" == SRR* ]]; then
    sra_runs "$runs"
    r1=$(cat $fq_folder/r1.txt)
    r2=$(cat $fq_folder/r2.txt)
else
    r1=$(download "$fastq_1")
    [[ "$layout" == "PAIRED" ]] && r2=$(download "$fastq_2")
fi
# Adapter trimming (not done for the original datasets): Smart-seq2 libraries (GSE118165) carry
# Nextera adapters and the ISPCR/TSO oligo at read starts, which otherwise leave ~25% of reads unmapped.
NEXTERA=CTGTCTCTTATACACATCT; TRUSEQ=AGATCGGAAGAGC; ISPCR=AAGCAGTGGTATCAACGCAGAGT
trim() {  # ','-separated R1 list, R2 list (may be empty) -> trimmed lists in t1.txt / t2.txt
    IFS=',' read -ra a <<< "$1"; IFS=',' read -ra b <<< "$2"
    local t1=() t2=()
    # minimum read length after trimming: 30, or 18 for short-read libraries (e.g. ImmGen 2x25)
    local readlen=$(head -2 ${a[0]} | tail -1 | tr -d '\n' | wc -c)
    local minlen=30
    (( readlen < 40 )) && minlen=18
    for k in "${!a[@]}"; do
        local o1=$fq_folder/trim_${k}_1.fastq o2=$fq_folder/trim_${k}_2.fastq
        if [[ "$layout" == "PAIRED" ]]; then
            $CUTADAPT -j 8 -q 20 -m $minlen -a $NEXTERA -a $TRUSEQ -g $ISPCR -A $NEXTERA -A $TRUSEQ -G $ISPCR \
                -o $o1 -p $o2 ${a[$k]} ${b[$k]} > $bam_folder/$sample.cutadapt_$k.log
            t2+=($o2)
        else
            $CUTADAPT -j 8 -q 20 -m $minlen -a $NEXTERA -a $TRUSEQ -g $ISPCR -o $o1 ${a[$k]} > $bam_folder/$sample.cutadapt_$k.log
        fi
        t1+=($o1)
        rm -f ${a[$k]} ${b[$k]:-}
    done
    (IFS=,; echo "${t1[*]}") > $fq_folder/t1.txt
    (IFS=,; echo "${t2[*]}") > $fq_folder/t2.txt
}
trim "$r1" "${r2:-}"
r1=$(cat $fq_folder/t1.txt)
r2=$(cat $fq_folder/t2.txt)

if [[ "$layout" == "PAIRED" ]]; then
    hisat2 -p 8 -x ${INDEX[$species]} --add-chrname --mm -1 $r1 -2 $r2 2> $bam_folder/$sample.hisat2.log \
        | samtools sort -@ 4 -o $bam_folder/$sample.sort.bam -
else
    hisat2 -p 8 -x ${INDEX[$species]} --add-chrname --mm -U $r1 2> $bam_folder/$sample.hisat2.log \
        | samtools sort -@ 4 -o $bam_folder/$sample.sort.bam -
fi
samtools index $bam_folder/$sample.sort.bam
samtools flagstat $bam_folder/$sample.sort.bam > $bam_folder/$sample.flagstat.txt

cd $bam_folder
$REGTOOLS junctions extract -s 0 -a 8 -m 50 -M 500000 $sample.sort.bam -o $junc_folder/$sample.sort.bam.junc
rm -rf $fq_folder
echo done
