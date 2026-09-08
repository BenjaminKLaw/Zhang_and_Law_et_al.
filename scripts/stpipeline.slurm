#!/bin/bash
sample=sample_id
FW=path/to/rawdata/sample_bar_umi.fastq.gz
RV=path/to/rawdata/sample_read-3.fastq.gz
MAP=path/to/stpipeline/map/STARindex_grcm39
ANN=path/to/stpipeline/ann/GRCm39_genes.gtf
CONT=path/to/stpipeline/cont/ncRNA/StarIndex
ID=path/to/stpipeline/barcode.xls
OUTPUT=path/to/output/$sample
mkdir -p path/to/output/$sample
TMP=path/to/output/$sample/tmp
mkdir -p path/to/output/$sample/tmp
EXP=$sample

st_pipeline_run.py --output-folder $OUTPUT --temp-folder $TMP --umi-start-position 16 --umi-end-position 26 --ids $ID --ref-map $MAP --ref-annotation $ANN --expName $EXP --htseq-no-ambiguous --verbose --threads 16 --log-file $OUTPUT/${EXP}_log.txt --star-two-pass-mode --no-clean-up --contaminant-index $CONT --disable-clipping --min-length-qual-trimming 30 $FW $RV

convertEnsemblToNames.py --annotation $ANN --output $OUTPUT/$sample\_stdata.updated.tsv $OUTPUT/$sample\_stdata.tsv

