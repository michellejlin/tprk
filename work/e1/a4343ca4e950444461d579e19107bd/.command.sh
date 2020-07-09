#!/bin/bash -ue
echo metadata_pacbio.csv
Rscript /Users/uwvirongs/Documents/Michelle/tprk_nextflow/tprk/og_files_to_all_reads.R -s /Users/uwvirongs/Documents/Michelle/tprk_nextflow/tprk -d example/ -m metadata_pacbio.csv --pacbio
