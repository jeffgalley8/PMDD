
#!/bin/bash

# PICRUSt2 functional prediction
# Input: QIIME2-exported representative sequences and BIOM feature table
# Output: picrust_out/

module load miniconda3/24.1.2-py310
conda activate picrust2

picrust2_pipeline.py -s exported-rep-seqs/dna-sequences.fasta -i exported-feature-table/feature-table.biom -o picrust_out -p 4