#!/bin/bash

# Usage: map_peaks_to_known_genes.py [options] <knownGene file> <knownGene xRef file> <peaks file>
#
#
# Map the peaks in <peaks file> to genes in <knownGene file>.  <knownGene file>
# is
# format is as specified in
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.sql.
# <peaks file> format is as produced by MACS.  If *auto* is chosen (default)
# file extension is examined for *.xls* for default MACS format or *.bed* for
# BED format.  If the --detailoption is provided, the following extra fields are
# appended to each row:
#
# peak loc, dist from feature, score, map type, map subtype
#
#
# Options:
#   -h, --help            show this help message and exit
#   --upstream-window=UPST_WIN
#                         window width in base pairs to consider promoter region
#                         [default: 5500]
#   --downstream-window=DNST_WIN
#                         window width in base pairs to consider downstream
#                         region [default: 2500]
#   --tss                 calculate downstream window from transcription start
#                         site instead of transcription end site
#   --map-output=PEAK_OUTPUT
#                         filename to output mapped peaks to [default: stdout]
#   --stats-output=STATS_OUTPUT
#                         filename to output summary stats in conversion
#                         [default: stderr]
b#   --peaks-format=PEAKS_FMT
#                         format of peaks input file [default: auto]
#   --detail              add extra fields to output, see description

ORG=mm9
KG_FN=$(org_settings.py $ORG known_gene_anno_path)
XREF_FN=$(org_settings.py $ORG known_gene_xref_path)
OPTS="--detail --tss --upstream-window=10000 --downstream-window=10000"
PEAKS_FN=test_peaks.xls

echo map_peaks_to_known_genes.py $OPTS $KG_FN $XREF_FN $PEAKS_FN
map_peaks_to_known_genes.py $OPTS $KG_FN $XREF_FN $PEAKS_FN