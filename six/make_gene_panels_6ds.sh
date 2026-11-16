#!/bin/bash
# Panel B of the gene figures (six-dataset events; cluster IDs of the id50_6ds event set) and the assembled figures.
set -e
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116/six
P=/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/figures/panels
PY=../.venv/bin/python
mkdir -p $P
$PY gene_panel_116.py clu_52748 $P/USP16_clu_52748 --title USP16 --colors '#2c5f8e,#a8c8e8' --no-legend
$PY gene_panel_116.py clu_86301 $P/PUF60_clu_86301 --title PUF60
$PY gene_panel_116.py clu_65288 $P/RPL22L1_clu_65288 --title RPL22L1 --colors '#a8c8e8,#2c5f8e'
$PY gene_panel_116.py clu_52538 $P/UCKL1_clu_52538 --title UCKL1
$PY gene_panel_116.py clu_91182 $P/DOCK8_clu_91182 --title DOCK8 --colors '#2c5f8e,#a8c8e8'
$PY gene_panel_116.py clu_44871 $P/PTPRC_clu_44871 --title PTPRC --colors '#1f77b4,#7f7f7f,#aec7e8,#ff7f0e,#ffbb78,#2ca02c,#9467bd,#98df8a,#d62728'
$PY gene_panel_116.py clu_63292 $P/CD47_clu_63292 --title CD47
$PY gene_diagram_116.py clu_12789 $P/EIF4B_clu_12789_A --colors '#2c5f8e,#a8c8e8' --below 0
$PY gene_panel_116.py clu_12789 $P/EIF4B_clu_12789 --title EIF4B --colors '#2c5f8e,#a8c8e8' --no-legend
$PY assemble_gene_figures_116.py
