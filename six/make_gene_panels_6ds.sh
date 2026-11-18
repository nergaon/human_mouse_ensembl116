#!/bin/bash
# Panel B of the gene figures (six-dataset events; cluster IDs of the id50_6ds event set) and the assembled figures.
set -e
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116/six
P=/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/figures/panels
PY=../.venv/bin/python
mkdir -p $P
$PY gene_panel_116.py clu_52830 $P/USP16_clu_52830 --title USP16 --colors '#2c5f8e,#a8c8e8' --no-legend
$PY gene_panel_116.py clu_86433 $P/PUF60_clu_86433 --title PUF60
$PY gene_panel_116.py clu_65387 $P/RPL22L1_clu_65387 --title RPL22L1 --colors '#a8c8e8,#2c5f8e'
$PY gene_panel_116.py clu_52620 $P/UCKL1_clu_52620 --title UCKL1
$PY gene_panel_116.py clu_91320 $P/DOCK8_clu_91320 --title DOCK8 --colors '#2c5f8e,#a8c8e8'
$PY gene_panel_116.py clu_44935 $P/PTPRC_clu_44935 --title PTPRC --colors '#1f77b4,#7f7f7f,#aec7e8,#ff7f0e,#ffbb78,#2ca02c,#9467bd,#98df8a,#d62728'
$PY gene_panel_116.py clu_63390 $P/CD47_clu_63390 --title CD47
$PY gene_diagram_116.py clu_12809 $P/EIF4B_clu_12809_A --colors '#2c5f8e,#a8c8e8' --below 0
$PY gene_panel_116.py clu_12809 $P/EIF4B_clu_12809 --title EIF4B --colors '#2c5f8e,#a8c8e8' --no-legend
# candidate figures (not in the manuscript)
C3='#2c5f8e,#a8c8e8,#e08a3c'
$PY gene_diagram_116.py clu_91320 $P/DOCK8_clu_91320_A --colors '#2c5f8e,#a8c8e8' >/dev/null || true
$PY gene_diagram_116.py clu_23540 $P/CIITA_clu_23540_A --colors '#2c5f8e,#a8c8e8' --below 1
$PY gene_panel_116.py clu_23540 $P/CIITA_clu_23540 --title CIITA --colors '#2c5f8e,#a8c8e8' --no-legend
$PY gene_diagram_116.py clu_23547 $P/CIITA_clu_23547_A --colors "$C3"
$PY gene_panel_116.py clu_23547 $P/CIITA_clu_23547 --title CIITA --colors "$C3" --no-legend
$PY gene_diagram_116.py clu_70891 $P/KLF3_clu_70891_A --colors "$C3"
$PY gene_panel_116.py clu_70891 $P/KLF3_clu_70891 --title KLF3 --colors "$C3" --no-legend
$PY gene_diagram_116.py clu_42544 $P/PSMD4_clu_42544_A --colors '#2c5f8e,#a8c8e8' --below 1
$PY gene_panel_116.py clu_42544 $P/PSMD4_clu_42544 --title PSMD4 --colors '#2c5f8e,#a8c8e8' --no-legend
$PY gene_diagram_116.py clu_11399 $P/ATF7IP_clu_11399_A --colors '#2c5f8e,#a8c8e8'
$PY gene_panel_116.py clu_11399 $P/ATF7IP_clu_11399 --title ATF7IP --colors '#2c5f8e,#a8c8e8' --no-legend
$PY assemble_gene_figures_116.py
