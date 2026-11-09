#!/usr/bin/env python3
"""
Fetch mouse->human orthologs from Ensembl BioMart, pinned to release 116.

Replaces get_homologous_table_HN1.py, which queried the live BioMart without
checking the release. Fails if the BioMart host is not serving release 116.

Outputs (ensembl116/all/):
  mouse_orthologs_v116.txt   all mouse genes with their human homologs
  good_ortologs_df.txt       high-confidence one-to-one orthologs
"""
import io
import os
import sys

import pandas as pd
import requests

RELEASE = 116
# Release-116 archive host (www.ensembl.org BioMart redirects here)
BIOMART_HOST = "https://jun2026.archive.ensembl.org"
WORKING_DIR = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/all/'

ATTRIBUTES = [
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_perc_id",
    "hsapiens_homolog_perc_id_r1",
    "hsapiens_homolog_goc_score",
    "hsapiens_homolog_wga_coverage",
    "hsapiens_homolog_orthology_confidence",
]


def check_release():
    r = requests.get(BIOMART_HOST + "/biomart/martservice?type=registry", timeout=60)
    r.raise_for_status()
    tag = f"Ensembl Genes {RELEASE}"
    if tag not in r.text:
        sys.exit(f"BioMart at {BIOMART_HOST} does not serve '{tag}'. Use the release-{RELEASE} archive host.")
    print(f"BioMart host {BIOMART_HOST} serves {tag}")


def query_biomart():
    attrs = "".join(f'<Attribute name="{a}"/>' for a in ATTRIBUTES)
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
           '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">'
           f'<Dataset name="mmusculus_gene_ensembl" interface="default">{attrs}</Dataset></Query>')
    r = requests.get(BIOMART_HOST + "/biomart/martservice", params={"query": xml}, timeout=1800)
    r.raise_for_status()
    if "ERROR" in r.text[:200]:
        sys.exit("BioMart error: " + r.text[:500])
    df = pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
    df.columns = ATTRIBUTES
    return df


def main():
    os.makedirs(WORKING_DIR, exist_ok=True)
    check_release()
    mouse_ortologs = query_biomart()
    mouse_ortologs.to_csv(WORKING_DIR + f'mouse_orthologs_v{RELEASE}.txt', index=False)
    print("there are", mouse_ortologs['ensembl_gene_id'].nunique(), "genes in mouse")
    filtered_df = mouse_ortologs[(mouse_ortologs['hsapiens_homolog_orthology_confidence'] == "1") &
                                 (mouse_ortologs['hsapiens_homolog_orthology_type'] == 'ortholog_one2one')]
    filtered_df = filtered_df.drop_duplicates(['ensembl_gene_id', 'hsapiens_homolog_ensembl_gene']).reset_index(drop=True)
    print("there are", filtered_df['ensembl_gene_id'].nunique(), "good orthologs genes in mouse")
    filtered_df.to_csv(WORKING_DIR + 'good_ortologs_df.txt', sep="\t")
    print("Data saved")


if __name__ == "__main__":
    main()
