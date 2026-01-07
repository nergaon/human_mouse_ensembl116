#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  4 14:09:34 2026

@author: nergaon
"""

import requests
import sys

def main():
    ensembl_gene_id = "ENSG00000112514"  # human CUTA
    # Allow overriding the Ensembl gene id from the command line
    if len(sys.argv) > 1:
        ensembl_gene_id = sys.argv[1]

    base_url = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{ensembl_gene_id}?expand=1"
    headers = {"Content-Type": "application/json"}

    response = requests.get(base_url + endpoint, headers=headers)
    gene_data = response.json()

    transcripts = gene_data.get("Transcript", [])
    print(f"Found {len(transcripts)} gene-assigned transcripts")  # should be 15

    url = f"{base_url}/overlap/id/{ensembl_gene_id}?feature=transcript"
    response = requests.get(url, headers=headers)
    overlap_data = response.json()

    # Build a dictionary: transcript_id -> tsl
    tsl_map = {t["id"]: t.get("transcript_support_level") for t in overlap_data}

    # Merge TSL info into gene-assigned transcripts
    gene_transcripts_with_tsl = []
    for t in transcripts:
        transcript_id = t["id"]
        tsl = tsl_map.get(transcript_id)
        gene_transcripts_with_tsl.append({
            "id": transcript_id,
            "display_name": t.get("display_name"),
            "biotype": t.get("biotype"),
            "is_canonical": t.get("is_canonical"),
            "tsl": tsl  # will be None if not assigned
        })

    # Print nicely
    for t in gene_transcripts_with_tsl:
        print(t["id"], t["display_name"], t["biotype"], t["is_canonical"], t["tsl"])
    return

if __name__ == "__main__":
    main()
