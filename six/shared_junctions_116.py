#!/usr/bin/env python3
"""
Shared orthologous junctions and AS clusters across several human and mouse datasets.

Generalises find_shared_junctions.py (which required presence in all 4 datasets)
to any number of datasets. Inputs are *_orthologs_junctions.tsv tables in human
coordinates (output of filter_ortholog_junctions.py). A dataset may consist of
several tables (the composite mouse M2 set), which are merged.

A junction is kept when
  - it is present (>= 1 read) in at least --min-human human and --min-mouse mouse datasets, and
  - it has >= --min-reads reads in >= --min-samples samples in >= --min-expressed datasets.
Kept junctions are clustered (shared start or end, assign_cluster_HN1_fixed) and
clusters with > 1 junction are written as the LeafCutter count table.
--extra tables (PBMC, fibroblasts) do not take part in the presence rule or the clustering:
their counts are added for the same junctions (0 where absent), so every comparison tests the
same events (as merge_fibroblast_clusters.py did for the fibroblasts).

Usage:
  shared_junctions_116.py --table H1:human:GSE60424_orthologs_junctions.tsv ... \
      [--extra H4:human:GSE107011_orthologs_junctions.tsv ...] --out-dir DIR
Outputs (DIR/):
  shared_junctions_filtered.tsv, shared_junctions_cluster_all.tsv,
  AS_clusters_value.txt (LeafCutter input; columns <sample>_<study>),
  junction_presence.tsv (per junction: present / expressed per dataset)
"""
import argparse
import os
from collections import defaultdict

import pandas as pd

from assign_cluster_HN1_fixed import assign_clusters_to_groups, write_leafcutter_value_file


def normalize_junction_id(junction_id):
    chrom, start, end = junction_id.split(':')[:3]
    a, b = sorted((int(start), int(end)))
    return f'{chrom}:{a}:{b}'


def load_table(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    df = df.drop(columns=[c for c in df.columns if c == 'original_junction_id'])
    df.index = [normalize_junction_id(j) for j in df.index]
    df = df.groupby(level=0).sum()  # junctions differing only by strand
    study = os.path.basename(path).split('_orthologs_junctions')[0]
    df.columns = [f'{study}_{c}' for c in df.columns]  # assign_cluster expects <dataset>_<sample>
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', action='append', required=True, help='LABEL:SPECIES:PATH (repeatable)')
    ap.add_argument('--extra', action='append', default=[], help='LABEL:SPECIES:PATH, counts added only')
    ap.add_argument('--min-human', type=int, default=2)
    ap.add_argument('--min-mouse', type=int, default=2)
    ap.add_argument('--min-reads', type=int, default=10)
    ap.add_argument('--min-samples', type=int, default=2)
    ap.add_argument('--min-expressed', type=int, default=2)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    species, tables = {}, defaultdict(list)
    for spec in args.table:
        label, sp, path = spec.split(':', 2)
        species[label] = sp
        tables[label].append(load_table(path))
        print(label, sp, path, tables[label][-1].shape, flush=True)
    datasets = {label: pd.concat(t, axis=1).fillna(0) for label, t in tables.items()}

    presence = pd.DataFrame({label: pd.Series(True, index=d.index) for label, d in datasets.items()})
    presence = presence.notna() & presence.fillna(False).astype(bool)
    expressed = pd.DataFrame({label: (d >= args.min_reads).sum(axis=1) >= args.min_samples
                              for label, d in datasets.items()}).reindex(presence.index)
    expressed = expressed.fillna(False).astype(bool)
    n_human = presence[[l for l in presence if species[l] == 'human']].sum(axis=1)
    n_mouse = presence[[l for l in presence if species[l] == 'mouse']].sum(axis=1)
    keep = (n_human >= args.min_human) & (n_mouse >= args.min_mouse) & (expressed.sum(axis=1) >= args.min_expressed)
    print(f'{len(presence)} orthologous junctions; {keep.sum()} kept', flush=True)
    pd.concat([presence.add_prefix('present_'), expressed.add_prefix('expressed_')], axis=1) \
        .assign(kept=keep).to_csv(os.path.join(args.out_dir, 'junction_presence.tsv'), sep='\t')

    merged = pd.concat([d.reindex(keep[keep].index) for d in datasets.values()], axis=1).fillna(0).astype(int)
    merged.index.name = 'junction_id'
    merged.to_csv(os.path.join(args.out_dir, 'shared_junctions_filtered.tsv'), sep='\t')

    all_clusters, as_clusters = assign_clusters_to_groups(merged, 'junction_id', 0, min_psi=0.0)
    as_clusters = as_clusters.reset_index().sort_values(['cluster', 'junction_id']).set_index('junction_id')
    all_clusters.to_csv(os.path.join(args.out_dir, 'shared_junctions_cluster_all.tsv'), sep='\t')
    for spec in args.extra:
        label, sp, path = spec.split(':', 2)
        extra = load_table(path).reindex(as_clusters.index).fillna(0).astype(int)
        print(label, sp, path, f'{(extra.sum(axis=1) > 0).sum()} of {len(extra)} AS junctions covered', flush=True)
        as_clusters = pd.concat([as_clusters, extra], axis=1)
    write_leafcutter_value_file(as_clusters, os.path.join(args.out_dir, 'AS_clusters_value.txt'))
    print(f"{as_clusters['cluster'].nunique()} AS clusters, {len(as_clusters)} junctions")


if __name__ == '__main__':
    main()
