#!/usr/bin/env python3
"""
Human-mouse orthologous exons, Ensembl release 116.

Replaces transcripts_exons_one_gene_HN30.py (one REST-based qsub job per gene).
Annotation comes from the release-116 GTF files and sequences from the local
GRCh38 / GRCm39 genome FASTA, so the run is pinned to one release and is
reproducible.

Logic kept from HN30:
  - exon sequence = its CDS part, taken from the first transcript containing the
    exon in the order: canonical, protein coding by length, others by length
  - global PairwiseAligner, match 2, mismatch -1, open gap -12, extend gap -4
  - %id = identical aligned positions / exon length; a human and a mouse exon
    are orthologous when %id >= 50 when normalised by each of the two lengths
Fixed relative to HN30:
  - REST FASTA line breaks were kept inside the aligned sequences
  - get_cds trimmed only one side of exons with UTR at both ends, and returned an
    arbitrary prefix (or nothing) for exons entirely in UTR; such exons now use
    their full sequence
  - only human x mouse pairs are aligned (HN30 also aligned within species)
Added:
  - all orthologous pairs (min %id >= 50) are saved with lengths, ranks and
    reciprocal-best flags, so stricter identity cutoffs can be applied later

Usage:
  exon_orthology_116.py [--threads N] [--genes i,j,...]
Outputs (ensembl116/orthology/):
  exon_pairs_116.tsv     every orthologous exon pair
  exons_h_116.tsv, exons_m_116.tsv   all exons of the ortholog genes
  genes_status_116.tsv   per gene: exon counts, orthologous exons, status
"""
import argparse
import gzip
import os
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner

MAIN_DIR = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
GTF = {'h': MAIN_DIR + 'ref/Homo_sapiens.GRCh38.116.gtf.gz',
       'm': MAIN_DIR + 'ref/Mus_musculus.GRCm39.116.gtf.gz'}
FASTA = {'h': '/gpfs0/tals/projects/data/Genomes/hg38/grch38/GRCh38.primary_assembly.genome.fa',
         'm': '/gpfs0/tals/projects/data/Genomes/mouse/GRCm39/GRCm39.genome.fa'}
ORTHOLOGS = MAIN_DIR + 'all/good_ortologs_df.txt'
OUT_DIR = MAIN_DIR + 'orthology/'
MIN_ID = 50  # min identity for orthologous exons

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


# ---------------------------------------------------------------- annotation
def parse_gtf(path, wanted_genes):
    """Return {gene_id: gene dict} for the wanted genes."""
    genes = {}
    tr = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            feature = p[2]
            if feature not in ('gene', 'transcript', 'exon', 'CDS', 'stop_codon'):
                continue
            gid = line.split('gene_id "', 1)[1].split('"', 1)[0]
            if gid not in wanted_genes:
                continue
            start, end = int(p[3]), int(p[4])
            if feature == 'gene':
                a = dict(ATTR_RE.findall(p[8]))
                genes[gid] = {'chrom': p[0], 'strand': 1 if p[6] == '+' else -1,
                              'name': a.get('gene_name', 'novel_gene'), 'transcripts': []}
                continue
            tid = line.split('transcript_id "', 1)[1].split('"', 1)[0]
            if feature == 'transcript':
                tags = re.findall(r'tag "([^"]*)"', p[8])
                biotype = re.search(r'transcript_biotype "([^"]*)"', p[8]).group(1)
                tr[tid] = {'id': tid, 'gene': gid, 'biotype': biotype,
                           'is_canonical': 'Ensembl_canonical' in tags,
                           'start': start, 'end': end, 'cds': None, 'exons': []}
            elif feature == 'exon':
                eid = re.search(r'exon_id "([^"]*)"', p[8]).group(1)
                tr[tid]['exons'].append((eid, start, end))
            else:  # CDS / stop_codon: translation span incl. stop, as in Ensembl Translation
                c = tr[tid]['cds']
                tr[tid]['cds'] = (start, end) if c is None else (min(c[0], start), max(c[1], end))
    for t in tr.values():
        g = genes[t['gene']]
        # transcription order of exons
        t['exons'].sort(key=lambda e: e[1], reverse=(g['strand'] == -1))
        t['length'] = sum(e[2] - e[1] + 1 for e in t['exons'])
        g['transcripts'].append(t)
    return genes


def reorder_transcripts(transcripts):
    canonical = [t for t in transcripts if t['is_canonical']]
    protein_coding = sorted([t for t in transcripts if not t['is_canonical'] and t['biotype'] == 'protein_coding'],
                            key=lambda t: t['length'], reverse=True)
    others = sorted([t for t in transcripts if not t['is_canonical'] and t['biotype'] != 'protein_coding'],
                    key=lambda t: t['length'], reverse=True)
    return canonical + protein_coding + others


# ---------------------------------------------------------------- sequence
class Fasta:
    """Minimal samtools-faidx style random access."""
    def __init__(self, path):
        self.f = open(path, 'rb')
        self.idx = {}
        with open(path + '.fai') as fai:
            for line in fai:
                name, length, offset, lbases, lbytes = line.split('\t')[:5]
                self.idx[name] = (int(length), int(offset), int(lbases), int(lbytes))

    def fetch(self, chrom, start, end):  # 1-based inclusive
        length, offset, lbases, lbytes = self.idx[chrom]
        s, e = start - 1, min(end, length)
        bs = offset + (s // lbases) * lbytes + s % lbases
        be = offset + ((e - 1) // lbases) * lbytes + (e - 1) % lbases + 1
        self.f.seek(bs)
        return self.f.read(be - bs).replace(b'\n', b'').decode().upper()


def fasta_chrom(chrom):
    if chrom == 'MT':
        return 'chrM'
    if chrom.isdigit() or chrom in ('X', 'Y'):
        return 'chr' + chrom
    return chrom


COMP = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(COMP)[::-1]


def exon_sequences(gene, fasta):
    """For each exon: the sequence used for alignment, its genomic span and seq type."""
    chrom = fasta_chrom(gene['chrom'])
    strand = gene['strand']
    seqs = {}
    for t in reorder_transcripts(gene['transcripts']):
        for eid, es, ee in t['exons']:
            if eid in seqs:
                continue  # each exon gets its sequence once, from the first transcript
            s, e, kind = es, ee, 'full'
            if t['cds'] is not None:
                cs, ce = t['cds']
                if es <= ce and ee >= cs:  # exon overlaps the CDS
                    s, e = max(es, cs), min(ee, ce)
                    kind = 'cds' if (s, e) != (es, ee) else 'full_cds'
                else:
                    kind = 'utr_full'
            seq = fasta.fetch(chrom, s, e)
            if strand == -1:
                seq = revcomp(seq)
            seqs[eid] = {'seq': seq, 'start': es, 'end': ee, 'seq_type': kind}
    return seqs


def canonical_ranks(gene):
    ranks = {}
    for t in gene['transcripts']:
        if t['is_canonical']:
            n = len(t['exons'])
            for i, (eid, _, _) in enumerate(t['exons']):
                ranks[eid] = f"E{i + 1}Last" if i == n - 1 else f"E{i + 1}"
    return ranks


# ---------------------------------------------------------------- alignment
def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -12
    aligner.extend_gap_score = -4
    return aligner


def identities(alignment, s1, s2):
    a = np.frombuffer(s1.encode(), dtype=np.uint8)
    b = np.frombuffer(s2.encode(), dtype=np.uint8)
    n = 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        n += int((a[a0:a1] == b[b0:b1]).sum())
    return n


def boundary_positions(chrom, strand, start, end):
    """Exon 5' and 3' boundaries, as 'chrN:pos' (HN30 junction_position)."""
    c = 'chr' + chrom
    return (f"{c}:{start}", f"{c}:{end}") if strand == 1 else (f"{c}:{end}", f"{c}:{start}")


def process_gene(args):
    i, row = args
    gh, gm = GENES['h'].get(row['hsapiens_homolog_ensembl_gene']), GENES['m'].get(row['ensembl_gene_id'])
    status = {'i': i, 'ensembl_h': row['hsapiens_homolog_ensembl_gene'], 'ensembl_m': row['ensembl_gene_id']}
    if gh is None or gm is None:
        status['status'] = 'not_in_gtf'
        return status, [], [], []
    symbol = gh['name']
    status.update(symbol=symbol, chrom_h=gh['chrom'], chrom_m=gm['chrom'])
    seq_h, seq_m = exon_sequences(gh, FASTAS['h']), exon_sequences(gm, FASTAS['m'])
    rank_h, rank_m = canonical_ranks(gh), canonical_ranks(gm)
    status.update(exons_h=len(seq_h), exons_m=len(seq_m),
                  transcripts_h=len(gh['transcripts']), transcripts_m=len(gm['transcripts']))
    exons_h = [dict(ensembl=status['ensembl_h'], exon=e, rank=rank_h.get(e, '*'), start=v['start'], end=v['end'],
                    seq_type=v['seq_type'], seq_len=len(v['seq'])) for e, v in seq_h.items()]
    exons_m = [dict(ensembl=status['ensembl_m'], exon=e, rank=rank_m.get(e, '*'), start=v['start'], end=v['end'],
                    seq_type=v['seq_type'], seq_len=len(v['seq'])) for e, v in seq_m.items()]
    if len(seq_h) < 2 or len(seq_m) < 2:
        status['status'] = 'less_than_2_exons'
        return status, exons_h, exons_m, []
    aligner = make_aligner()
    pairs = []
    for eh, vh in seq_h.items():
        sh = vh['seq']
        for em, vm in seq_m.items():
            sm = vm['seq']
            lh, lm = len(sh), len(sm)
            # identical positions <= shorter length, so the pair cannot reach MIN_ID of the longer one
            if lh == 0 or lm == 0 or min(lh, lm) * 100 < MIN_ID * max(lh, lm):
                continue
            aln = aligner.align(sh, sm)[0]
            n = identities(aln, sh, sm)
            id_h, id_m = round(n / lh * 100), round(n / lm * 100)
            if id_h < MIN_ID or id_m < MIN_ID:
                continue
            avg = (id_h + id_m) / 2
            start_h, end_h = boundary_positions(gh['chrom'], gh['strand'], vh['start'], vh['end'])
            start_m, end_m = boundary_positions(gm['chrom'], gm['strand'], vm['start'], vm['end'])
            pairs.append(dict(ensembl=status['ensembl_h'], ensembl_m=status['ensembl_m'], symbol=symbol,
                              exon_h=eh, rank_h=rank_h.get(eh, '*'), start_h=start_h, end_h=end_h,
                              exon_m=em, rank_m=rank_m.get(em, '*'), start_m=start_m, end_m=end_m,
                              len_h=lh, len_m=lm, seq_type_h=vh['seq_type'], seq_type_m=vm['seq_type'],
                              identities=n, align_score=aln.score,
                              **{'%id_human': id_h, '%id_mouse': id_m, 'Average %id': avg,
                                 'Ratio_based_score': avg * min(id_h, id_m) / max(id_h, id_m)}))
    status['ortholog_pairs'] = len(pairs)
    status['exons_with_ortholog_h'] = len({p['exon_h'] for p in pairs})
    status['exons_with_ortholog_m'] = len({p['exon_m'] for p in pairs})
    status['status'] = 'OK' if len(pairs) > 1 else 'less_than_2_ortholog_exons'
    return status, exons_h, exons_m, pairs


def mark_reciprocal_best(pairs):
    """Flag pairs where each exon is the other's best match by Ratio_based_score."""
    if pairs.empty:
        pairs['reciprocal_best'] = []
        return pairs
    best_h = pairs.groupby('exon_h')['Ratio_based_score'].transform('max')
    best_m = pairs.groupby('exon_m')['Ratio_based_score'].transform('max')
    pairs['reciprocal_best'] = (pairs['Ratio_based_score'] == best_h) & (pairs['Ratio_based_score'] == best_m)
    return pairs


GENES, FASTAS = {}, {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--genes', help='comma separated row indices of good_ortologs_df (for testing)')
    ap.add_argument('--out-dir', default=OUT_DIR)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    orthologs = pd.read_csv(ORTHOLOGS, sep='\t', index_col=0, dtype=str)
    if args.genes:
        orthologs = orthologs.loc[[int(x) for x in args.genes.split(',')]]
    for sp, col in (('h', 'hsapiens_homolog_ensembl_gene'), ('m', 'ensembl_gene_id')):
        print('parsing', GTF[sp], flush=True)
        GENES[sp] = parse_gtf(GTF[sp], set(orthologs[col]))
        FASTAS[sp] = Fasta(FASTA[sp])

    statuses, all_h, all_m, all_pairs = [], [], [], []
    with Pool(args.threads) as pool:  # GENES/FASTAS are inherited by fork
        for k, (st, eh, em, pr) in enumerate(pool.imap_unordered(process_gene, orthologs.iterrows(), chunksize=4)):
            statuses.append(st)
            all_h += eh
            all_m += em
            all_pairs += pr
            if k % 500 == 0:
                print(k, 'genes done', flush=True)

    pairs = pd.DataFrame(all_pairs)
    pairs = pairs.groupby('ensembl', group_keys=False).apply(mark_reciprocal_best) if len(pairs) else pairs
    pairs.to_csv(args.out_dir + 'exon_pairs_116.tsv', sep='\t', index=False)
    pd.DataFrame(all_h).to_csv(args.out_dir + 'exons_h_116.tsv', sep='\t', index=False)
    pd.DataFrame(all_m).to_csv(args.out_dir + 'exons_m_116.tsv', sep='\t', index=False)
    st = pd.DataFrame(statuses).sort_values('i')
    st.to_csv(args.out_dir + 'genes_status_116.tsv', sep='\t', index=False)
    print(st['status'].value_counts())


if __name__ == '__main__':
    main()
