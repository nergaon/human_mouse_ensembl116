#!/usr/bin/env python3
"""
Supplementary tables (Excel) for the revised manuscript, from the Ensembl 116 / id50 results
(orthology sensitivity tables also use id70, id80 and id50_rbh).

  Table_S1  orthology map per identity cutoff
  Table_S7  junction-level LeafCutter results, one sheet per comparison
  Table_S8  differential splicing output of every comparison x cell type
  Table_S9  high- and low-confidence counts
  Table_S10  event status per cell type (confidence level) and in fibroblasts
  Table_S11 enrichment of high-confidence DS-only / equivalent-only genes, with coverage check
  Table_S5 sensitivity analyses (orthology cutoff, dPSI cutoff, common events, splice-site
            correspondence, matched populations)
  Table_S6 coverage bias model
Output: ensembl116/revision_separation/supplementary/
"""
import os
import re
import shutil
import sys
import zipfile

import numpy as np
import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
V = 'id50'
ED = BASE + f'events/{V}_sep/'
OUT = BASE + 'revision_separation/supplementary/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
DATASETS = {'H1': 'GSE60424', 'H2': 'GSE115736', 'H3': 'GSE118165', 'M1': 'GSE116177',
            'M2': 'GSE287527 / GSE261448 / ENCODE4 / GSE313371', 'FH': 'E-MTAB-5919 human',
            'FM': 'E-MTAB-5919 mouse'}


def write(name, title, note, sheets):
    """sheets: {sheet name: DataFrame}; row 1 title, row 2 note, header in row 3."""
    path = OUT + f'{name}.xlsx'
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        for sheet, df in sheets.items():
            df.to_excel(xw, sheet_name=sheet[:31], index=False, startrow=2)
            ws = xw.sheets[sheet[:31]]
            ws['A1'] = title
            ws['A2'] = note
            ws.freeze_panes = 'A4'
            for col in ws.iter_cols(min_row=3, max_row=min(ws.max_row, 200)):
                width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)
    times_new_roman(path)
    print('wrote', path)


def times_new_roman(path):
    """All fonts of the workbook in Times New Roman (instead of the Calibri theme font of openpyxl)."""
    tmp = path + '.tmp'
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == 'xl/styles.xml':
                x = data.decode()
                x = re.sub(r'<name val="[^"]*"\s*/>', '<name val="Times New Roman"/>', x)
                x = re.sub(r'<family val="\d+"\s*/>', '<family val="1"/>', x)
                x = re.sub(r'<scheme val="\w+"\s*/>', '', x)
                data = x.encode()
            zout.writestr(item, data)
    shutil.move(tmp, path)


def descriptions():
    h = pd.read_csv(BASE + 'all/hgnc_complete_set.txt', sep='\t', usecols=['symbol', 'name', 'ensembl_gene_id'],
                    dtype=str)
    return h.dropna(subset=['ensembl_gene_id']).drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id')['name']


def kind(c):
    if c.startswith('matched_'):
        return 'matched ' + kind(c[len('matched_'):])
    if c == 'FH_FM':
        return 'fibroblasts'
    a, b = c.split('_')
    return 'human-mouse' if a[0] != b[0] else ('human-human' if a[0] == 'H' else 'mouse-mouse')


def table_s4():
    s = pd.concat([pd.read_csv(BASE + f'events/{v}/../../orthology/orthology_summary_116.tsv', sep='\t')
                   for v in ['id50']]).drop_duplicates('variant')
    ev = []
    for v in ['id50', 'id70', 'id80', 'id50_rbh']:
        a = pd.read_csv(BASE + f'events/{v}_sep/classification/events_annotation.tsv', sep='\t')
        ev.append(dict(variant=v, AS_events=len(a), genes=a['ensembl'].nunique(),
                       **a['event_type'].value_counts().to_dict()))
    s = s.merge(pd.DataFrame(ev), on='variant')
    s.insert(1, 'definition', s['variant'].map({
        'id50': 'min(%id human, %id mouse) >= 50 (main analysis)', 'id70': 'min %id >= 70', 'id80': 'min %id >= 80',
        'id50_rbh': '>= 50 and reciprocal best match'}))
    write('Table_S1', 'Supplementary Table S1: Human-mouse exon orthology map (Ensembl 116) at different stringencies',
          '%id = identical aligned positions / exon length, in each species; exon_pairs = human-mouse exon pairs passing '
          'the identity cutoff (many-to-many, because overlapping exon variants with alternative ends are separate '
          'Ensembl exons, so one exon can match several exons of the other species); boundaries = unique orthologous exon '
          'boundaries after selecting the best match per boundary; AS events = orthologous alternative splicing events. '
          'rank_consistency: exon rank (canonical transcript) of the orthologous boundaries used by the events; rank is '
          'not used to define orthology.',
          {'orthology_map': s,
           'rank_consistency': pd.concat([pd.read_csv(BASE + f'orthology/rank_consistency_116_{v}.tsv', sep='\t')
                                          for v in ['id50', 'id70', 'id80', 'id50_rbh']])})


def table_s5(points, ann):
    info = ann.set_index('cluster')
    h2m = dict(zip(points['position_h'], points['position_m']))
    rank_h = dict(zip(points['position_h'], points['rank_h']))
    rank_m = dict(zip(points['position_h'], points['rank_m']))
    sheets = {}
    # only the pairs of the article (matched_ pairwise LeafCutter runs are not reported; the matched populations are
    # analysed by the five-dataset classification, Table S5)
    for subset in ('original',):
        for pair in sorted(os.listdir(ED + f'leafcutter/{subset}')):
            rows = None
            for cell in CELLS + ['Fibro']:
                d = ED + f'leafcutter/{subset}/{pair}/{cell}/'
                if not os.path.exists(d + 'leafcutter_ds_effect_sizes.txt'):
                    continue
                eff = pd.read_csv(d + 'leafcutter_ds_effect_sizes.txt', sep='\t')
                sig = pd.read_csv(d + 'leafcutter_ds_cluster_significance.txt', sep='\t')
                a, b = eff.columns[2], eff.columns[3]
                eff['cluster'] = eff['intron'].str.split(':').str[3]
                sig['cluster'] = sig['cluster'].str.split(':').str[1]
                eff = eff.merge(sig[['cluster', 'p.adjust']], on='cluster', how='left')
                part = eff.set_index('intron')[[a, b, 'deltapsi', 'p.adjust']]
                part.columns = [f'{cell}_PSI_{a}', f'{cell}_PSI_{b}', f'{cell}_deltaPSI_{b}-{a}', f'{cell}_p.adjust']
                rows = part if rows is None else rows.join(part, how='outer')
            if rows is None:
                continue
            idx = rows.index.to_series()
            parts = idx.str.split(':', expand=True)
            first = pd.DataFrame({
                'Human junction': parts[0] + ':' + parts[1] + ':' + parts[2],
                'Mouse junction (orthologous)': [
                    f"{h2m.get(f'{c}:{s}', 'NA')}-{str(h2m.get(f'{c}:{e}', 'NA')).split(':')[-1]}"
                    for c, s, e in zip(parts[0], parts[1], parts[2])],
                'Human exon ranks (junction start-end)': [
                    f"{rank_h.get(f'{c}:{s}', 'NA')}-{rank_h.get(f'{c}:{e}', 'NA')}"
                    for c, s, e in zip(parts[0], parts[1], parts[2])],
                'Mouse exon ranks (junction start-end)': [
                    f"{rank_m.get(f'{c}:{s}', 'NA')}-{rank_m.get(f'{c}:{e}', 'NA')}"
                    for c, s, e in zip(parts[0], parts[1], parts[2])],
                'Alternative splicing event ID': parts[0] + ':' + parts[3],
                'Human gene symbol': parts[3].map(info['symbol']).values,
                'Alternative splicing type': parts[3].map(info['event_type']).values,
                'Well-supported splice-site correspondence': parts[3].map(info['splice_site_correspondence']).values},
                index=rows.index)
            sheets[('matched_' if subset == 'matched' else '') + pair] = pd.concat([first, rows], axis=1)
    sheets = {'summary': pair_summary(ann), **sheets}
    write('Table_S7', 'Supplementary Table S7: Junction-level LeafCutter differential splicing results',
          'Sheet summary: one row per junction of the events testable in at least one cell type. For each cell type: '
          'PSI of the junction in each dataset (mean over the samples of the dataset; junction reads / event reads, '
          'samples with >= 20 event reads), gap = human-mouse separation of the junction (lowest human minus highest '
          'mouse PSI, or lowest mouse minus highest human, whichever is larger; negative = the species overlap), and the '
          'status of the event in that cell type (Table 1; DS = differentially spliced, equivalent, inconclusive or '
          'not testable; high and low confidence combined). An event is differentially spliced when the gap of at least '
          'one junction is >= 0.1 (high confidence >= 0.2), and equivalent when the PSI range of every junction is '
          '<= 0.1 (high confidence <= 0.05). Empty = not testable in that cell type. '
          'Other sheets: one sheet per comparison (dataset labels as in Table S3). PSI = LeafCutter fitted junction usage per group; p.adjust is per event (cluster). SE = '
          'skipped exon; A3SS / A5SS = alternative 3\' / 5\' splice site; well-supported splice-site correspondence: '
          'identical distance between the alternative sites (A3SS/A5SS) or identical skipped-exon length (SE) in '
          'human and mouse (empty for complex events). Exon ranks: rank, in the Ensembl canonical transcript, of the exon '
          'whose boundary is at the junction start and end coordinates (start = lower genomic coordinate, so for '
          'minus-strand genes the start is the downstream exon); Last = last exon; * = exon not in the canonical '
          'transcript.', sheets)


def pair_summary(ann):
    """Per junction: dataset PSI of each cell type (as in the five-dataset classification), the human-mouse gap of
    the junction and the status of the event in the cell type."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import classify_separation_116 as cs
    counts = pd.read_csv(cs.BASE + 'events/id50/AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(cs.BASE + 'events/id50/leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t', usecols=['cell', 'cluster', 'status_high_low'])
    st = st.pivot(index='cluster', columns='cell', values='status_high_low').replace(
        {'divergent': 'DS', 'conserved': 'equivalent'})
    cols = []
    for cell in CELLS:
        psi, _ = cs.dataset_psi(counts, design, cell)
        part = psi[cs.DATASETS].add_prefix(f'{cell} ')
        part[f'{cell} gap'] = cs.gap(psi, cs.HUMAN)
        cols.append(part)
    j = pd.concat(cols, axis=1)
    j['cluster'] = j.index.str.split(':').str[3]
    for cell in CELLS:
        j.insert(j.columns.get_loc(f'{cell} gap') + 1, f'{cell} status', j['cluster'].map(st[cell]))
    info = ann.set_index('cluster').loc[j['cluster']]
    parts = j.index.to_series().str.split(':', expand=True)
    first = pd.DataFrame({'Alternative splicing event ID': parts[0] + ':' + parts[3],
                          'Human gene symbol': info['symbol'].to_numpy(),
                          'Alternative splicing type': info['event_type'].to_numpy(),
                          'Human junction': parts[0] + ':' + parts[1] + '-' + parts[2]}, index=j.index)
    n_ds = (st == 'DS').sum(axis=1)
    out = pd.concat([first, j.drop(columns='cluster').round(3)], axis=1)
    out['_n'], out['_s'] = j['cluster'].map(n_ds).to_numpy(), parts[1].astype(int).to_numpy()
    out = out.sort_values(['_n', 'Human gene symbol', 'Alternative splicing event ID', '_s'],
                          ascending=[False, True, True, True]).drop(columns=['_n', '_s'])
    return out.reset_index(drop=True)


def table_s7():
    s = pd.read_csv(ED + 'effects/comparison_summary.tsv', sep='\t')
    s['kind'] = s['comparison'].map(kind)
    # matched pairs were replaced by the five-dataset classification (Table S5)
    s = s[~s['comparison'].str.startswith('matched_')].copy()
    s = s.rename(columns={'testable': 'testable events', 'pct_DS': '% differentially spliced',
                          'pct_equivalent_m0.025': '% equivalent (+-0.025)', 'pct_equivalent_m0.05': '% equivalent (+-0.05)',
                          'pct_equivalent_m0.1': '% equivalent (+-0.1)', 'median_max_abs_dpsi': 'median event max|dPSI|',
                          'pearson': 'Pearson r of junction PSI'})
    write('Table_S8', 'Supplementary Table S8: Differential splicing output of all comparisons',
          'Differentially spliced: p.adjust <= 0.05 and max|dPSI| >= 0.1; equivalent: 90% confidence interval of every '
          'junction dPSI within the margin (two one-sided tests at 0.05); testable = enough covered samples. '
          'Human-mouse pairs (Figure 2C) are a step before the classification of the events from all five datasets '
          '(Table 1).',
          {'comparisons': s.round(2)})


def table_s8():
    s = pd.read_csv(ED + 'classification/summary_by_cell.tsv', sep='\t', index_col=0)
    s.columns = [c.replace('high_low_', 'high + low: ').replace('high_', 'high: ').replace('divergent', 'DS')
                 .replace('conserved', 'equivalent').replace('_DS_within_species', ', DS within species')
                 for c in s.columns]
    write('Table_S9', 'Supplementary Table S9: Differentially spliced and equivalent events by confidence level',
          'Testable: every junction of the event has a PSI (mean over the samples with >= 20 cluster reads) in each of '
          'the three human and two mouse datasets. Differentially spliced (DS): for at least one junction, the PSI of '
          'every human dataset differs from the PSI of every mouse dataset by at least the gap (high confidence: gap '
          '>= 0.2; high + low: gap >= 0.1). Equivalent: the PSI of every junction lies within a range of at most 0.1 '
          'over all five datasets (high confidence: at most 0.05). high + low includes the high-confidence events; '
          '"DS within species" = events of the class that are differentially spliced (p.adjust <= 0.05 and max|dPSI| '
          '>= 0.1) in at least one within-species pair; "DS_expected" = mean number of DS events over the nine '
          'splits of the five datasets into three and two that mix the species (expected by chance).',
          {'counts': s.reset_index()})


def table_s9(status, desc, ann):
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    fib = ev[ev['comparison'] == 'FH_FM'].set_index('cluster')['status_m0.1']
    wide = status.pivot_table(index='cluster', columns='cell', values='status_confidence', aggfunc='first')
    wide = wide[[c for c in CELLS if c in wide.columns]].replace(
        {'divergent high': 'DS, high', 'divergent low': 'DS, low', 'conserved high': 'equivalent, high',
         'conserved low': 'equivalent, low'})
    info = status.drop_duplicates('cluster').set_index('cluster')[['symbol', 'ensembl', 'event_type']]
    t = info.join(wide)
    for col, lab in (('species_gap', 'gap'), ('psi_range', 'range')):
        v = status.pivot_table(index='cluster', columns='cell', values=col, aggfunc='first').round(3)
        t = t.join(v[[c for c in CELLS if c in v.columns]].add_suffix(f' {lab}'))
    t['Fibroblasts (human vs mouse)'] = t.index.map(fib).fillna('not testable')
    t['Gene description (HGNC)'] = t['ensembl'].map(desc)
    t.index = t.index.map(ann.set_index('cluster')['chrom']) + ':' + t.index
    t = t.reset_index().rename(columns={'cluster': 'Alternative splicing event', 'symbol': 'Gene',
                                        'ensembl': 'ENSEMBL ID', 'event_type': 'Type'})
    write('Table_S10', 'Supplementary Table S10: Classification of each orthologous alternative splicing event',
          'Per immune cell type: differentially spliced (DS) or equivalent with high or low confidence, inconclusive, or not testable '
          '(definitions in Table S9); gap = smallest PSI difference between a human and a mouse dataset (negative: the '
          'human and mouse PSI ranges overlap), range = PSI range over the five datasets (max over junctions). Fibroblasts: single comparison, no confidence level (DS / equivalent at +-0.1 / '
          'inconclusive).',
          {'events': t})


def relabel(df):
    """Class names as in the manuscript: divergent -> DS (differentially spliced), conserved -> equivalent."""
    def r(x):
        return (x.replace('species-divergent', 'DS').replace('divergent', 'DS').replace('conserved', 'equivalent')
                if isinstance(x, str) else x)
    df = df.copy()
    df.columns = [r(c) for c in df.columns]
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(r)
    return df


def table_s10():
    sheets = {}
    for name in ('conserved_only', 'divergent_only'):
        for lib in ('Reactome_2022', 'GO_Biological_Process_2025', 'hallmark'):
            f = ED + f'revision/enrichment_{name}_{lib}.tsv'
            d = pd.read_csv(f, sep='\t')
            label = {'conserved_only': 'equivalent_only', 'divergent_only': 'DS_only'}[name]
            sheets[f'{label}_{lib}'[:31]] = d[d['Adjusted P-value'] <= 0.1] if len(d) else d
    sheets['coverage_check'] = relabel(pd.read_csv(ED + 'revision/enrichment_coverage_check.tsv', sep='\t'))
    write('Table_S11', 'Supplementary Table S11: Functional enrichment of high-confidence genes',
          'MSigDB 2024.1.Hs; background = genes with at least one testable event; terms with FDR <= 0.1. '
          'coverage_check: logistic regression of list membership on term membership adjusted for event coverage.',
          sheets)


def table_s11():
    sheets = {}
    rows = []
    for v in ['id50', 'id70', 'id80', 'id50_rbh']:
        s = pd.read_csv(BASE + f'events/{v}/effects/comparison_type_summary.tsv', sep='\t')
        c = pd.read_csv(BASE + f'events/{v}_sep/classification/summary_by_cell.tsv', sep='\t')
        s.insert(0, 'orthology', v)
        rows.append(s)
        c.insert(0, 'orthology', v)
        sheets.setdefault('_hc', []).append(c)
    # human-mouse events are classified only from all five datasets; single human-mouse pairs give effect sizes only
    sheets = {'orthology_cutoff_high_conf': pd.concat(sheets['_hc'])}
    d = pd.read_csv(ED + 'revision/dpsi_cutoff_sensitivity.tsv', sep='\t')
    sheets['dPSI_cutoff_within_species'] = d[d['kind'].isin(['HH', 'MM', 'fibro'])]
    ce = pd.read_csv(ED + 'revision/common_events_summary.tsv', sep='\t')
    sheets['common_events'] = ce[['cell', 'comparison', 'kind', 'common_events', 'median_max_abs_dpsi', 'pearson']]
    sheets['common_events_tests'] = pd.read_csv(ED + 'revision/common_events_tests.tsv', sep='\t')
    sheets['splice_site_subset'] = pd.read_csv(ED + 'revision/splice_site_subset_summary.tsv', sep='\t') \
        .drop(columns=['pct_cross_DS', 'pct_within_DS'])
    sheets['minor_isoform_usage'] = pd.read_csv(ED + 'revision/minor_usage_summary.tsv', sep='\t')
    sheets['small_shift_A3SS_A5SS'] = pd.read_csv(ED + 'revision/small_shift_test.tsv', sep='\t')
    b = pd.read_csv(ED + 'revision/single_cell_type_divergent_breakdown.tsv', sep='\t')
    sheets['single_cell_type_DS'] = b[b['measure'] != 'human-mouse pairs in other cell types']
    sheets['matched_populations'] = pd.read_csv(ED + 'revision/matched_separation.tsv', sep='\t')
    sheets = {k: relabel(v) for k, v in sheets.items()}
    write('Table_S5', 'Supplementary Table S5: Sensitivity analyses',
          'Orthology stringency (id50 main, id70, id80, reciprocal best), dPSI cutoff for differential splicing within '
          'species, effect sizes on the events testable in every comparison of a cell type, events with '
          'well-supported splice-site correspondence, the classification with mouse M1 restricted to naive T cells '
          'and blood monocytes (matched; all = all M1 samples), classification by '
          'minor-isoform usage, small (<= 6 nt) alternative splice-site shifts, and the status in the other cell types '
          'of events that are differentially spliced with high confidence in only one cell type.', sheets)


def table_s12():
    write('Table_S6', 'Supplementary Table S6: Effect of coverage on testability and classification',
          'Logistic regressions over event x cell type (cell type as covariate); coverage = log10 mean reads per '
          'sample of the event across immune samples.',
          {'model': pd.read_csv(ED + 'revision/coverage_bias_model.tsv', sep='\t')})


def main():
    os.makedirs(OUT, exist_ok=True)
    ann = pd.read_csv(ED + 'classification/events_annotation.tsv', sep='\t')
    status = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
    points = pd.read_csv(BASE + f'orthology/unique_points_116_{V}.txt', sep='\t', index_col=0)
    table_s4()
    table_s5(points, ann)
    table_s7()
    table_s8()
    table_s9(status, descriptions(), ann)
    table_s10()
    table_s11()
    table_s12()


if __name__ == '__main__':
    main()
