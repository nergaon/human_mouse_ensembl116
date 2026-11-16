#!/usr/bin/env python3
"""
Supplementary tables (Excel) for the revised manuscript, from the Ensembl 116 / id50 results
(orthology sensitivity tables also use id70, id80 and id50_rbh).

  Table_S1  orthology map per identity cutoff
  Table_S8  junction-level LeafCutter results, one sheet per comparison
  Table_S9  differential splicing output of every comparison x cell type
  Table_S10  high- and low-confidence counts
  Table_S11  event status per cell type (confidence level) and in fibroblasts
  Table_S12 enrichment of high-confidence DS-only / equivalent-only genes, with coverage check
  Table_S6 sensitivity analyses (orthology cutoff, dPSI cutoff, common events, splice-site
            correspondence)
  Table_S7 coverage bias model
Output: ensembl116/revision_6ds/supplementary/
"""
import itertools
import os
import re
import shutil
import sys
import zipfile

import numpy as np
import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
V = 'id50_6ds'
ED = BASE + f'events/{V}_anova/'
OUT = BASE + 'revision_6ds/supplementary/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
DATASETS = {'H1': 'GSE60424', 'H2': 'GSE115736', 'H3': 'GSE118165', 'M1': 'GSE287527 / GSE215775 / GSE313371 / GSE261448',
            'M2': 'GSE307667 / GSE184496 / GSE263132 / GSE236819 / GSE254202',
            'M3': 'GSE231772 / GSE134322 / GSE284441 / GSE180053 / PRJNA797236', 'FH': 'E-MTAB-5919 human',
            'FM': 'E-MTAB-5919 mouse'}


def write(name, title, note, sheets):
    """sheets: {sheet name: DataFrame}; row 1 title, row 2 note, header in row 3. Internal table names are
    renumbered for publication (table_numbers_116); the removed coverage table (internal S7) is not written."""
    import table_numbers_116 as tn
    if tn.name(name) is None:
        print('not written (removed table):', name)
        return
    name, title, note = tn.name(name), tn.renum(title), tn.renum(note)
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
        a = pd.read_csv(BASE + f'events/{v}_6ds_anova/classification/events_annotation.tsv', sep='\t')
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
           'rank_consistency': pd.concat([pd.read_csv(BASE + f'orthology/rank_consistency_116_{v}_6ds.tsv', sep='\t')
                                          for v in ['id50', 'id70', 'id80', 'id50_rbh']])})


def table_s5(points, ann):
    info = ann.set_index('cluster')
    h2m = dict(zip(points['position_h'], points['position_m']))
    rank_h = dict(zip(points['position_h'], points['rank_h']))
    rank_m = dict(zip(points['position_h'], points['rank_m']))
    sheets = {}
    # only the pairs of the article (matched_ pairwise LeafCutter runs are not reported; the matched populations are
    # analysed by the five-dataset classification, Table S6)
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
    # ANOVA version: no split-based summary sheet (per-event classification in Table S11)
    write('Table_S8', 'Supplementary Table S8: Junction-level LeafCutter differential splicing results',
          'One sheet per comparison (dataset labels as in Table S4). PSI = LeafCutter fitted junction usage per group; p.adjust is per event (cluster). SE = '
          'skipped exon; A3SS / A5SS = alternative 3\' / 5\' splice site (both sites on one exon); AFE / ALE = '
          'alternative first / last exon; well-supported splice-site correspondence: '
          'identical distance between the alternative sites (A3SS/A5SS) or identical skipped-exon length (SE) in '
          'human and mouse (empty for AFE, ALE and complex events). Exon ranks: rank, in the Ensembl canonical transcript, of the exon '
          'whose boundary is at the junction start and end coordinates (start = lower genomic coordinate, so for '
          'minus-strand genes the start is the downstream exon); Last = last exon; * = exon not in the canonical '
          'transcript.', sheets)


def pair_summary(ann):
    """Per junction: dataset PSI of each cell type (as in the five-dataset classification), the human-mouse gap of
    the junction and the status of the event in the cell type."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import classify_separation_116 as cs
    counts = pd.read_csv(cs.BASE + 'events/id50_6ds/AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(cs.BASE + 'events/id50_6ds/leafcutter/sample_design.tsv', sep='\t').set_index('sample')
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
        # the nine splits that mix the species: event gap (max over junctions) of each, listed when >= 0.1
        psi, clu = cs.dataset_psi(counts, design, cell)
        g = pd.DataFrame({s: cs.gap(psi, s).groupby(clu).max() for s in itertools.combinations(cs.DATASETS, 3)
                          if s != cs.HUMAN})
        name = {s: ' '.join(s) + ' vs ' + ' '.join(d for d in cs.DATASETS if d not in s) for s in g.columns}
        n_low, n_high = g.ge(cs.GAP_LOW).sum(axis=1), g.ge(cs.GAP_HIGH).sum(axis=1)
        listed = g.apply(lambda r: '; '.join(f'{name[s]} ({v:.2f})' for s, v in r.sort_values(ascending=False).items()
                                             if v >= cs.GAP_LOW), axis=1)
        k = j.columns.get_loc(f'{cell} status') + 1
        j.insert(k, f'{cell} mixed splits gap>=0.1', j['cluster'].map(n_low))
        j.insert(k + 1, f'{cell} mixed splits gap>=0.2', j['cluster'].map(n_high))
        j.insert(k + 2, f'{cell} mixed splits separating the event', j['cluster'].map(listed))
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
    # matched pairs were replaced by the five-dataset classification (Table S6)
    s = s[~s['comparison'].str.startswith('matched_')].copy()
    s = s.drop(columns=['pct_equivalent_m0.025', 'pct_equivalent_m0.05'])   # same status at every margin
    s = s.rename(columns={'testable': 'testable events', 'pct_DS': '% differentially spliced',
                          'pct_equivalent_m0.1': '% splicing conserved', 'median_max_abs_dpsi': 'median event max|dPSI|',
                          'pearson': 'Pearson r of junction PSI'})
    write('Table_S9', 'Supplementary Table S9: Differential splicing output of all comparisons',
          'Differentially spliced: LeafCutter p.adjust <= 0.05 and max|dPSI| >= 0.1; splicing conserved: LeafCutter P '
          '(not adjusted) > 0.1 and max|dPSI| <= 0.05; inconclusive: other testable events; testable = enough covered '
          'samples. '
          'Human-mouse pairs (Figure 2C) are a step before the classification of the events from all six datasets '
          '(Table 1).',
          {'comparisons': s.round(2)})


def table_s8():
    s = pd.read_csv(ED + 'classification/summary_by_cell.tsv', sep='\t', index_col=0)
    keep = ['high_not testable', 'high_divergent', 'high_conserved', 'high_inconclusive',
            'high_low_divergent_DS_within_species', 'high_low_conserved_DS_within_species']
    s = s[keep]
    s.insert(0, 'testable', s[['high_divergent', 'high_conserved', 'high_inconclusive']].sum(axis=1))
    s.columns = ['testable', 'not testable', 'DS', 'splicing conserved', 'inconclusive',
                 'DS, also DS within species', 'splicing conserved, also DS within species']
    write('Table_S10', 'Supplementary Table S10: Human-mouse classification of the alternative splicing events',
          'Per immune cell type. Testable: a PSI in at least two samples of each of the six datasets. ANOVA on the '
          'junction with the largest range of dataset PSIs. Splicing conserved: ANOVA P > 0.1 and dataset PSIs within '
          'a range of 0.05 in every junction. Differentially spliced (DS): ANOVA Benjamini-Hochberg adjusted P <= 0.05 '
          'and, in at least one junction, |mean of the human dataset PSIs - mean of the mouse dataset PSIs| >= 0.1 and '
          '1.5 x the largest within-species |dPSI| < the smallest human-mouse |dPSI| (dataset PSIs). Inconclusive: all '
          'other testable events. "also DS '
          'within species": differentially spliced (LeafCutter adjusted P <= 0.05 and max|dPSI| >= 0.1) in at least '
          'one within-species dataset pair.', {'counts': s.reset_index()})


def table_s9(status, desc, ann):
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    fib = ev[ev['comparison'] == 'FH_FM'].set_index('cluster')['status_m0.1']
    tests = pd.read_csv(ED + 'classification/event_tests.tsv', sep='\t')
    wide = status.pivot_table(index='cluster', columns='cell', values='status_high', aggfunc='first')
    wide = wide[[c for c in CELLS if c in wide.columns]].replace({'divergent': 'DS', 'conserved': 'splicing conserved'})
    info = status.drop_duplicates('cluster').set_index('cluster')[['symbol', 'ensembl', 'event_type']]
    t = info.join(wide)
    for col, lab in (('junction', 'ANOVA junction'), ('anova_p', 'ANOVA P'), ('anova_p_adjust', 'ANOVA adjusted P'),
                     ('species_junction', 'species junction'), ('gap', 'human-mouse gap'),
                     ('max_within_dpsi', 'largest within-species |dPSI|'),
                     ('min_HM_dpsi', 'smallest human-mouse |dPSI|')):
        v = tests.pivot_table(index='cluster', columns='cell', values=col, aggfunc='first')
        if not col.endswith('junction'):
            v = v.astype(float).map(lambda x: float(f'{x:.3g}') if pd.notna(x) else x)
        else:
            v = v.map(lambda x: ':'.join(x.split(':')[:3]) if isinstance(x, str) else x)
        t = t.join(v[[c for c in CELLS if c in v.columns]].add_suffix(f' {lab}'))
    t['Fibroblasts (human vs mouse)'] = pd.Series(t.index.map(fib), index=t.index).fillna('not testable').replace(
        {'equivalent': 'splicing conserved', 'not_testable': 'not testable'})
    t['Gene description (HGNC)'] = t['ensembl'].map(desc)
    t.index = t.index.map(ann.set_index('cluster')['chrom']) + ':' + t.index
    t = t.reset_index().rename(columns={'cluster': 'Alternative splicing event', 'symbol': 'Gene',
                                        'ensembl': 'ENSEMBL ID', 'event_type': 'Type'})
    write('Table_S11', 'Supplementary Table S11: Classification of each orthologous alternative splicing event',
          'Per immune cell type: differentially spliced (DS), splicing conserved, inconclusive or not testable '
          '(definitions in Table S10); the ANOVA junction (largest range of dataset PSIs) with the ANOVA P and adjusted '
          'P; the species junction (a junction that meets the DS rule, with the largest human-mouse gap, else the '
          'junction with the largest gap) with its gap (|mean human - mean mouse dataset PSI|), the largest |dPSI| over '
          'the six within-species dataset pairs and the smallest over the nine human-mouse pairs. Fibroblasts: single '
          'comparison (DS / splicing conserved / inconclusive, as in Table S9).',
          {'events': t})


def relabel(df):
    """Class names as in the manuscript: divergent -> DS (differentially spliced), conserved -> equivalent."""
    def r(x):
        return (x.replace('species-divergent', 'DS').replace('divergent', 'DS').replace('conserved', 'splicing conserved')
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
            label = {'conserved_only': 'conserved_only', 'divergent_only': 'DS_only'}[name]
            sheets[f'{label}_{lib}'[:31]] = d[d['Adjusted P-value'] <= 0.1] if len(d) else d
    write('Table_S12', 'Supplementary Table S12: Functional enrichment of high-confidence genes',
          'MSigDB 2024.1.Hs; background = genes with at least one testable event; terms with FDR <= 0.1. ',
          sheets)


def table_s11():
    """Table S6 (ANOVA version): one class per event (no confidence levels), so only one set of columns."""
    sheets = {}
    orth = []
    for v in ['id50', 'id70', 'id80', 'id50_rbh']:
        c = pd.read_csv(BASE + f'events/{v}_6ds_anova/classification/summary_by_cell.tsv', sep='\t', index_col=0)
        st = pd.read_csv(BASE + f'events/{v}_6ds_anova/classification/event_status_by_cell.tsv', sep='\t')
        t = c['high_divergent'] + c['high_conserved'] + c['high_inconclusive']
        o = pd.DataFrame({'orthology': v, 'testable': t, 'DS': c['high_divergent'],
                          'pct_DS': (100 * c['high_divergent'] / t).round(1), 'splicing conserved': c['high_conserved'],
                          'pct_splicing_conserved': (100 * c['high_conserved'] / t).round(1),
                          'inconclusive': c['high_inconclusive'], 'not testable': c['high_not testable']})
        o['DS events (any cell type)'] = st.loc[st['status_high'] == 'divergent', 'cluster'].nunique()
        o['splicing conserved events (any cell type)'] = st.loc[st['status_high'] == 'conserved', 'cluster'].nunique()
        orth.append(o.rename_axis('cell').reset_index())
    sheets['orthology_cutoff'] = pd.concat(orth)[['orthology', 'cell', 'testable', 'DS', 'pct_DS', 'splicing conserved',
                                                  'pct_splicing_conserved', 'inconclusive', 'not testable',
                                                  'DS events (any cell type)', 'splicing conserved events (any cell type)']]
    d = pd.read_csv(ED + 'revision/dpsi_cutoff_sensitivity.tsv', sep='\t')
    sheets['dPSI_cutoff_within_species'] = d[d['kind'].isin(['HH', 'MM', 'fibro'])]
    ce = pd.read_csv(ED + 'revision/common_events_summary.tsv', sep='\t')
    sheets['common_events'] = ce[['cell', 'comparison', 'kind', 'common_events', 'median_max_abs_dpsi', 'pearson']]
    sheets['common_events_tests'] = pd.read_csv(ED + 'revision/common_events_tests.tsv', sep='\t')
    ss = pd.read_csv(ED + 'revision/splice_site_subset_summary.tsv', sep='\t')
    sheets['splice_site_subset'] = ss[['subset', 'cell', 'testable', 'pct_high_divergent', 'pct_high_conserved']] \
        .rename(columns={'pct_high_divergent': 'pct_DS', 'pct_high_conserved': 'pct_splicing_conserved'})
    mu = pd.read_csv(ED + 'revision/minor_usage_summary.tsv', sep='\t')
    sheets['minor_isoform_usage'] = mu[mu['definition'] == 'high'].drop(columns='definition')
    sh = pd.read_csv(ED + 'revision/small_shift_test.tsv', sep='\t')
    sheets['small_shift_A3SS_A5SS'] = sh[sh['definition'] == 'high'].drop(columns='definition') \
        .rename(columns=lambda c: c.replace('_divergent', '_DS'))
    b = pd.read_csv(ED + 'revision/single_cell_type_divergent_breakdown.tsv', sep='\t')
    b = b[(b['definition'] == 'high') & (b['measure'] == 'high-confidence status in other cell types')]
    sheets['single_cell_type_DS'] = b.drop(columns=['definition', 'measure']).rename(
        columns={'divergent_events': 'DS_events', 'events_high_confidence_in_one_cell_type': 'DS_in_one_cell_type',
                 'event_x_other_cell_type': 'event_x_other_cell_type', 'category': 'status_in_other_cell_types'})
    # the classification of every split of the six datasets into two groups of three (the human-mouse split and the
    # nine splits that mix the species): DSA option 1 (dsa_analysis_116.py, events/id50_6ds)
    sp = pd.read_csv(BASE + 'revision_6ds/conserved_options_DS_split_calibration.tsv', sep='\t')
    sheets['species_vs_mixed_splits'] = sp.rename(columns={'split': 'split (DS events per cell type)'})
    sheets = {k: relabel(v) for k, v in sheets.items()}
    write('Table_S6', 'Supplementary Table S6: Sensitivity analyses',
          'orthology_cutoff: the human-mouse classification (DS, splicing conserved, inconclusive; Methods 2.5) with orthologous '
          'exons defined by 50% (main), 70% or 80% identity, or 50% identity and reciprocal best matches. '
          'dPSI_cutoff_within_species: % differentially spliced in single dataset pairs at |dPSI| thresholds 0.05, 0.1, '
          '0.2. common_events(_tests): effect sizes on the events testable in all 15 dataset pairs of a cell type, and '
          'the paired Wilcoxon test that human-mouse max|dPSI| > within-species. splice_site_subset: the classification '
          'restricted to events with well-supported splice-site correspondence. minor_isoform_usage: classification by '
          'the usage of the minor isoform. small_shift_A3SS_A5SS: % DS for alternative splice sites up to 6 nt apart vs '
          'more distant (Fisher test). single_cell_type_DS: status in the other cell types of the events that are '
          'differentially spliced in only one cell type. species_vs_mixed_splits: the classification steps applied to '
          'each split of the six datasets into two groups of three (ANOVA adjusted P <= 0.05 and, in a junction, gap '
          '>= 0.1 and 1.5 x the largest within-group |dPSI| < the smallest between-group |dPSI|); the human-mouse split and the nine splits that mix the species. There are no confidence levels.', sheets)


def table_s12():
    write('Table_S7', 'Supplementary Table S7: Effect of coverage on testability and classification',
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
