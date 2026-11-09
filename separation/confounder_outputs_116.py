#!/usr/bin/env python3
"""
Figure S7 and Table S12: sex, age and tissue control comparisons vs the human-mouse and within-species
comparisons of the same cell type, on the same events (confounders_116.py summary).
Outputs: revision_separation/figures/Figure_S7.pdf/.tif, revision_separation/supplementary/Table_S12.xlsx
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from build_supp_tables_116 import times_new_roman

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50_sep/'
FIG = B + 'revision_separation/figures/'
SUP = B + 'revision_separation/supplementary/'
matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
CONTROLS = [  # (comparison, label, dataset)
    ('human_sex', 'Human sex', 'GSE107011'), ('human_age', 'Human age', 'GSE143213'),
    ('human_tissue', 'Human tissue', 'GSE133383'), ('mouse_sex', 'Mouse sex', 'ENCODE'),
    ('mouse_age', 'Mouse age', 'GSE184496'), ('mouse_tissue', 'Mouse tissue', 'GSE233347'),
    ('mouse_tissue_m1', 'Mouse tissue', 'GSE116177')]
GROUPS = {'human_sex': 'male vs female', 'human_age': 'young vs old', 'human_tissue': 'blood vs spleen',
          'mouse_sex': 'male vs female', 'mouse_age': 'young vs old', 'mouse_tissue': 'spleen vs liver',
          'mouse_tissue_m1': 'bone marrow vs blood'}
COLORS = {'control': '#7b3294', 'HH': '#2a78d6', 'MM': '#1baf7a', 'cross': '#eb6834'}
LABELS = {'control': 'Control (sex / age / tissue): two sample groups, same species and same dataset', 'HH': 'Human-Human dataset pairs',
          'MM': 'Mouse-Mouse dataset pairs', 'cross': 'Human-Mouse dataset pairs'}


def main():
    s = pd.read_csv(ED + 'revision/confounders_summary.tsv', sep='\t')
    ctrl = pd.read_csv(ED + 'effects_confounders/event_effects.tsv', sep='\t')
    spec = pd.read_csv(B + 'metadata/confounder_comparisons.tsv', sep='\t')

    # ---------- Table S12
    tests = s[s['kind'] == 'test'].set_index(['control', 'cell'])['wilcoxon_p_cross_gt_control']
    rows = []
    for comp, label, ds in CONTROLS:
        for cell, g in s[(s['control'] == comp) & (s['kind'] != 'test')].groupby('cell'):
            c = ctrl[(ctrl['comparison'] == comp) & (ctrl['cell'] == cell)]
            grp = spec[(spec['comparison'] == comp) & (spec['cell'] == cell)].groupby('group')['sample'] \
                .apply(lambda x: ', '.join(v.rsplit('_', 1)[0] for v in x))
            row = dict(comparison=f'{label} ({GROUPS[comp]})', dataset=ds, cell_type=cell,
                       group_samples='; '.join(f'{k}: {v}' for k, v in grp.items()),
                       testable_events=int((c['status_m0.05'] != 'not_testable').sum()),
                       DS_events=int((c['status_m0.05'] == 'DS').sum()),
                       common_events=int(g['common_events'].iloc[0]))
            for k in ('control', 'cross', 'HH', 'MM'):
                r = g[g['kind'] == k]
                if len(r):
                    row[f'pct_DS_{k}'] = round(r['pct_DS_median'].iloc[0], 2)
                    row[f'median_abs_dpsi_{k}'] = round(r['median_max_abs_dpsi'].iloc[0], 3)
                    row[f'pearson_{k}'] = round(r['pearson_median'].iloc[0], 3)
            row['wilcoxon_p_human_mouse_gt_control'] = tests.get((comp, cell))
            rows.append(row)
    t = pd.DataFrame(rows)
    with pd.ExcelWriter(SUP + 'Table_S12.xlsx', engine='openpyxl') as xw:
        t.to_excel(xw, sheet_name='controls', index=False, startrow=2)
        ws = xw.sheets['controls']
        ws['A1'] = 'Supplementary Table S12: Sex, age and tissue control comparisons'
        ws['A2'] = ('Each control varies one factor within a species (same cell type, laboratory and protocol). '
                    'Comparisons with human-mouse (cross), human-human (HH) and mouse-mouse (MM) dataset pairs use '
                    'the events testable in the control and in every human-mouse pair of the cell type (common '
                    'events); values for HH, MM and cross are medians over dataset pairs. DS: p.adjust <= 0.05 and '
                    'max|dPSI| >= 0.1. Wilcoxon: paired test that the mean human-mouse max|dPSI| is larger than in '
                    'the control.')
        ws.freeze_panes = 'A4'
    times_new_roman(SUP + 'Table_S12.xlsx')
    print('wrote Table_S12.xlsx')

    # ---------- Figure S7
    # A: dataset pairs as in Table 1 (median % DS over the pairs, all events testable in each pair)
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ev = ev[ev['status_m0.1'] != 'not_testable']
    a, b = ev['comparison'].str[0], ev['comparison'].str.split('_').str[1].str[0]
    ev['kind'] = np.where(a != b, 'cross', np.where(a == 'H', 'HH', 'MM'))
    ev = ev[~ev['comparison'].str.startswith('matched_') & (ev['comparison'] != 'FH_FM')]
    table1 = (ev.groupby(['kind', 'cell', 'comparison'])['status_m0.1'].apply(lambda x: 100 * (x == 'DS').mean())
              .groupby(['kind', 'cell']).median())
    main_ctrl = CONTROLS
    labels, species, data = [], [], {k: [] for k in COLORS}
    for comp, label, ds in main_ctrl:
        for cell, g in s[(s['control'] == comp) & (s['kind'] != 'test')].groupby('cell'):
            labels.append(f"{label.split()[1]}\n{CELL_LABEL[cell]}\n{GROUPS[comp].replace(' vs ', '/').replace('bone marrow', 'BM')}\n{ds}")
            species.append(label.split()[0])
            for k in COLORS:
                r = g[g['kind'] == k]
                v = (r[['pct_DS_median', 'median_max_abs_dpsi', 'pearson_median']].iloc[0].values.astype(float)
                     if len(r) else np.full(3, np.nan))
                if k != 'control':
                    v[0] = table1.get((k, cell), np.nan)
                data[k].append(v)
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True, gridspec_kw={'hspace': 0.15})
    x = np.arange(len(labels))
    w = 0.2
    for j, (ylab, lim) in enumerate([('% differentially spliced', None), ('median event max |ΔPSI|', None),
                                     ('Pearson r of junction PSI', (0.75, 1.0))]):
        ax = axes[j]
        for i, k in enumerate(COLORS):
            vals = [v[j] for v in data[k]]
            if j == 0 and k == 'cross':      # human-mouse events are classified only from all five datasets
                continue
            ax.bar(x + (i - 1.5) * w, vals, width=w, color=COLORS[k], label=LABELS[k])
        ax.set_ylabel(ylab)
        if lim:
            ax.set_ylim(*lim)
        ax.text(-0.08, 1.02, 'ABC'[j], transform=ax.transAxes, fontsize=10, fontweight='bold')
    from matplotlib.patches import Patch
    axes[0].legend([Patch(color=COLORS[k]) for k in COLORS],
                   [LABELS[k] + (' (B, C)' if k == 'cross' else '') for k in COLORS],
                   frameon=False, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.09), fontsize=6.5)
    split = species.index('Mouse') - 0.5
    for ax in axes:
        ax.axvline(split, color='#9e9e9e', lw=0.8, ls='--')
    for xs, name in ((-0.5 + (split + 0.5) / 2, 'Human controls: within one species and one dataset'),
                     ((split + len(labels) - 0.5) / 2, 'Mouse controls: within one species and one dataset')):
        axes[0].text(xs, 1.0, name, transform=axes[0].get_xaxis_transform(), ha='center', va='bottom',
                     fontsize=7, fontweight='bold')
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, fontsize=5.5)
    fig.text(0.02, 0.0, 'Figure S7', fontsize=8)
    fig.savefig(FIG + 'Figure_S7.pdf', bbox_inches='tight')
    fig.savefig(FIG + 'Figure_S7.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    print('wrote Figure_S7')


if __name__ == '__main__':
    main()
