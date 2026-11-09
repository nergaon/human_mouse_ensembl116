#!/usr/bin/env python3
"""
Figure S7 and Table S13: sex, age and tissue control comparisons vs the human-mouse and within-species
comparisons of the same cell type, on the same events (confounders_116.py summary).
Outputs: revision_files/figures/Figure_S7.pdf/.tif, revision_files/supplementary/Table_S13.xlsx
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from build_supp_tables_116 import times_new_roman

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
FIG = B + 'revision_files/figures/'
SUP = B + 'revision_files/supplementary/'
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
LABELS = {'control': 'Control (sex / age / tissue)', 'HH': 'Human-Human datasets', 'MM': 'Mouse-Mouse datasets',
          'cross': 'Human-Mouse datasets'}


def main():
    s = pd.read_csv(ED + 'revision/confounders_summary.tsv', sep='\t')
    ctrl = pd.read_csv(ED + 'effects_confounders/event_effects.tsv', sep='\t')
    spec = pd.read_csv(B + 'metadata/confounder_comparisons.tsv', sep='\t')

    # ---------- Table S13
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
    with pd.ExcelWriter(SUP + 'Table_S13.xlsx', engine='openpyxl') as xw:
        t.to_excel(xw, sheet_name='controls', index=False, startrow=2)
        ws = xw.sheets['controls']
        ws['A1'] = 'Supplementary Table S13: Sex, age and tissue control comparisons'
        ws['A2'] = ('Each control varies one factor within a species (same cell type, laboratory and protocol). '
                    'Comparisons with human-mouse (cross), human-human (HH) and mouse-mouse (MM) dataset pairs use '
                    'the events testable in the control and in every human-mouse pair of the cell type (common '
                    'events); values for HH, MM and cross are medians over dataset pairs. DS: p.adjust <= 0.05 and '
                    'max|dPSI| >= 0.1. Wilcoxon: paired test that the mean human-mouse max|dPSI| is larger than in '
                    'the control.')
        ws.freeze_panes = 'A4'
    times_new_roman(SUP + 'Table_S13.xlsx')
    print('wrote Table_S13.xlsx')

    # ---------- Figure S7
    main_ctrl = CONTROLS
    labels, data = [], {k: [] for k in COLORS}
    for comp, label, ds in main_ctrl:
        for cell, g in s[(s['control'] == comp) & (s['kind'] != 'test')].groupby('cell'):
            labels.append(f"{label.replace('Human', 'H').replace('Mouse', 'M')}\n{CELL_LABEL[cell]}" +
                          ('\n(Mo BM/PB)' if comp == 'mouse_tissue_m1' else ''))
            for k in COLORS:
                r = g[g['kind'] == k]
                data[k].append(r[['pct_DS_median', 'median_max_abs_dpsi', 'pearson_median']].iloc[0].values
                               if len(r) else [np.nan] * 3)
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.6), sharex=True, gridspec_kw={'hspace': 0.15})
    x = np.arange(len(labels))
    w = 0.2
    for j, (ylab, lim) in enumerate([('% differentially spliced', None), ('median event max |ΔPSI|', None),
                                     ('Pearson r of junction PSI', (0.75, 1.0))]):
        ax = axes[j]
        for i, k in enumerate(COLORS):
            vals = [v[j] for v in data[k]]
            ax.bar(x + (i - 1.5) * w, vals, width=w, color=COLORS[k], label=LABELS[k])
        ax.set_ylabel(ylab)
        if lim:
            ax.set_ylim(*lim)
        ax.text(-0.08, 1.02, 'ABC'[j], transform=ax.transAxes, fontsize=10, fontweight='bold')
    axes[0].legend(frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, 1.05), fontsize=6.5)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, fontsize=6.5)
    fig.text(0.02, 0.0, 'Figure S7', fontsize=8)
    fig.savefig(FIG + 'Figure_S7.pdf', bbox_inches='tight')
    fig.savefig(FIG + 'Figure_S7.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    print('wrote Figure_S7')


if __name__ == '__main__':
    main()
