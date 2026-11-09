import csv, re, os
import openpyxl
from openpyxl.styles import Font, Alignment

MD = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/metadata'
BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons'
GEO = f'{MD}/geo_cache'  # GEO / ENA / ENCODE records downloaded when the table was built

rows = list(csv.reader(open(f'{MD}/sample_metadata_used_datasets.tsv'), delimiter='\t'))
hdr, body = rows[0], rows[1:]
H = {h: i for i, h in enumerate(hdr)}
GRP = H['analysis_celltype_group(s)']
LAB = {'GSE60424': 'H1', 'GSE115736': 'H2', 'GSE116177': 'M1'}

def is_neut(r):
    s = r[H['sample_name_used']]
    return s.startswith('Neut') or s.startswith('Neutrophils')

M1_MATCHED = ('NveCD4T-', 'NveCd8T-', 'MonoPB-')  # M1-matched sets (Table S4); for B and NK cells M1 is used as is
CELL_LABEL = {'NveB': 'B', 'fibroblast': 'Fibroblasts'}  # sample set names as in Table S4
out = []
# flagstat values that were stale in the *_flagstat_sum_HN4.txt summaries (re-checked on the BAMs)
FIX = {(f[0], f[1]): f[2:4] for f in (l.rstrip('\n').split('\t') for l in list(open(f'{MD}/flagstat_corrections.tsv'))[1:])}
for r in body:
    ds = r[0]
    if ds == 'GSE180020' or is_neut(r):
        continue
    r = list(r)
    if (ds, r[H['sample_name_used']]) in FIX:
        r[H['flagstat_total_alignments(incl_secondary)']], r[H['flagstat_mapped']] = FIX[(ds, r[H['sample_name_used']])]
    if ds == 'E-MTAB-5919':
        r[GRP] = 'FH:Fibroblasts' if r[1] == 'Homo sapiens' else 'FM:Fibroblasts'
    else:
        cell = r[GRP].split(';')[0].split(':')[1]
        cell = CELL_LABEL.get(cell, cell)
        lab = f'{LAB[ds]}:{cell}'
        if ds == 'GSE116177' and r[H['sample_name_used']].startswith(M1_MATCHED):
            lab += f', M1-matched:{cell}'
        r[GRP] = lab
    out.append(r)

# ---------------- new datasets ----------------
ena = {}
for l in open(f'{GEO}/ena.tsv'):
    f = l.rstrip('\n').split('\t')
    ena[f[0]] = dict(gsm=f[1], title=f[3], instr=f[4], layout=f[5], reads=int(f[9]), bases=int(f[10]))
# SRR36402895: ENA shows 0 reads; NCBI SRA runinfo: 19,353,794 spots, 6,115,798,904 bases, avgLength 316 (I1+I2+R1+R2)
ena['SRR36402895'].update(reads=19353794, bases=19353794 * 300)

def flagstat(study, s):
    t = open(f'{BASE}/{study}/bam_files/{s}.flagstat.txt').read()
    tot = re.search(r'^(\d+) \+ \d+ in total', t, re.M).group(1)
    mp = re.search(r'^(\d+) \+ \d+ mapped', t, re.M).group(1)
    return tot, mp

def hisat(study, s):
    L = open(f'{BASE}/{study}/bam_files/{s}.hisat2.log').read().strip().split('\n')
    n = int(L[0].split()[0]); rate = L[-1].split('%')[0]
    unit = 'reads' if 'unpaired' in L[1] else 'read pairs'
    return n, unit, rate

TRIM = ('reads adapter/quality-trimmed with cutadapt before alignment (3\' Nextera CTGTCTCTTATACACATCT + TruSeq AGATCGGAAGAGC, '
        '5\' Smart-seq2 ISPCR AAGCAGTGGTATCAACGCAGAGT; -q 20, -m 30)')

G118 = dict(celltype={'CD4T': ('Naive_Teffs', 'naive conventional (effector-lineage, non-Treg) CD4+ T cells'),
                      'CD8T': ('Naive_CD8_T', 'naive CD8+ T cells'),
                      'Mono': ('Monocytes', 'resting blood monocytes (subset not specified in GEO)'),
                      'NK': ('Mature_NK', 'mature NK cells, resting'),
                      'NveB': ('Naive_B', 'naive B cells, resting')})
M2CELL = {
    'GSE261448': {
        'NK': ('NK cells (CD45+CD11-Nk1.1+)', 'FACS: CD45+ CD11b- NK1.1+ (GEO characteristic written as CD11-)',
               'Spleen', 'splenic NK cells, CD11b- gate (CD11b- NK cells are the less mature NK subset); no CD3 exclusion listed'),
        'Mono': ('Monocytes (CD45+CD11b+Ly6C+)', 'FACS: CD45+ CD11b+ Ly6C+',
                 'Bone Marrow', 'bone-marrow Ly6C+ monocytes (resting; classical Ly6C+ subset)')},
}

sheet = list(csv.DictReader(open(f'{MD}/sample_sheet_new_datasets.tsv'), delimiter='\t'))
keep = [r for r in sheet if r['study'] in ('GSE118165', 'GSE287527', 'GSE261448', 'GSE313371')
        or r['sample'] in ('CD4T_1_ENCODE4', 'CD4T_2_ENCODE4')]
order = ['GSE118165', 'GSE287527', 'GSE261448', 'GSE313371', 'ENCODE4']
keep.sort(key=lambda r: order.index(r['study']))

ENC = {
    'CD4T_1_ENCODE4': dict(exp='ENCSR968ZOA', lib='ENCLB414OYA', bs='ENCBS958MXN', gsm='GSM8694419', sex='female', age='57 days (adult)',
                           reads=38445313,
                           title='MouseF1 of C57BL/6 (f) x CAST Eij (m)_m.CD4 naive_spleen_F_57days RNA-seq'),
    'CD4T_2_ENCODE4': dict(exp='ENCSR835KFT', lib='ENCLB558ULM', bs='ENCBS997RQF', gsm='GSM8694049', sex='male', age='76 days (adult)',
                           reads=41200944,
                           title='MouseF1 of C57BL/6 (f) x CAST Eij (m)_m.CD4 naive_spleen_M_76days RNA-seq'),
}

for r in keep:
    st, s, cell = r['study'], r['sample'], r['cell']
    row = {h: '' for h in hdr}
    row['dataset'] = st
    row['species'] = 'Homo sapiens' if r['species'] == 'human' else 'Mus musculus'
    row['SRR_runs'] = r['runs']
    row['sample_name_used'] = s
    row['analysis_celltype_group(s)'] = f"{r['dataset']}:{CELL_LABEL.get(cell, cell)}"
    tot, mp = flagstat(st, s)
    row['flagstat_total_alignments(incl_secondary)'] = tot
    row['flagstat_mapped'] = mp
    n, unit, rate = hisat(st, s)
    notes = []
    if st != 'ENCODE4':
        e = ena[r['runs']]
        row['GSM'] = e['gsm']; row['GEO_title'] = e['title']; row['instrument'] = e['instr']
        row['SRA_layout'] = e['layout']
        rl = e['bases'] / e['reads'] / 2
        row['SRA_total_spots(reads_or_pairs)'] = str(e['reads'])
        row['processed_layout_in_our_pipeline'] = 'paired-end (cutadapt-trimmed)'
    if st == 'GSE118165':
        lab, state = G118['celltype'][cell]
        donor = e['title'].split('-')[0]
        row['original_celltype_label'] = f'{lab} (lineage: ' + {'CD4T': 'CD4', 'CD8T': 'CD8', 'Mono': 'MYELOID', 'NK': 'NK', 'NveB': 'B'}[cell] + ')'
        row['sorting_definition_markers'] = 'FACS-sorted; gating markers not given in GEO (see Calderon et al. 2019 Nat Genet, PMID 31570894)'
        row['tissue_source'] = 'peripheral blood (healthy adult donors; GEO source: primary blood)'
        row['cell_state'] = state
        row['treatment_condition'] = 'unstimulated (-U; GEO treatment: no_treament); stimulated (-S) samples of the same donors not used'
        row['donor_or_mouse_id'] = f'donor {donor}'
        row['sex'] = 'not given in GEO'; row['age'] = 'adult (not given in GEO)'; row['strain_or_ethnicity'] = 'not given in GEO'
        row['library_type'] = 'polyA RNA; Smart-seq2 with custom Nextera dual-index PCR primers'
        row['read_length_nt'] = '2x75' if r['runs'] == 'SRR7647808' else f'2x~{rl:.0f} (variable-length reads in SRA; mean)'
        notes.append('unstimulated (-U) samples only; Smart-seq2 libraries')
        if cell == 'NveB' and donor == '1010':
            notes.append('GEO has no unstimulated naive-B sample for donor 1003, so donor 1010 is the 4th naive-B replicate (HiSeq 4000 run, 2x75)')
    elif st == 'GSE261448':
        lab, mk, tis, state = M2CELL[st][cell]
        mouse = re.match(r'(\d+)', e['title']).group(1)
        row['original_celltype_label'] = lab; row['sorting_definition_markers'] = mk
        row['tissue_source'] = tis; row['cell_state'] = state
        row['treatment_condition'] = 'PBS control (100 ul PBS s.c., sacrificed day 3); BCG-vaccinated mice not used'
        row['donor_or_mouse_id'] = f'control mouse {mouse} (NK_{mouse} and Mono_{mouse} are from the same mouse)'
        row['sex'] = 'not given in GEO'; row['age'] = 'not given in GEO'; row['strain_or_ethnicity'] = 'not given in GEO; genotype: wt'
        row['library_type'] = 'rRNA-depleted RNA (NEBNext rRNA Depletion v2; GEO molecule field says polyA RNA); NEBNext Ultra II Directional, UDI-UMI adaptors'
        row['read_length_nt'] = f'2x{rl:.0f}'
        notes.append('PBS-control mice of a BCG vaccination study (PMID 39252057); NK from spleen, Ly6C+ monocytes from bone marrow')
    elif st == 'GSE287527':
        rep = e['title'].split('rep')[1]
        row['original_celltype_label'] = 'Follicular B cell'
        row['sorting_definition_markers'] = 'FACS: B220+ CD93- CD138- GL7- FAS- CD21+ CD23+'
        row['tissue_source'] = 'Spleen'
        row['cell_state'] = 'mature naive follicular B cells (non-GC, non-plasma, non-MZ)'
        row['treatment_condition'] = 'none; wild-type controls of a Mek1/Mek2 deletion study (Mek1-flox/Mek2+/- Vav1-icre mutants not used)'
        row['donor_or_mouse_id'] = f'wt replicate {rep} (individual mice not specified)'
        row['sex'] = 'not given in GEO'; row['age'] = '3-5 months'; row['strain_or_ethnicity'] = 'not given in GEO; genotype: wt'
        row['library_type'] = 'polyA RNA; Illumina TruSeq RNA'
        row['read_length_nt'] = f'2x{rl:.0f}'
        notes.append('wild-type spleen follicular B cells of a Mek1/Mek2 study')
    elif st == 'GSE313371':
        m = e['title'].split(',')[-1].strip()
        row['original_celltype_label'] = 'CD8 T cells'
        row['sorting_definition_markers'] = 'splenic CD8 T cells; isolation method/markers not given in GEO'
        row['tissue_source'] = 'Spleen'
        row['cell_state'] = 'total CD8 T cells (naive+memory mixed; no naive gate)'
        row['treatment_condition'] = 'uninfected (no Mtb infection); WT genotype (Nr4a1-KO and infected mice not used)'
        row['donor_or_mouse_id'] = f'WT uninfected mouse {m}'
        row['sex'] = 'not given in GEO'; row['age'] = 'not given in GEO'; row['strain_or_ethnicity'] = 'B6129 (B6;129 background), WT'
        row['library_type'] = 'total RNA (2 ng); modified Smart-seq2 + Nextera XT'
        row['read_length_nt'] = '2x150'
        notes.append('uninfected wild-type B6129 mice, spleen total CD8 T cells, from an Nr4a1/Mtb study')
        if r['runs'] == 'SRR36402895':
            notes.append('SRA run stores 4 reads per spot (I1, I2, R1, R2; avgLength 316); biological pair = SRA reads 3-4; ENA lists 0 reads, spot count from NCBI SRA runinfo')
    elif st == 'ENCODE4':
        d = ENC[s]
        row['GSM'] = f"{d['exp']} / {d['lib']} / {d['bs']} (GEO {d['gsm']})"
        row['GEO_title'] = d['title']
        row['original_celltype_label'] = 'CD4-positive naive resting alpha-beta T cell'
        row['sorting_definition_markers'] = 'not given in ENCODE metadata (biosample source: Mark Groudine lab, Fred Hutch)'
        row['tissue_source'] = 'Spleen'
        row['cell_state'] = 'naive resting CD4+ alpha-beta T cells'
        row['treatment_condition'] = 'none (no treatments in ENCODE biosample)'
        row['donor_or_mouse_id'] = f"biosample {d['bs']} (ENCODE donor strain record ENCDO509HIY)"
        row['sex'] = d['sex']; row['age'] = d['age']
        row['strain_or_ethnicity'] = 'B6CASTF1/J (F1 of C57BL/6J female x CAST/EiJ male)'
        row['instrument'] = 'Illumina HiSeq 2500'
        row['library_type'] = 'total RNA-seq; NEBNext Ultra II Directional (reverse-stranded)'
        row['SRA_layout'] = 'SINGLE (ENCODE FASTQ)'
        row['read_length_nt'] = '101'
        row['processed_layout_in_our_pipeline'] = 'single-end (cutadapt-trimmed)'
        row['SRA_total_spots(reads_or_pairs)'] = f"{d['reads']} (ENCODE file read_count)"
        notes.append('ENCODE4 naive CD4 T cells from spleen of B6CASTF1/J F1 hybrid mice (Wold lab, Caltech); F1 hybrid (B6 x CAST) reads aligned to the standard mouse reference genome')
    notes.append(TRIM)
    notes.append(f'HISAT2: {n} {unit} after trimming, overall alignment {rate}%')
    notes.append('flagstat counts include secondary alignments')
    row['notes'] = ' | '.join(notes)
    out.append([row[h] for h in hdr])

hdr = ['group_name' if h == 'analysis_celltype_group(s)' else h for h in hdr]  # group names of Table S4
with open(f'{MD}/sample_metadata_revision.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(hdr); w.writerows(out)

wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'samples'
ws.append(['Supplementary Table S2: Sample-level metadata of all samples used'])
ws['A1'].font = Font(bold=True, size=12)
ws.append(hdr)
for c in ws[2]:
    c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical='top')
NUM = {'SRA_total_spots(reads_or_pairs)', 'flagstat_total_alignments(incl_secondary)', 'flagstat_mapped'}
for r in out:
    vals = [int(v) if hdr[i] in NUM and v.isdigit() else v for i, v in enumerate(r)]
    ws.append(vals)
ws.freeze_panes = 'A3'
for i, h in enumerate(hdr, 1):
    L = max([len(h)] + [len(str(r[i - 1])) for r in out])
    ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = max(10, min(L + 2, 60))
ws.auto_filter.ref = f'A2:{openpyxl.utils.get_column_letter(len(hdr))}{len(out) + 2}'
wb.save(f'{MD}/sample_metadata_revision.xlsx')

from collections import Counter
print(Counter(r[0] + ' ' + r[1] for r in out))
print(len(out))
