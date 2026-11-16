"""
Supplementary Table S2 (full): one row per sample of the 3 human x 3 mouse dataset design (+ fibroblasts).
Samples = the 'sample' column of events/id50_6ds/sample_design.tsv (checked at the end).
H1 GSE60424, H2 GSE115736, H3 GSE118165; M1 GSE287527 (B) / GSE215775 (CD4 T) / GSE313371 (CD8 T) / GSE261448 (NK,
monocytes); M2 GSE307667 / GSE184496 (young) / GSE263132 / GSE236819 / GSE254202; M3 GSE231772 / GSE134322 /
GSE284441 / GSE180053 / PRJNA797236 (blood); FH/FM E-MTAB-5919. GSE116177 and ENCODE4 naive CD4 T (two replicates per
cell type) are no longer used in the main analysis and are not in the table.
Metadata of the new studies: GEO SOFT family files, ENA run reports and NCBI BioSample records cached in metadata/geo_cache/.
Output: revision_6ds/supplementary/Table_S2.xlsx (moved to full_tables/ by make_tables_b_116.py) and
revision_6ds/sample_metadata_separation.tsv
"""
import csv, re, os
from collections import defaultdict
import openpyxl
from openpyxl.styles import Font, Alignment

import build_supp_tables_116 as bst

MD = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/metadata'
BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons'
SEP = f'{BASE}/ensembl116/revision_6ds'  # outputs of this copy (the shared metadata/ files are not changed)
GEO = f'{MD}/geo_cache'  # GEO / ENA / ENCODE records downloaded when the table was built

rows = list(csv.reader(open(f'{MD}/sample_metadata_used_datasets.tsv'), delimiter='\t'))
hdr, body = rows[0], rows[1:]
H = {h: i for i, h in enumerate(hdr)}
GRP = H['analysis_celltype_group(s)']
LAB = {'GSE60424': 'H1', 'GSE115736': 'H2'}
# dataset label (H1..M3, FH, FM) and cell type of every sample used: '<sample>_<study>'
DESIGN = {l.split('\t')[0]: l.split('\t')[1] for l in list(open(f'{BASE}/ensembl116/events/id50_6ds/sample_design.tsv'))[1:]}

def is_neut(r):
    s = r[H['sample_name_used']]
    return s.startswith('Neut') or s.startswith('Neutrophils')

CELL_LABEL = {'NveB': 'B', 'fibroblast': 'Fibroblasts'}  # sample set names as in Table S5
out = []
# flagstat values that were stale in the *_flagstat_sum_HN4.txt summaries (re-checked on the BAMs)
FIX = {(f[0], f[1]): f[2:4] for f in (l.rstrip('\n').split('\t') for l in list(open(f'{MD}/flagstat_corrections.tsv'))[1:])}
for r in body:
    ds = r[0]
    if ds in ('GSE180020', 'GSE116177') or is_neut(r):  # GSE116177: two replicates per cell type, no longer used
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
        r[GRP] = lab
    out.append(r)

# ---------------- new datasets ----------------
ena = {}
for l in open(f'{GEO}/ena.tsv'):
    f = l.rstrip('\n').split('\t')
    ena[f[0]] = dict(gsm=f[1], title=f[3], instr=f[4], layout=f[5], reads=int(f[9]), bases=int(f[10]))
# SRR36402895: ENA shows 0 reads; NCBI SRA runinfo: 19,353,794 spots, 6,115,798,904 bases, avgLength 316 (I1+I2+R1+R2)
ena['SRR36402895'].update(reads=19353794, bases=19353794 * 300)

M3DIR = f'{BASE}/ensembl116/m3'  # BAMs of the studies added for the 3 x 3 design (except GSE184496, PRJNA797236)

def bamdir(study):
    return f'{M3DIR}/{study}/bam_files' if os.path.isdir(f'{M3DIR}/{study}') else f'{BASE}/{study}/bam_files'

def flagstat(study, s):
    t = open(f'{bamdir(study)}/{s}.flagstat.txt').read()
    tot = re.search(r'^(\d+) \+ \d+ in total', t, re.M).group(1)
    mp = re.search(r'^(\d+) \+ \d+ mapped', t, re.M).group(1)
    return tot, mp

def hisat(study, s):
    L = open(f'{bamdir(study)}/{s}.hisat2.log').read().strip().split('\n')
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
keep = [r for r in sheet if r['study'] in ('GSE118165', 'GSE287527', 'GSE261448', 'GSE313371')]
order = ['GSE118165', 'GSE287527', 'GSE261448', 'GSE313371']
keep.sort(key=lambda r: order.index(r['study']))

for r in keep:
    st, s, cell = r['study'], r['sample'], r['cell']
    row = {h: '' for h in hdr}
    row['dataset'] = st
    row['species'] = 'Homo sapiens' if r['species'] == 'human' else 'Mus musculus'
    row['SRR_runs'] = r['runs']
    row['sample_name_used'] = s
    row['analysis_celltype_group(s)'] = f"{DESIGN[f'{s}_{st}']}:{CELL_LABEL.get(cell, cell)}"
    tot, mp = flagstat(st, s)
    row['flagstat_total_alignments(incl_secondary)'] = tot
    row['flagstat_mapped'] = mp
    n, unit, rate = hisat(st, s)
    notes = []
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
    notes.append(TRIM)
    notes.append(f'HISAT2: {n} {unit} after trimming, overall alignment {rate}%')
    notes.append('flagstat counts include secondary alignments')
    row['notes'] = ' | '.join(notes)
    out.append([row[h] for h in hdr])

# ---------------- mouse studies added for the 3 human x 3 mouse design ----------------
# ENA run reports (ena_6ds.tsv), GEO SOFT family files (<GSE>.soft) and, for PRJNA797236 (no GEO series), the NCBI
# BioSample record (PRJNA797236_biosample.txt) and the ENA experiment record (PRJNA797236_SRX13805443_ena.xml).
ena6 = {r['run_accession']: r for r in csv.DictReader(open(f'{GEO}/ena_6ds.tsv'), delimiter='\t')}


def soft_samples(study):
    """GSM -> {'title': ..., 'ch': {characteristic: value}} from the GEO SOFT family file."""
    S, cur = {}, None
    for l in open(f'{GEO}/{study}.soft', encoding='utf-8'):
        l = l.rstrip('\n')
        if l.startswith('^SAMPLE'):
            cur = l.split('=')[1].strip(); S[cur] = {'title': '', 'ch': {}}
        elif l.startswith('^'):
            cur = None
        elif cur and l.startswith('!Sample_title = '):
            S[cur]['title'] = l.split(' = ', 1)[1]
        elif cur and l.startswith('!Sample_characteristics_ch1 = '):
            k, _, v = l.split(' = ', 1)[1].partition(': ')
            S[cur]['ch'][k.strip().lower()] = v.strip()
    return S


NG = 'not given in GEO'
NEW = {  # study-specific curation from the GEO series/sample records (summary, overall design, protocols)
    'GSE215775': dict(
        label='naive CD4+ T cells', markers='freshly prepared naive CD4+ T (TN) cells; isolation method and markers not given in GEO',
        tissue='Spleen', state='naive CD4+ T cells, resting (freshly isolated)',
        treatment='no treatment (GEO); wild-type (GEO genotype: WT); mesenteric lymph node TN samples of the same study not used',
        strain='C57BL/6J (GEO: C57BL6/J), wild type', library='polyA RNA; NEBNext Ultra II Directional RNA',
        donor=lambda t, ch: f"spleen TN replicate {t.split('rep')[-1]} (individual mice not specified)",
        note='wild-type spleen naive CD4 T cells of a study comparing spleen and mesenteric lymph node TN cells (spleen only used)'),
    'GSE307667': dict(
        label='Follicular B cells', markers='FACS-sorted follicular B cells; gating markers not given in GEO',
        tissue='Spleen', state='mature naive follicular B cells, resting',
        treatment='none; Mb1-Cre (Cre-only) controls of an Arid2 conditional knockout study (GEO genotype: Mb1-Cre; the only '
                  'control group of the study); Mb1-Cre Arid2-F/F knockouts and the bone-marrow pro-B and pre-B samples not used',
        strain=f'{NG}; genotype: Mb1-Cre (Cre knock-in into Cd79a, no floxed allele)',
        library='total RNA (10 ng); SMARTer Ultra Low RNA (full-length cDNA), Covaris fragmentation, Illumina adapters',
        donor=lambda t, ch: f"Mb1-Cre control mouse, biological replicate {t.split('rep')[-1].strip()}",
        note='Cre-only (Mb1-Cre) control follicular B cells of an Arid2 deletion study'),
    'GSE184496': dict(
        label='CD4+ Naive', markers='EasySep CD4+ T cell negative selection, then FACS: CD3+ CD4+ CD25- CD62L+ CD44-',
        tissue='pooled lymph nodes (axillary, brachial, cervical, inguinal) and spleen',
        state='naive CD4+ T cells, resting',
        treatment='none; young mice only (7-11 weeks); old mice (83-86 weeks) and effector memory CD4 T samples of the same '
                  'study not used here (the old vs young naive CD4 T comparison is a control analysis, Supplementary Note 1)',
        strain='C57BL/6', library='polyA mRNA (NEBNext Poly(A) mRNA Magnetic Isolation); NEBNext Ultra II DNA Library Prep',
        donor=lambda t, ch: f"young {ch['sex']} mouse {t.split('#')[-1]}",
        note='young naive CD4 T cells of a study of old vs young and male vs female mice'),
    'GSE263132': dict(
        label='Naive CD8 T cells', markers='magnetic CD8 T cell enrichment (BioLegend), then FACS: CD3+ CD8+ B220- CD44low live single cells',
        tissue='Spleen', state='naive (CD44low) CD8+ T cells, resting',
        treatment='none (GEO treatment: none); wild-type C57BL/6J mice; hemizygous MD4 transgenic mice of the same study not used',
        strain='C57BL/6J, wild type', library='polyA mRNA; Tecan Universal Plus mRNA-Seq (200 ng total RNA)',
        donor=lambda t, ch: f"wild-type mouse {t.split('_')[-1]}",
        note='wild-type C57BL/6J naive CD8 T cells (controls of MD4 B-cell-receptor transgenic mice)'),
    'GSE236819': dict(
        label='NK cells', markers='NK cell magnetic negative selection (MagniSort), then FACS sort of 50,000 NK cells per replicate; '
                                  'gating markers not given in GEO',
        tissue='Spleen', state='splenic NK cells, steady state',
        treatment='steady state, uninfected (GEO treatment: Steady-state); GEO genotype WT = littermate controls of the NK '
                  'cell-specific Xbp1 knockout (XBP1-KO) and Ire1/Xbp1 double knockout (DKO) mice (control genotype, e.g. '
                  'Cre-only or flox-only, not given in GEO); knockout samples not used',
        strain=f'{NG}; genotype: WT (littermate controls)',
        library='mRNA (GEO library type); standard Illumina protocol (VIB Nucleomics Core)',
        donor=lambda t, ch: f"WT littermate control mouse, replicate {t.split('Rep')[-1]}",
        note='WT littermate-control splenic NK cells of an Ire1/Xbp1 knockout study'),
    'GSE254202': dict(
        label='Ly6Chi monocytes', markers='sort-purified Ly6Chi monocytes (40,000-50,000 cells per sample); gating markers not given in GEO',
        tissue='Spleen', state='classical (Ly6Chi) monocytes, resting',
        treatment='untreated (GEO treatment: sort-purified (untreated)); wild type (GEO genotype: WT); Cebpa +37 enhancer '
                  'knockout mice and Ly6Clo monocytes of the same study not used',
        strain='C57BL/6, wild type',
        library='total RNA (10 ng); SMARTer Ultra Low RNA (full-length cDNA), Covaris fragmentation, Illumina adapters',
        donor=lambda t, ch: f"wild-type mouse {t.split('WT')[-1]}",
        note='wild-type splenic Ly6Chi monocytes of a C/EBPa enhancer knockout study'),
    'GSE231772': dict(
        label='FoB (follicular B cells)', markers=NG,
        tissue=NG, state='follicular B cells of naive (uninfected, day 0) mice',
        treatment='naive, day 0 (no Brucella infection); BCL6fl/fl (flox-only, no Cre) controls of FoxP3-Cre BCL6fl/fl '
                  '(T follicular regulatory cell-deficient) mice, the control group of the study; FoxP3-Cre-BCL6fl/fl mice '
                  'and infected mice not used',
        strain=f'{NG}; genotype: BCL6fl/fl', library='polyA mRNA; stranded mRNA library (University of Missouri Genomics Core)',
        donor=lambda t, ch: f"BCL6fl/fl naive mouse {t.split('[')[-1].rstrip(']')}",
        note='flox-only (BCL6fl/fl) follicular B cells of uninfected mice'),
    'GSE134322': dict(
        label='Naive CD4 T cells', markers='negatively isolated ex vivo (STEMCELL Technologies kit; not FACS); markers not given in GEO',
        tissue='Spleen', state='naive CD4+ T cells, resting',
        treatment='no special treatment (GEO); wild-type littermates of PIK3CD gain-of-function (heterozygous) mice; '
                  'PIK3CD GOF mice not used',
        strain='C57BL/6, wild type', library='polyA RNA; Illumina TruSeq RNA',
        donor=lambda t, ch: f"wild-type littermate mouse {t.split(' ')[-1]}",
        note='wild-type littermate naive CD4 T cells of a PIK3CD gain-of-function study; sex from the GEO series design '
             '(8 week, male); GEO series summary names HiSeq 2500, the sample and SRA records HiSeq X Ten'),
    'GSE284441': dict(
        label='Naive CD8+ T cells', markers='naive CD8+ T cell isolation kit (magnetic; not FACS); markers not given in GEO',
        tissue='Spleen', state='naive CD8+ T cells, resting',
        treatment='untreated cells of uninfected, specific-pathogen-free wild-type mice (GEO genotype: WT); virtual memory '
                  '(TVM) samples and Aiolos-deficient (Ikzf3-/-) mice not used',
        strain=f'{NG}; genotype: WT', library='total RNA; standard library preparation by Azenta Life Sciences',
        donor=lambda t, ch: f"wild-type mouse {t.split(' ')[-1]}",
        note='wild-type naive CD8 T cells of an Aiolos virtual-memory T cell study'),
    'GSE180053': dict(
        label='NK', markers='FACS: CD45.2+ CD3- NK1.1+ NKp46+',
        tissue='Spleen', state='splenic NK cells, resting',
        treatment='none; Mettl3fl/fl (flox-only, no Cre; GEO: WT) controls of Ncr1-Cre Mettl3fl/fl knockout mice, the '
                  'control group of the study; knockout mice not used',
        strain='C57BL/6; genotype: Mettl3fl/fl', library='polyA RNA; NEBNext Ultra RNA Library Prep',
        donor=lambda t, ch: f"Mettl3fl/fl mouse {t.split('-')[-1]}",
        note='flox-only (Mettl3fl/fl) control splenic NK cells of a METTL3 knockout study'),
    'PRJNA797236': dict(
        label='Ly6Chi Monocyte (BioSample cell type)',
        markers='FACS-sorted on an Aria III (3,000-20,000 cells into TRIzol); gating markers not given in the BioSample or SRA records',
        tissue='Peripheral blood (BioSample tissue); spleen Ly6Chi monocytes of the same study not used here (blood vs spleen '
               'is a control analysis, Supplementary Note 1)',
        state='classical (Ly6Chi) monocytes, resting',
        treatment='none (no treatment in the BioSample or SRA records)',
        strain='not given (BioSample strain: missing)',
        library='SMART-Seq HT (full-length cDNA, 12 cycles) + TruePrep library kit (Vazyme)',
        donor=lambda t, ch: f"replicate {int(t.split('_')[-1])} (one BioSample, SAMN24975549, for the four blood replicates; "
                            f"individual mice not specified)",
        note='blood Ly6Chi monocytes of a dendritic cell type 3 ontogeny study (no GEO series; metadata from NCBI BioSample '
             'SAMN24975549 and SRA experiment SRX13805443)'),
}
# PRJNA797236 (BioSample SAMN24975549): sex male, age 6-8 weeks, tissue Peripheral blood, cell type Ly6Chi Monocyte
PRJ_BS = {'sex': 'male', 'age': '6-8 weeks'}

new_sheet = [r for r in csv.DictReader(open(f'{M3DIR}/sample_sheet_m3.tsv'), delimiter='\t')]
new_sheet += [r for r in sheet if r['study'] in ('GSE184496', 'PRJNA797236')]
new_sheet = [r for r in new_sheet if f"{r['sample']}_{r['study']}" in DESIGN]
soft = {st: soft_samples(st) for st in NEW if st.startswith('GSE')}
for r in new_sheet:
    st, s, cell = r['study'], r['sample'], r['cell']
    c = NEW[st]
    runs = r['runs'].split(',')
    E = [ena6[x] for x in runs]
    gsm = E[0]['sample_alias'] if st.startswith('GSE') else E[0]['sample_accession']
    assert all((e['sample_alias'] if st.startswith('GSE') else e['sample_accession']) == gsm for e in E)
    reads = sum(int(e['read_count']) for e in E); bases = sum(int(e['base_count']) for e in E)
    row = {h: '' for h in hdr}
    row['dataset'] = st; row['species'] = 'Mus musculus'; row['SRR_runs'] = r['runs']; row['sample_name_used'] = s
    row['analysis_celltype_group(s)'] = f"{DESIGN[f'{s}_{st}']}:{CELL_LABEL.get(cell, cell)}"
    if st.startswith('GSE'):
        g = soft[st][gsm]; title, ch = g['title'], g['ch']
        row['GSM'] = gsm
    else:
        title, ch = E[0]['library_name'], PRJ_BS
        row['GSM'] = f'{gsm} (BioSample; no GEO record)'
    row['GEO_title'] = title if st.startswith('GSE') else f"{E[0]['sample_title']} (library {title})"
    row['original_celltype_label'] = c['label']; row['sorting_definition_markers'] = c['markers']
    row['tissue_source'] = c['tissue']; row['cell_state'] = c['state']; row['treatment_condition'] = c['treatment']
    row['donor_or_mouse_id'] = c['donor'](title, ch)
    sex = ch.get('sex', '')
    row['sex'] = {'F': 'female', 'M': 'male'}.get(sex, sex) or ('male (GEO series design)' if st == 'GSE134322' else NG)
    row['age'] = (ch.get('age', '').replace('-week old', ' weeks').replace('73 days/10.4 weeks old', '73 days (10.4 weeks)')
                  or NG)
    row['strain_or_ethnicity'] = c['strain']
    row['instrument'] = E[0]['instrument_model']; row['library_type'] = c['library']
    row['SRA_layout'] = E[0]['library_layout']
    row['SRA_total_spots(reads_or_pairs)'] = str(reads)
    notes = [c['note']]
    if st == 'GSE236819':
        row['read_length_nt'] = '76 (read 1 only)'
        row['processed_layout_in_our_pipeline'] = 'single-end, read 1 only (cutadapt-trimmed)'
        notes.append('SRA spots hold read 1 (76 nt) and an index read (ENA layout: PAIRED); only read 1 was aligned, single-end')
    elif E[0]['library_layout'] == 'SINGLE':
        row['read_length_nt'] = f'{bases / reads:.0f}'
        row['processed_layout_in_our_pipeline'] = 'single-end (cutadapt-trimmed)'
    else:
        row['read_length_nt'] = f'2x{bases / reads / 2:.0f}'
        row['processed_layout_in_our_pipeline'] = 'paired-end (cutadapt-trimmed)'
    if len(runs) > 1:
        notes.append(f'{gsm} sequenced in {len(runs)} SRA runs, merged')
    tot, mp = flagstat(st, s)
    row['flagstat_total_alignments(incl_secondary)'] = tot; row['flagstat_mapped'] = mp
    n, unit, rate = hisat(st, s)
    notes.append(TRIM)
    notes.append(f'HISAT2: {n} {unit} after trimming, overall alignment {rate}%')
    notes.append('flagstat counts include secondary alignments')
    row['notes'] = ' | '.join(notes)
    out.append([row[h] for h in hdr])

# order: dataset label, then cell type, as in the design; the table must hold exactly the samples of the design
LORD = ['H1', 'H2', 'H3', 'M1', 'M2', 'M3', 'FH', 'FM']
CORD = ['B', 'CD4T', 'CD8T', 'NK', 'Mono', 'Fibroblasts']
out.sort(key=lambda r: (LORD.index(r[GRP].split(':')[0]), CORD.index(r[GRP].split(':')[1])))


def design_key(r):
    ds, s = r[0], r[H['sample_name_used']]
    if ds == 'E-MTAB-5919':
        return f"{s}_EMTAB5919{'H' if r[1] == 'Homo sapiens' else 'M'}"
    return f'{s}_{ds}'


keys = [design_key(r) for r in out]
assert len(keys) == len(set(keys)) and set(keys) == set(DESIGN), (set(DESIGN) ^ set(keys))
assert all(DESIGN[k] == r[GRP].split(':')[0] for k, r in zip(keys, out))

hdr = ['group_name' if h == 'analysis_celltype_group(s)' else h for h in hdr]  # group names of Table S5
os.makedirs(f'{SEP}/supplementary', exist_ok=True)
with open(f'{SEP}/sample_metadata_separation.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(hdr); w.writerows(out)

wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'samples'
ws.append(['Supplementary Table S2: Sample-level metadata of all samples used'])
ws['A1'].font = Font(name='Times New Roman', bold=True, size=12)
ws.append(hdr)
for c in ws[2]:
    c.font = Font(name='Times New Roman', bold=True); c.alignment = Alignment(wrap_text=True, vertical='top')
NUM = {'SRA_total_spots(reads_or_pairs)', 'flagstat_total_alignments(incl_secondary)', 'flagstat_mapped'}
for r in out:
    vals = [int(v) if hdr[i] in NUM and v.isdigit() else v for i, v in enumerate(r)]
    ws.append(vals)
ws.freeze_panes = 'A3'
for i, h in enumerate(hdr, 1):
    L = max([len(h)] + [len(str(r[i - 1])) for r in out])
    ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = max(10, min(L + 2, 60))
ws.auto_filter.ref = f'A2:{openpyxl.utils.get_column_letter(len(hdr))}{len(out) + 2}'
wb.save(f'{SEP}/supplementary/Table_S2.xlsx')
bst.times_new_roman(f'{SEP}/supplementary/Table_S2.xlsx')

from collections import Counter
print(Counter(r[0] + ' ' + r[GRP] for r in out))
print(len(out))
