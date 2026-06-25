import re
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd

import config

MP = config.MODELING_PARAMS


# ── heuristic theme helpers (gsea_heuristic_signatures) ────────────────────
def load_gsea(ph, gsea_dir=None):
    gsea_dir = gsea_dir or config.GSEA_DIR
    f = gsea_dir / f"gsea_result_{ph.replace('/', '_')}.csv"
    if not f.exists():
        f = gsea_dir / f"gsea_result_{ph}.csv"
    return pd.read_csv(f) if f.exists() else None


def _theme_mask(terms, subs):
    rx = re.compile('|'.join(r'(?<![A-Za-z])' + re.escape(s) for s in subs), re.IGNORECASE)
    return terms.apply(lambda t: bool(rx.search(t)))


def _strip_code(name):
    for sep in (' R-HSA-', ' (GO:', ' WP'):
        name = name.split(sep)[0]
    return name.strip()


@dataclass
class SigContext:
    dis_pheno: object
    sym_to_idx: dict
    phenos_u: list
    meanZ: dict
    samp_n: dict


def make_context(dd):
    """DiseaseData -> signature 시각화 컨텍스트 (sym_to_idx, meanZ, samp_n)."""
    sym_to_idx = {}
    for i, s in enumerate(dd.gene_syms):
        if s not in sym_to_idx:
            sym_to_idx[s] = i
    phenos_u = sorted(np.unique(dd.dis_pheno))
    meanZ = {p: dd.Z_dis[dd.dis_pheno == p].mean(axis=0) for p in phenos_u}
    samp_n = {p: int((dd.dis_pheno == p).sum()) for p in phenos_u}
    return SigContext(dd.dis_pheno, sym_to_idx, phenos_u, meanZ, samp_n)


def theme_rows(ph, ctx, themes=None, gsea_dir=None, cap_per_theme=None):
    """fig_phenotype 데이터 빌드: (rows, bands, lead_pool, yc) 반환 (시각화 분리용)."""
    themes = (themes or THEMES)[ph]
    cap = MP['sig_cap_per_theme'] if cap_per_theme is None else cap_per_theme
    g = load_gsea(ph, gsea_dir).copy()
    g['NES'] = pd.to_numeric(g['NES'], errors='coerce')
    g['clean'] = g['Term'].str.split('__').str[1].fillna(g['Term']).map(_strip_code)
    g['tag'] = g['Tag %'].astype(str).str.split('/').str[0].replace('nan', np.nan).astype(float)

    rows, bands, lead_pool, yc = [], [], [], 0
    for ti, (lab, subs, nov) in enumerate(themes):
        sub = g[_theme_mask(g['clean'], subs)].copy()
        ntot = len(sub)
        if ntot == 0:
            continue
        sub = sub.reindex(sub['NES'].abs().sort_values(ascending=False).index).head(cap)
        sub = sub.sort_values('NES')
        y0 = yc
        for _, r in sub.iterrows():
            rows.append(dict(y=yc, nes=r['NES'], tag=r['tag'] if pd.notna(r['tag']) else 5,
                             term=r['clean'], theme=ti))
            if isinstance(r['Lead_genes'], str):
                lead_pool += r['Lead_genes'].split(';')
            yc += 1
        bands.append(dict(ti=ti, y0=y0 - 0.5, y1=yc - 0.5,
                          label=f"{lab} (n={ntot})" + (' *' if nov else '')))
        yc += 0.6
    return rows, bands, lead_pool, yc


# ── EnrichmentMap leading-edge clustering (gsea_enrichmentmap_signatures) ───
def get_sig(ph, direction='up', gsea_dir=None):
    gsea_dir = gsea_dir or config.GSEA_DIR
    f = gsea_dir / f"gsea_result_{ph.replace('/', '_')}.csv"
    g = pd.read_csv(f)
    g['NES'] = pd.to_numeric(g['NES'], errors='coerce')
    g['clean'] = g['Term'].str.split('__').str[1].fillna(g['Term'])
    g = g[g['NES'].notna() & g['Lead_genes'].notna()].copy()
    if direction == 'up':
        g = g[g['NES'] > 0]
    elif direction == 'down':
        g = g[g['NES'] < 0]
    return g.reset_index(drop=True)


def _lead_sets(df):
    return [set(str(s).split(';')) for s in df['Lead_genes']]


def emap_sim(a, b, min_size=3):
    if len(a) < min_size or len(b) < min_size:
        return 0.0
    inter = len(a & b)
    if inter == 0:
        return 0.0
    return 0.5 * (inter / len(a | b)) + 0.5 * (inter / min(len(a), len(b)))


def cluster_pathways(df, sim_thr=None, resolution=None):
    """leading-edge 유사도 그래프 + Leiden. resolution=None이면 modularity 동적 탐색."""
    sim_thr = MP['emap_sim_thr'] if sim_thr is None else sim_thr
    sets = _lead_sets(df)
    n = len(df)
    edges, w = [], []
    for i, j in combinations(range(n), 2):
        s = emap_sim(sets[i], sets[j])
        if s >= sim_thr:
            edges.append((i, j))
            w.append(s)
    G = ig.Graph(n=n, edges=edges)
    if w:
        G.es['weight'] = w
    if resolution is not None:
        best_part = la.find_partition(
            G, la.RBConfigurationVertexPartition,
            weights='weight' if w else None,
            resolution_parameter=resolution, seed=42)
    else:
        best_part, best_mod = None, -1.0
        for res in [0.3, 0.5, 0.8, 1.0, 1.2]:
            part = la.find_partition(
                G, la.RBConfigurationVertexPartition,
                weights='weight' if w else None,
                resolution_parameter=res, seed=42)
            mod = G.modularity(part.membership, weights='weight' if w else None)
            if mod > best_mod:
                best_mod, best_part = mod, part
    return np.array(best_part.membership), edges, w


_STOP = {'of', 'the', 'and', 'to', 'in', 'by', 'via', 'an', 'a', 'signaling', 'pathway', 'pathways',
         'process', 'regulation', 'positive', 'negative', 'response', 'complex', 'activity',
         'cell', 'protein', 'mediated', 'dependent', 'formation', 'involved', 'its', 'for', 'on'}


def label_cluster(terms, k=3):
    toks = []
    for t in terms:
        for wd in re.split(r'[\s,\-/():]+', t.lower()):
            if (wd and wd not in _STOP and len(wd) > 2
                    and not wd.replace('.', '').isdigit() and not wd.startswith('hsa')):
                toks.append(wd)
    top = [w for w, _ in Counter(toks).most_common(k)]
    return ' / '.join(top) if top else '(misc)'


def cluster_summary(ph, direction='up', sim_thr=None, resolution=None,
                    save=True, save_network=True, cluster_dir=None):
    """군집 요약표 + (옵션) Cytoscape용 GraphML 저장."""
    cluster_dir = cluster_dir or (config.GSEA_DIR / 'Clusters')
    cluster_dir.mkdir(parents=True, exist_ok=True)
    df = get_sig(ph, direction)
    if len(df) == 0:
        return df, pd.DataFrame()
    labels, edges, w = cluster_pathways(df, sim_thr, resolution)
    df['cluster'] = labels
    rows = []
    for c, sub in df.groupby('cluster'):
        rep = sub.loc[sub['NES'].abs().idxmax()]
        rows.append(dict(cluster=int(c), size=len(sub),
                         auto_label=label_cluster(list(sub['clean'])),
                         representative=rep['clean'], rep_NES=round(rep['NES'], 2),
                         mean_NES=round(sub['NES'].mean(), 2)))
    summ = pd.DataFrame(rows).sort_values('size', ascending=False).reset_index(drop=True)
    if save:
        fn_csv = f"clusters_{ph.replace(' ', '_').replace('/', '_')}_{direction}.csv"
        summ.to_csv(cluster_dir / fn_csv, index=False)
    if save_network and len(edges) > 0:
        G_export = ig.Graph(n=len(df), edges=edges)
        if w:
            G_export.es['weight'] = w
        G_export.vs['name'] = df['clean'].tolist()
        G_export.vs['NES'] = df['NES'].astype(float).tolist()
        G_export.vs['cluster'] = [int(lbl) for lbl in labels]
        fn_graphml = f"network_{ph.replace(' ', '_').replace('/', '_')}_{direction}.graphml"
        G_export.write_graphml(str(cluster_dir / fn_graphml))
    return df, summ


# ── per-phenotype theme registry (heuristic) — verbatim ─────────────────────
THEMES = {
 'HIV':[
   ('IFN / antiviral ISG',['Interferon','Antiviral','ISG15','Defense Response To Virus','Viral Genome Replication','Defense Response To Symbiont'],False),
   ('Neutrophil / NET',['Neutrophil','Extracellular Trap'],False),
   ('Senescence / SASP',['Senescence','SASP','Telomere'],False),
   ('HIV Rev–host (pos ctrl)',['Import Of Rev','Interactions Of Rev','Rev-mediated','Rev Protein'],True),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S'],False),
   ('Mitochondrial (↓)',['Mitochondrial'],False)],
 'HIV + Tuberculosis':[
   ('SLE / NET / histone',['Lupus','Neutrophil','Extracellular Trap','Histone'],False),
   ('Pyroptosis / inflammasome',['Pyroptosis','Inflammasome'],False),
   ('DNA replication / cell cycle',['DNA Replication','Cell Cycle','Checkpoint','Pre-Replicative','ORC'],False),
   ('B cell (↓)',['B Cell'],True),
   ('Translation / GCN2 (↓)',['Translation','Ribosom','rRNA','Peptide Chain','EIF2AK4','40S','60S'],False)],
 'Tuberculosis':[
   ('Phagocytosis / Mtb',['Phagocyt','Mycobacterium','Phagosom','FCGR'],False),
   ('Interferon (α/β + γ)',['Interferon'],False),
   ('Glycolysis / metabolism',['Glycolytic','Carbohydrate Catabolic'],False),
   ('Autophagy',['Autophagy'],False),
   ('Spliceosome / splicing (↓)',['Splic','Spliceosome','mRNA Processing'],True),
   ('B cell (↓)',['B Cell'],False)],
 'Liver Cancer (Chen)':[
   ('Complement / coagulation',['Complement','Coagulation','Fibrinolysis','Platelet'],False),
   ('Lipoprotein / cholesterol',['Lipoprotein','Cholesterol','HDL','Apolipoprotein'],False),
   ('IGF / secreted protein',['IGF Transport','IGFBP','Post-translational Protein Phosphorylation'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S'],False)],
 'Liver Cancer (Roskams-Hieter)':[
   ('Complement / coagulation',['Complement','Coagulation','Fibrinolysis'],False),
   ('Cell cycle / G1-S',['G1/S','Cell Cycle','Senescence'],False),
   ('YAP / Hippo',['YAP','Hippo','TEAD','WWTR1'],False),
   ('GTPase / cytoskeleton / EMT',['GTPase','Adherens','Actin','RHO'],False),
   ('Antiviral (HBV/HCV bg)',['Viral Genome Replication','Antiviral','Interferon'],False),
   ('OXPHOS (↓)',['Oxidative phosphorylation','Respiratory Electron','Mitochondrial Respiratory'],False)],
 'Colorectal Cancer':[
   ('Hippo / YAP',['Hippo','YAP','WWTR1','TEAD'],True),
   ('Integrin / ECM',['Integrin','Adhesion'],False),
   ('RTK (ERBB2/KIT/EGFR)',['ERBB2','KIT','EGFR','STAT'],False),
   ('Insulin / glucose',['Insulin','Glucose'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S','NMD'],False)],
 'Esophagus Cancer (Chen)':[
   ('Embryonic / developmental',['Embryonic','Mammary','Morphogenesis','Differentiation'],False),
   ('YAP / Hippo',['YAP','Hippo','TEAD','WWTR1'],False),
   ('FGFR2 signaling',['FGFR'],True),
   ('IGF1R signaling',['IGF1R','IGF'],False),
   ('Neuro / ion channel',['GABA','Nicotine','Synaptic','Ion Transport','Ion Channel','Neurotransmitter','Glutamate Receptor'],False),
   ('OXPHOS + translation (↓)',['Oxidative','Respiratory Electron','Translation','Ribosom','rRNA','Mitochondrial'],False)],
 'Stomach Cancer':[
   ('Cell-cell junction / desmosome',['Junction','Adherens','Desmosom','Cadherin'],False),
   ('Integrin / ECM',['Integrin','Extracellular Matrix','Matrix'],False),
   ('RET signaling',['RET Signaling'],False),
   ('Developmental morphogenesis',['Morphogenesis','Mammary','Gonad'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S'],False)],
 'Pancreatic Cancer (Moore)':[
   ('Complement / coagulation / platelet',['Complement','Coagulation','Platelet'],False),
   ('ECM / desmoplasia',['Extracellular Matrix','Collagen','ECM','Non-integrin'],False),
   ('Rho GTPase / invasion',['RHO','GTPase'],False),
   ('Hippo / YAP / RTK',['Hippo','YAP','ERBB','EGFR','IGF'],False),
   ('TCR / T cell (↓)',['TCR','ZAP-70','CD3','Immunological Synapse','Primary Immunodeficiency'],True),
   ('OXPHOS / NF-κB (↓)',['Oxidative phosphorylation','Respiratory Electron','NF-kB','NF-kappaB'],False)],
 'Pancreatitis':[
   ('Hippo / YAP',['Hippo','YAP','WWTR1'],False),
   ('Rho GTPase / focal adhesion',['RHO','GTPase','Focal Adhesion','Cell-Substrate'],False),
   ('NOTCH',['NOTCH'],False),
   ('CFTR (↓)',['CFTR','Cystic Fibrosis'],True),
   ('OXPHOS / mito translation (↓)',['Oxidative','Respiratory Electron','Mitochondrial'],False)],
 'Lung Cancer':[
   ('Integrin / adhesion',['Integrin','Adhesion','Elastic'],False),
   ('RTK (ERBB2/EGFR/RET)',['ERBB2','EGFR','RET Signaling','GAB1'],False),
   ('Hippo / YAP',['Hippo','YAP','WWTR1'],False),
   ('Angiogenesis',['Angiogenesis','Vascular','Endothelial'],False),
   ('Insulin / glycogen',['Insulin','Glycogen'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S'],False)],
 'ICI-m':[
   ('Lysosome / proteostasis',['Lysosome','Proteasome','Vacuolar'],False),
   ('OXPHOS / mito (↑)',['Oxidative phosphorylation','Respiratory','Proton Motive','Aerobic','ATP Synthesis'],False),
   ('Antigen cross-presentation',['Cross Presentation','Cross-presentation','ER-Phagosome','Antigen Processing'],True),
   ('ER stress / ERAD',['ERAD','Endoplasmic Reticulum','Protein processing'],True),
   ('Neutrophil',['Neutrophil'],False),
   ('Translation (↓)',['Translation','Ribosom','Peptide Chain','40S','60S','EIF2AK4'],False)],
 'ICI-treated Cancer':[
   ('Lysosome / glycan degradation',['Lysosome','Vacuolar','Glycan','Glycosaminoglycan'],False),
   ('Antigen presentation / MHC-I',['Antigen Processing','MHC Class I','Presentation'],False),
   ('Alloimmune / autoimmune',['Allograft','Graft','Autoimmune','Thyroid','Type I Diabetes'],False),
   ('Mucosal IgA',['IgA','Intestinal Immune'],True),
   ('Innate / defensin',['Defensin','Antifungal','ROS And RNS'],False),
   ('Translation / SLIT-ROBO (↓)',['Translation','Ribosom','Peptide Chain','SLIT','ROBO'],False)],
 'MM':[
   ('Cell cycle / mitosis / kinesin',['Kinesin','Cell Cycle','G1/S','Polo-like','Chromosome','Spindle','Sister Chromatid','Mitotic'],False),
   ('Erythrocyte / Hb',['Gas Transport','Carbon Dioxide','Oxygen Transport','Erythrocyte'],False),
   ('Epigenetic / chromatin',['HDAC','Histone','Nucleosome','DNA Methylation','Telomere'],False),
   ('Immunoglobulin / B cell (↓)',['Immunoglobulin','B Cell','MHC Class II'],True),
   ('Translation / OXPHOS (↓)',['Translation','Ribosom','Peptide Chain','Oxidative phosphorylation'],False)],
 'MGUS':[
   ('Epigenetic / chromatin',['DNA Methylation','SIRT1','PRC2','HDAC','Histone','Chromosome Condensation','Telomere'],False),
   ('SLE / NET',['Lupus','Neutrophil','Extracellular Trap'],False),
   ('YAP / NOTCH',['YAP','WWTR1','NOTCH'],False),
   ('Senescence',['Senescence'],False),
   ('NK / antigen presentation (↓)',['Natural Killer','NK Cell','Antigen Processing','MHC Class I','Lysosome'],False)],
 'ME/CFS':[
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S','NMD'],False),
   ('FoxO signaling (↓)',['FoxO'],True),
   ('Insulin / glucose (↓)',['Insulin','Glucose'],True),
   ('Integrated stress / L13a (↓)',['EIF2AK4','L13a','Selenocysteine','Selenoamino'],False),
   ('Antiviral / coronavirus (↓)',['Coronavirus','Interferon','Antiviral'],False),
   ('Estrogen (↓)',['Estrogen'],False)],
 'CAD_HF+':[
   ('FGF / FGFR1 axis (↑)',['FGFR'],True),
   ('Cardiac fibrosis / ECM (↑)',['Collagen','Non-integrin membrane-ECM','Glycosaminoglycan','Aminoglycan'],True),
   ('Ion-channel remodeling (↑)',['Calcium Ion','Sodium Ion','Potassium Ion','Sarcoplasmic Reticulum','Ankyrin'],True),
   ('Pro-coagulant / hemostasis (↑)',['Coagulation','Fibrinolysis','Hemostasis'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S','NMD'],False),
   ('OXPHOS / mito (↓)',['Oxidative','Respiratory Electron','Complex I','TCA','Mitochondrial'],False)],
 'CAD_HF-':[
   ('Sensory / pain (↑)',['Sensory Perception','Pain'],False),
   ('Innate β-defensin (↑)',['Defensin'],False),
   ('Mucin O-glycosylation (↑)',['GALNT','C1GALT','Mucin'],False),
   ('Purinergic / GAG (↑)',['Purinergic','Nucleotide-like','Glycosaminoglycan','GAG'],False),
   ('Coagulation regulation (↑)',['Coagulation','Fibrinolysis'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S','NMD'],False),
   ('OXPHOS / mito (↓)',['Oxidative','Respiratory Electron','Complex I','TCA','Mitochondrial'],False)],
 'Pre-eclampsia':[
   ('Neurexin / synapse (↑)',['Neurexin','Neuroligin','Postsynaptic','Synaptic'],True),
   ('DNA methylation (↑)',['DNA Methylation'],False),
   ('Ribosome / translation (↓)',['Ribosom','Translation','rRNA','Peptide Chain','40S','60S'],False),
   ('Innate immune / IFN (↓)',['Interferon','Antiviral','ISG15','NLR','Inflammasome','Lectin'],False),
   ('Splicing / chromatin (↓)',['Splic','Spliceosome','Chromatin','Histone Deacetyl'],False),
   ('Lipoprotein metabolism (↓)',['Lipoprotein','Fatty Acid'],False)],
 'Other Cancer':[
   ('Rho GTPase (full set)',['RHO','GTPase','RAC1','CDC42'],False),
   ('Focal adhesion / cell-matrix',['Focal Adhesion','Cell-Substrate','Cell-Matrix','Adhesion'],False),
   ('VEGF / angiogenesis',['VEGF','Vascular Endothelial','Fibroblast Migration'],False),
   ('Translation / ribosome (↓)',['Translation','Ribosom','rRNA','Peptide Chain','40S','60S','L13a'],False),
   ('OXPHOS (↓)',['Oxidative','Respiratory Electron'],False),
   ('SLIT-ROBO / proteasome (↓)',['SLIT','ROBO','Proteasome'],False)],
}