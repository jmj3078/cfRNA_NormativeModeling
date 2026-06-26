import json
import re
from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

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
    """Build signature visualization context (sym_to_idx, meanZ, samp_n) from a DiseaseData."""
    sym_to_idx = {}
    for i, s in enumerate(dd.gene_syms):
        if s not in sym_to_idx:
            sym_to_idx[s] = i
    phenos_u = sorted(np.unique(dd.dis_pheno))
    meanZ = {p: dd.Z_dis[dd.dis_pheno == p].mean(axis=0) for p in phenos_u}
    samp_n = {p: int((dd.dis_pheno == p).sum()) for p in phenos_u}
    return SigContext(dd.dis_pheno, sym_to_idx, phenos_u, meanZ, samp_n)


def theme_rows(ph, ctx, themes=None, gsea_dir=None, cap_per_theme=None):
    """Build per-theme lollipop rows for a phenotype; return (rows, bands, lead_pool, yc)."""
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
    """Build leading-edge similarity graph and Leiden clusters; return (membership, edges, weights)."""
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
                    save=True, save_network=True, gsea_dir=None, cluster_dir=None):
    """Compute Leiden cluster summary; optionally save CSV and GraphML network."""
    gsea_dir = gsea_dir or config.GSEA_DIR
    cluster_dir = cluster_dir or (gsea_dir / 'Clusters')
    cluster_dir.mkdir(parents=True, exist_ok=True)
    df = get_sig(ph, direction, gsea_dir=gsea_dir)
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


with open(Path(__file__).parent / 'themes.json', encoding='utf-8') as _f:
    THEMES = json.load(_f)

