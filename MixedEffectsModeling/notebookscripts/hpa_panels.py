"""Cell-type marker panels from HPA 25.1 (proteinatlas.tsv).

panel_table: official single-cell type specificity ('Cell type enriched' / 'Group enriched').
tct_panel_table: cell type taken from the HPA tissue cell type section (enrichment within the
tissue), kept only when the gene is also body-wide specific to that organ (tissue specificity) or
to that cell type (single-cell specificity), enhanced included. Single-cell enriched genes of the
panel stay in even when the tissue cell type section misses them (e.g. CTRB2), so the panel is a
superset of panel_table before the cap. Genes are ranked enriched before
enhanced, then by HPA specificity score, and capped per panel. Panels HPA has no tissue cell type
for fall back to the single-cell enriched rule under the same cap. Shared by
6_disease_tissue_signal and 7_pancreatic_cancer_cohort_deepdive.
"""
import pandas as pd

import MixedEffectsModeling.config as config

ENRICHED = ('Cell type enriched', 'Group enriched')
TISSUE_ENRICHED = ('Tissue enriched', 'Group enriched')

# panel: (HPA tissue cell type labels, organ in 'RNA tissue specific nTPM'); a label starting
# with '- ' matches that cell type in every tissue
TCT = {
    'hepatocyte': (['Liver - Hepatocytes'], 'liver'),
    'cholangiocyte': (['Liver - Cholangiocytes'], 'liver'),
    'hepatic stellate': (['Liver - Hepatic stellate cells'], 'liver'),
    'kupffer': (['Liver - Kupffer cells'], 'liver'),
    'pancreatic acinar': (['Pancreas - Exocrine glandular cells'], 'pancreas'),
    'pancreatic duct': (['Pancreas - Ductal cells'], 'pancreas'),
    'pancreatic islet': (['Pancreas - Alpha cells', 'Pancreas - Beta cells'], 'pancreas'),
    'colonocyte': (['Colon - Colon enterocytes'], 'intestine'),
    'gastric epithelium': (['Stomach - Chief cells', 'Stomach - Gastric mucous cells', 'Stomach - Parietal cells'],
                           'stomach 1'),
    'alveolar type 2': (['Lung - Alveolar cells type 2'], 'lung'),
    'alveolar type 1': (['Lung - Alveolar cells type 1'], 'lung'),
    'airway': (['Lung - Respiratory ciliated cells'], 'lung'),
    'cardiomyocyte': (['Heart muscle - Cardiomyocytes'], 'heart muscle'),
    'proximal tubule': (['Kidney - Proximal tubular cells'], 'kidney'),
    'distal tubule / loop of Henle': (['Kidney - Ascending Loop of Henle cells'], 'kidney'),
    'skeletal muscle': (['Skeletal muscle - Skeletal myocytes'], 'skeletal muscle'),
    'macrophage': (['- Macrophages', 'Liver - Kupffer cells'], None),
    'vascular endothelium': (['- Endothelial cells', 'Liver - Vascular endothelial cells',
                              'Liver - Sinusoid Endothelial cells'], None),
    'fibroblast': (['- Fibroblasts', '- Fibroblasts_1', '- Fibroblasts_2'], None),
    'smooth muscle': (['- Smooth muscle cells'], None),
    'plasma cell': (['- Plasma cells'], None),
    'T cell': (['- T-cells'], None),
    'B / NK cell': (['Spleen - B-cells (Spleen)', 'Lung - B-cells (Lung)', '- NK-cells (Liver)', '- NK-cells (Lung)',
                     '- NK-cells (Spleen)'], None),
    'neutrophil': (['- Neutrophils', '- Neutrophils_1', '- Neutrophils_2'], None),
    'erythroid': (['Liver - Erythroid cells', 'Spleen - Erythrocytes'], None),
    'mesothelium': (['Adipose visceral - Mesothelial cells'], None),
}


def load_specificity():
    a = pd.read_csv(config.HPA_ATLAS_ZIP, sep='\t', usecols=['Gene', 'Ensembl', 'RNA single cell type specificity',
                                                             'RNA single cell type specific nCPM'])
    a.columns = ['symbol', 'ens', 'category', 'specific_ncpm']
    a = a[a['category'].isin(ENRICHED)].set_index('ens')
    a['types'] = a['specific_ncpm'].map(
        lambda s: {t.split(': ')[0].lower(): float(t.split(': ')[1]) for t in s.split(';')})
    return a


def panel_table(panels, atlas=None):
    """panels: {name: [HPA cell types, lower case]} -> one row per (panel, gene)."""
    a = load_specificity() if atlas is None else atlas
    panel_of_type = {t: p for p, ts in panels.items() for t in ts}
    rows = [dict(panel=p, ens=e, symbol=r.symbol, hpa_category=r.category, hpa_types=';'.join(r.types),
                 ncpm_in=max(v for t, v in r.types.items() if panel_of_type.get(t) == p))
            for r in a.itertuples() for e in [r.Index]
            for p in sorted({panel_of_type[t] for t in r.types if t in panel_of_type})]
    return pd.DataFrame(rows)


def load_atlas():
    cols = {'Gene': 'symbol', 'Ensembl': 'ens', 'RNA tissue specificity': 't_cat',
            'RNA tissue specificity score': 't_score', 'RNA tissue specific nTPM': 't_ntpm',
            'RNA tissue cell type enrichment': 'tct', 'RNA single cell type specificity': 's_cat',
            'RNA single cell type specificity score': 's_score', 'RNA single cell type specific nCPM': 's_ncpm'}
    a = pd.read_csv(config.HPA_ATLAS_ZIP, sep='\t', usecols=list(cols)).rename(columns=cols)
    split = lambda s, sep, kv: {t.split(kv)[0].lower() for t in s.split(sep) if t}
    a['organs'] = a['t_ntpm'].fillna('').map(lambda s: split(s, ';', ':'))
    a['types'] = a['s_ncpm'].fillna('').map(lambda s: split(s, ';', ': '))
    a['tct'] = a['tct'].fillna('').map(lambda s: [t for t in s.split(', ') if t])
    return a


def tct_panel_table(panels, cap=200, atlas=None):
    """panels: {name: [HPA single-cell types, lower case]} -> one row per (panel, gene), at most cap per panel."""
    a = load_atlas() if atlas is None else atlas
    rows = []
    for p, types in panels.items():
        if not types:
            continue
        insc = a['types'].map(lambda s: bool(s & set(types)))
        sc_enr = insc & a['s_cat'].isin(ENRICHED)
        if p in TCT:
            labs, organ = TCT[p]
            base = a['tct'].map(lambda l: any(x == y or (y.startswith('- ') and x.endswith(y)) for x in l for y in labs))
            inorg = a['organs'].map(lambda s: organ in s) if organ else pd.Series(False, index=a.index)
            enr = sc_enr | (base & inorg & a['t_cat'].isin(TISSUE_ENRICHED))
            keep = enr | (base & (inorg | insc))
        else:
            inorg = pd.Series(False, index=a.index)
            enr = keep = sc_enr
        score = pd.concat([a['t_score'].where(inorg), a['s_score'].where(insc)], axis=1).max(axis=1)
        t = a[keep].assign(panel=p, tier=enr[keep].map({True: 'enriched', False: 'enhanced'}), score=score[keep])
        t = t.sort_values(['tier', 'score'], ascending=[False, False], kind='stable').head(cap)
        rows.append(t[['panel', 'ens', 'symbol', 'tier', 'score']])
    return pd.concat(rows, ignore_index=True)


if __name__ == '__main__':
    t = panel_table({'hepatocyte': ['hepatocytes'], 'acinar': ['pancreatic acinar cells']})
    assert {'APOH', 'FGB'} <= set(t.loc[t.panel == 'hepatocyte', 'symbol'])
    assert 'CTRB2' in set(t.loc[t.panel == 'acinar', 'symbol'])
    assert set(t['hpa_category']) <= set(ENRICHED)
    u = tct_panel_table({'hepatocyte': ['hepatocytes'], 'pancreatic acinar': ['pancreatic acinar cells'],
                         'plasma cell': ['plasma cells'], 'neuron': ['brain excitatory neurons']})
    n = u.groupby('panel').size()
    assert n.max() <= 200
    assert set(t.loc[t.panel == 'acinar', 'symbol']) <= set(u.loc[u.panel == 'pancreatic acinar', 'symbol'])
    h = u[u.panel == 'hepatocyte']
    assert (h['tier'] == 'enriched').all() and {'APOH', 'FGB'} <= set(h['symbol'])
    print(n)
