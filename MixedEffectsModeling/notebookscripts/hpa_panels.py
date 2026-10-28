"""Cell-type marker panels from the official HPA 25.1 single-cell type specificity categories.

A gene belongs to a panel when HPA calls it 'Cell type enriched' or 'Group enriched' and one of
its enriched cell types is in the panel. Shared by 6_disease_tissue_signal and
7_pancreatic_cancer_cohort_deepdive so both use one marker definition.
"""
import pandas as pd

import MixedEffectsModeling.config as config

ENRICHED = ('Cell type enriched', 'Group enriched')


def load_specificity():
    a = pd.read_csv(config.HPA_ATLAS_ZIP, sep='\t', usecols=['Gene', 'Ensembl', 'RNA single cell type specificity',
                                                             'RNA single cell type specific nCPM'])
    a = a[a['RNA single cell type specificity'].isin(ENRICHED)].set_index('Ensembl')
    a['types'] = a['RNA single cell type specific nCPM'].map(
        lambda s: {t.split(': ')[0].lower(): float(t.split(': ')[1]) for t in s.split(';')})
    return a


def panel_table(panels, atlas=None):
    """panels: {name: [HPA cell types, lower case]} -> one row per (panel, gene)."""
    a = load_specificity() if atlas is None else atlas
    panel_of_type = {t: p for p, ts in panels.items() for t in ts}
    rows = [dict(panel=p, ens=e, symbol=r.Gene, hpa_category=r._3, hpa_types=';'.join(r.types),
                 ncpm_in=max(v for t, v in r.types.items() if panel_of_type.get(t) == p))
            for r in a.itertuples() for e in [r.Index]
            for p in sorted({panel_of_type[t] for t in r.types if t in panel_of_type})]
    return pd.DataFrame(rows)


if __name__ == '__main__':
    t = panel_table({'hepatocyte': ['hepatocytes'], 'acinar': ['pancreatic acinar cells']})
    assert {'APOH', 'FGB'} <= set(t.loc[t.panel == 'hepatocyte', 'symbol'])
    assert 'CTRB2' in set(t.loc[t.panel == 'acinar', 'symbol'])
    print(t.groupby('panel').size())
