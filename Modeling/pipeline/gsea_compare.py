import glob
import json
import os

import numpy as np
import pandas as pd

import config

REF_DIR = config.BENCHMARK_DIR / 'disease_reference'

FDR = config.MODELING_PARAMS['gsea_fdr_thr']
OUT_DIR = config.BENCHMARK_DIR / 'gsea_compare'


def _load_dir(d):
    out = {}
    for f in sorted(glob.glob(str(d / 'gsea_result_*.csv'))):
        ph = os.path.basename(f)[len('gsea_result_'):-4].strip()
        df = pd.read_csv(f)
        out[ph] = df[df['FDR q-val'] < FDR].copy()
    return out


def load_sets():
    """Three FDR-filtered term-set dicts keyed by canonical phenotype name:
    no_filter (rare-excluded normative), with_rare (rare-included), deseq2 (DESeq2 rank GSEA)."""
    return (_load_dir(config.GSEA_DIR / 'no_filter'),
            _load_dir(config.GSEA_DIR / 'with_rare'),
            _load_dir(config.DESEQ2_GSEA_DIR))


def _lib(term):
    return term.split('__', 1)[0]


def term_nes(df):
    """{Term: NES} keeping the strongest-|NES| row per term."""
    d = (df.assign(_a=df['NES'].abs())
           .sort_values('_a', ascending=False)
           .drop_duplicates('Term'))
    return dict(zip(d['Term'], d['NES']))


def overlap(a_df, b_df):
    """Term-set overlap between two GSEA result frames. Returns a stats dict; 'right'
    is the second set, so right_only = terms significant in b but not a."""
    a, b = term_nes(a_df), term_nes(b_df)
    sa, sb = set(a), set(b)
    common = sa & sb
    agree = sum(1 for t in common if np.sign(a[t]) == np.sign(b[t]))
    return {'n_left': len(sa), 'n_right': len(sb), 'n_common': len(common),
            'left_only': len(sa - sb), 'right_only': len(sb - sa),
            'jaccard': len(common) / len(sa | sb) if (sa | sb) else np.nan,
            'sign_agree': agree / len(common) if common else np.nan}


def diff_terms(a_df, b_df, which):
    """Terms unique to one side, as a frame with NES + FDR + library + lead genes.
    which='right_only' -> in b not a; 'left_only' -> in a not b."""
    a, b = set(term_nes(a_df)), set(term_nes(b_df))
    src, keep = (b_df, b - a) if which == 'right_only' else (a_df, a - b)
    sub = (src[src['Term'].isin(keep)]
           .assign(_a=src.loc[src['Term'].isin(keep), 'NES'].abs())
           .sort_values('_a', ascending=False)
           .drop_duplicates('Term'))
    sub = sub.assign(library=sub['Term'].map(_lib))
    cols = [c for c in ['Term', 'library', 'NES', 'FDR q-val', 'NOM p-val', 'Lead_genes']
            if c in sub.columns]
    return sub[cols].reset_index(drop=True)


COMPARISONS = {
    'rare_vs_norare': ('no_filter', 'with_rare'),
    'deseq2_vs_norare': ('deseq2', 'no_filter'),
    'deseq2_vs_withrare': ('deseq2', 'with_rare'),
}


def compare_all(save=True):
    """Per-disease term-overlap stats for all three comparisons + per-disease diff term
    lists. Returns (stats_df, diffs) where diffs[comparison][phenotype][which] = frame."""
    nf, wr, dq = load_sets()
    src = {'no_filter': nf, 'with_rare': wr, 'deseq2': dq}
    rows, diffs = [], {}
    for name, (lkey, rkey) in COMPARISONS.items():
        L, R = src[lkey], src[rkey]
        diffs[name] = {}
        for ph in sorted(set(L) & set(R)):
            st = overlap(L[ph], R[ph])
            st.update({'comparison': name, 'left': lkey, 'right': rkey, 'phenotype': ph})
            rows.append(st)
            diffs[name][ph] = {'right_only': diff_terms(L[ph], R[ph], 'right_only'),
                               'left_only': diff_terms(L[ph], R[ph], 'left_only')}
    stats = pd.DataFrame(rows)[['comparison', 'left', 'right', 'phenotype', 'n_left',
                                'n_right', 'n_common', 'left_only', 'right_only',
                                'jaccard', 'sign_agree']]
    if save:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        stats.to_csv(OUT_DIR / 'overlap_stats.csv', index=False)
        for name, phd in diffs.items():
            for ph, d in phd.items():
                stem = ph.replace('/', '_')
                for which, frame in d.items():
                    frame.to_csv(OUT_DIR / f'{name}__{which}__{stem}.csv', index=False)
    return stats, diffs


def rare_gene_symbols():
    """HGNC symbols of the 559 rare-branch genes, for attributing new-in-with_rare terms to
    the rare covariate GLM. Uses adata.var GeneName (same ensg->symbol map GSEA ranked on),
    not the DESeq2 csv whose gene universe excludes most very-low-expression rare genes."""
    import anndata as ad
    ts = pd.read_csv(config.ENGINE_DIR / 'training_summary.csv')
    rare = ts.loc[ts['branch'].astype(str).str.startswith('rare'), 'gene'].tolist()
    var = ad.read_h5ad(config.PATHS['merged_qc'], backed='r').var
    gname = var['GeneName']
    vn = set(var.index)
    return {gname[g] for g in rare if g in vn and pd.notna(gname[g])}


def load_reference():
    """{phenotype: set(reference gene symbols)} from the Open Targets reference DB."""
    out = {}
    for f in glob.glob(str(REF_DIR / '*.json')):
        r = json.load(open(f))
        out[r['phenotype']] = {g for g, _ in r['genes']}
    return out


def _leads(s):
    return set() if (isinstance(s, float) and np.isnan(s)) else set(str(s).split(';'))


def _term_leads(df):
    """{Term: set(lead genes)} keeping the strongest-|NES| row per term."""
    d = (df.assign(_a=df['NES'].abs()).sort_values('_a', ascending=False)
           .drop_duplicates('Term'))
    return {t: _leads(l) for t, l in zip(d['Term'], d.get('Lead_genes', pd.Series(index=d.index)))}


def deseq2_coverage(save=True):
    """For each phenotype, how well DESeq2-rank GSEA recovers the DB-supported significant
    pathways that the normative model finds. A normative term is db_supported if its lead
    genes intersect the Open Targets reference; captured if DESeq2 also calls it significant.
    Computed against both normative variants (no_filter, with_rare). Returns a stats frame."""
    nf, wr, dq = load_sets()
    ref = load_reference()
    rows = []
    for nkey, norm in [('no_filter', nf), ('with_rare', wr)]:
        for ph in sorted(set(dq) & set(norm)):
            dref = ref.get(ph, set())
            nl = _term_leads(norm[ph])
            nn = term_nes(norm[ph])
            dn = term_nes(dq[ph])
            dset = set(dn)
            db_terms = {t for t, lead in nl.items() if lead & dref} if dref else set()
            captured = {t for t in db_terms if t in dset}
            common = set(nn) & dset
            disagree = sum(1 for t in common if np.sign(nn[t]) != np.sign(dn[t]))
            rows.append({
                'normative': nkey, 'phenotype': ph,
                'n_norm_sig': len(nn), 'n_deseq2_sig': len(dn),
                'n_norm_db': len(db_terms), 'n_db_captured': len(captured),
                'db_capture_rate': round(len(captured) / len(db_terms), 3) if db_terms else np.nan,
                'n_common': len(common), 'n_sign_disagree': disagree,
                'sign_disagree_rate': round(disagree / len(common), 3) if common else np.nan,
                'has_ot_ref': len(dref) > 0})
    cov = pd.DataFrame(rows)
    if save:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        cov.to_csv(OUT_DIR / 'deseq2_coverage.csv', index=False)
    return cov


def validate_rare_novel(diffs, save=True):
    """Cross-validate the new-in-with_rare terms (rare_vs_norare right_only) against the
    reference DB and the rare-branch gene set. Returns (long, summary). A new term is
    db_supported if its lead genes hit the disease reference; rare_driven if they hit the
    rare-branch genes."""
    ref = load_reference()
    rare_syms = rare_gene_symbols()
    rows = []
    for ph, d in diffs['rare_vs_norare'].items():
        dref = ref.get(ph, set())
        for _, r in d['right_only'].iterrows():
            lead = _leads(r.get('Lead_genes'))
            db_hit = sorted(lead & dref)
            rare_hit = sorted(lead & rare_syms)
            rows.append({'phenotype': ph, 'Term': r['Term'], 'library': r.get('library'),
                         'NES': r['NES'], 'FDR q-val': r.get('FDR q-val'), 'n_lead': len(lead),
                         'n_db': len(db_hit), 'db_hits': ';'.join(db_hit),
                         'n_rare': len(rare_hit), 'rare_hits': ';'.join(rare_hit),
                         'db_supported': len(db_hit) > 0, 'rare_driven': len(rare_hit) > 0})
    long = pd.DataFrame(rows)
    g = long.groupby('phenotype')
    summary = pd.DataFrame({
        'n_new_terms': g.size(),
        'n_db_supported': g['db_supported'].sum(),
        'n_rare_driven': g['rare_driven'].sum(),
        'n_db_and_rare': g.apply(lambda x: (x['db_supported'] & x['rare_driven']).sum()),
        'has_ot_ref': g['phenotype'].first().map(lambda p: len(ref.get(p, set())) > 0),
    })
    summary['frac_db_supported'] = (summary['n_db_supported'] / summary['n_new_terms']).round(3)
    summary = summary.reset_index().sort_values('n_new_terms', ascending=False)
    if save:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        long.to_csv(OUT_DIR / 'rare_novel_validated.csv', index=False)
        summary.to_csv(OUT_DIR / 'rare_novel_summary.csv', index=False)
    return long, summary
