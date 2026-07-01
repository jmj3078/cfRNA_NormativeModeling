import datetime
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

import config

OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
TOP_N = 300

# phenotype (GSEA file stem) -> Open Targets disease search query.
# Heterogeneous / undefined-cohort phenotypes map to None (no OT gene reference; literature-only).
PHENO_QUERY = {
    'CAD_HF+ (Ward)': 'coronary artery disease',
    'CAD_HF- (Ward)': 'coronary artery disease',
    'Colorectal Cancer (Chen)': 'colorectal carcinoma',
    'Esophagus Cancer (Chen)': 'esophageal carcinoma',
    'HIV (Chang)': 'HIV infection',
    'HIV + Tuberculosis (Chang)': 'tuberculosis',
    'Tuberculosis (Chang)': 'tuberculosis',
    'Liver Cancer (Chen)': 'hepatocellular carcinoma',
    'Liver Cancer (Roskams-Hieter)': 'hepatocellular carcinoma',
    'Lung Cancer (Chen)': 'lung carcinoma',
    'ME_CFS (Gardella)': 'chronic fatigue syndrome',
    'MGUS (Roskams-Hieter)': 'monoclonal gammopathy of undetermined significance',
    'MM (Roskams-Hieter)': 'multiple myeloma',
    'Pancreatic Cancer (Moore)': 'pancreatic carcinoma',
    'Pancreatitis (Moore)': 'pancreatitis',
    'Pre-eclampsia (Moufarrej)': 'pre-eclampsia',
    'Stomach Cancer (Chen)': 'gastric carcinoma',
    'Other Cancer (Moore)': 'cancer',
    'ICI-treated Cancer (Raissadati)': None,
    'ICI-m (Raissadati)': 'myocarditis',
}


def gql(query, variables=None):
    body = json.dumps({'query': query, 'variables': variables or {}}).encode()
    req = urllib.request.Request(OPEN_TARGETS, data=body,
                                 headers={'Content-Type': 'application/json'})
    err = None
    for _ in range(4):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=30).read())
        except Exception as e:
            err = e
            time.sleep(2)
    raise err


def resolve(q):
    d = gql('query($q:String!){search(queryString:$q,entityNames:["disease"],'
            'page:{index:0,size:1}){hits{id name}}}', {'q': q})
    h = d['data']['search']['hits']
    return (h[0]['id'], h[0]['name']) if h else (None, None)


def assoc_targets(efo, size=TOP_N):
    q = ('query($id:String!,$p:Pagination!){disease(efoId:$id){name '
         'associatedTargets(page:$p,orderByScore:"score"){count '
         'rows{target{approvedSymbol} score}}}}')
    d = gql(q, {'id': efo, 'p': {'index': 0, 'size': size}})
    dd = d['data']['disease']
    if dd is None:
        return None, []
    at = dd['associatedTargets']
    return at['count'], [(r['target']['approvedSymbol'], round(r['score'], 4)) for r in at['rows']]


def main():
    outdir = config.BENCHMARK_DIR / 'disease_reference'
    outdir.mkdir(parents=True, exist_ok=True)
    today = datetime.date.today().isoformat()
    for ph, q in PHENO_QUERY.items():
        stem = ph.replace('/', '_')
        if q is None:
            rec = {'phenotype': ph, 'query': None, 'efo': None, 'ot_disease': None,
                   'n_assoc': 0, 'genes': [], 'note': 'literature-only (heterogeneous)',
                   'source': 'Open Targets Platform GraphQL v4', 'retrieved': today}
        else:
            efo, name = resolve(q)
            cnt, genes = assoc_targets(efo) if efo else (0, [])
            rec = {'phenotype': ph, 'query': q, 'efo': efo, 'ot_disease': name,
                   'n_assoc': cnt, 'top_n': TOP_N, 'genes': genes,
                   'source': 'Open Targets Platform GraphQL v4', 'retrieved': today}
        json.dump(rec, open(outdir / f'{stem}.json', 'w'), indent=1)
        print(f'{ph:32s} {str(rec["efo"]):16s} n_genes={len(rec["genes"])}')
        time.sleep(0.3)
    print(f'\nsaved {len(PHENO_QUERY)} reference files to {outdir}')


if __name__ == '__main__':
    main()
