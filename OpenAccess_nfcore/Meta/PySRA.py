import pandas as pd
from pysradb.sraweb import SRAweb
from tqdm import tqdm
import time

with open("./SRA_ids.txt") as f:
    srr_list = [line.strip() for line in f if line.strip()]
db = SRAweb()
all_results = []
chunk_size = 50  
for i in tqdm(range(0, len(srr_list), chunk_size), desc="pysradb Metadata Fetching"):
    chunk = srr_list[i:i + chunk_size]
    try:
        df_chunk = db.sra_metadata(chunk, detailed=True)
        
        if df_chunk is not None and not df_chunk.empty:
            all_results.append(df_chunk)
        
        time.sleep(1) 
        
    except Exception as e:
        print(f"\n[Error Occured]: Index {i} : {e}")
        

def make_columns_unique(df):
    new_cols = []
    col_count = {}
    for col in df.columns:
        if col in col_count:
            col_count[col] += 1
            new_cols.append(f"{col}_{col_count[col]}")
        else:
            col_count[col] = 0
            new_cols.append(col)
    df.columns = new_cols
    return df

sanitized_results = [make_columns_unique(df.copy()) for df in all_results]
final_df = pd.concat(sanitized_results, axis=0, ignore_index=True, sort=False)

print(f"{final_df.shape}")
final_df = final_df.sort_values(by=['bioproject','run_accession'])
final_df.drop_duplicates().to_csv("./meta_raw.tsv", sep='\t')
author = pd.read_table("/project/cfRNA_Disentaglement/Data/OpenAccess_nfcore/Meta/meta_author.tsv", index_col=0)
meta_with_author = final_df.merge(author, left_on='bioproject', right_on='BioProject', how='left')
meta_with_author.to_csv("./meta_raw_with_author.tsv", sep='\t', index=False)