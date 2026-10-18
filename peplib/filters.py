import os
import re
import glob
from typing import Dict, Tuple
import pandas as pd

CANONICAL_REGEX = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
AMBIGUOUS_CHARS = set('XBZJ')
RARE_CHARS = set('UO')

REJECTION_CAUSES = [
    'manquante',
    'vide',
    'lettres_ambigues',
    'lettres_rares',
    'autres_caracteres',
    'trop_courte',
    'trop_longue',
]

def classify_sequence(seq: str, min_len: int = 10, max_len: int = 50) -> str | None:
    if seq == '':
        return 'vide'
    if any(c in AMBIGUOUS_CHARS for c in seq):
        return 'lettres_ambigues'
    if any(c in RARE_CHARS for c in seq):
        return 'lettres_rares'
    if not CANONICAL_REGEX.match(seq):
        return 'autres_caracteres'
    if len(seq) < min_len:
        return 'trop_courte'
    if len(seq) > max_len:
        return 'trop_longue'
    return None

def filter_dataframe(df: pd.DataFrame, seq_col: str = 'sequence', min_len: int = 10, max_len: int = 50) -> Tuple[pd.DataFrame, Dict[str, int]]:
    total_raw = len(df)
    rejections: Dict[str, int] = {cause: 0 for cause in REJECTION_CAUSES}
    na_mask = df[seq_col].isna()
    rejections['manquante'] = int(na_mask.sum())
    df_non_na = df.dropna(subset=[seq_col]).copy()
    df_non_na[seq_col] = df_non_na[seq_col].astype(str).str.strip().str.upper()
    kept_indices = []
    for idx, seq in df_non_na[seq_col].items():
        cause = classify_sequence(seq, min_len=min_len, max_len=max_len)
        if cause is None:
            kept_indices.append(idx)
        else:
            rejections[cause] += 1
    df_kept = df_non_na.loc[kept_indices].copy().reset_index(drop=True)
    total_rejected = sum(rejections.values())
    if len(df_kept) + total_rejected != total_raw:
        raise ValueError(f'Bilan invalide: {len(df_kept)} + {total_rejected} != {total_raw}')
    return df_kept, rejections

def standardize_raw_tsv(raw_path: str, source_name: str, out_csv_path: str | None = None) -> pd.DataFrame:
    df_raw = pd.read_csv(raw_path, sep='\t', dtype=str)
    id_col = 'Entry' if 'Entry' in df_raw.columns else df_raw.columns[0]
    seq_col = 'Sequence' if 'Sequence' in df_raw.columns else 'sequence'
    df_t11 = pd.DataFrame({'sequence': df_raw[seq_col], 'source': source_name, 'id': df_raw[id_col]})
    if out_csv_path:
        os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)
        df_t11.to_csv(out_csv_path, index=False)
    return df_t11

def run_filtering_pipeline(raw_dir: str = 'data/raw', processed_dir: str = 'data/processed') -> pd.DataFrame:
    os.makedirs(processed_dir, exist_ok=True)
    raw_std_dir = os.path.join(processed_dir, 'raw_standardized')
    os.makedirs(raw_std_dir, exist_ok=True)
    prefixes = {'uniprot': 'uniprot_reviewed_', 'amp': 'amp_database_', 'bactibase': 'bactibase_bacteriocins_'}
    rejection_rows = []
    all_cleaned_dfs = []
    for source_name, prefix in prefixes.items():
        matches = sorted(glob.glob(os.path.join(raw_dir, f'{prefix}*.tsv')))
        if not matches:
            raise FileNotFoundError(f'Introuvable: {prefix}')
        raw_path = matches[-1]
        t11_csv = os.path.join(raw_std_dir, f'{source_name}_t11.csv')
        df_t11 = standardize_raw_tsv(raw_path, source_name, out_csv_path=t11_csv)
        df_kept, rejections = filter_dataframe(df_t11, seq_col='sequence', min_len=10, max_len=50)
        df_kept.to_csv(os.path.join(processed_dir, f'{source_name}.csv'), index=False)
        all_cleaned_dfs.append(df_kept)
        rejection_rows.append({'source': source_name, 'total_brut': len(df_t11), 'conservees': len(df_kept), 'total_rejetees': sum(rejections.values()), **rejections})
    df_rejections = pd.DataFrame(rejection_rows)
    df_rejections.to_csv(os.path.join(processed_dir, 'rejection_table_t14.csv'), index=False)
    pd.concat(all_cleaned_dfs, ignore_index=True).to_csv(os.path.join(processed_dir, 'corpus_filtered_all.csv'), index=False)
    return df_rejections

if __name__ == '__main__':
    print(run_filtering_pipeline().to_string(index=False))
