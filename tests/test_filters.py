import os
import numpy as np
import pandas as pd
from peplib.filters import CANONICAL_REGEX, classify_sequence, filter_dataframe, run_filtering_pipeline

def test_classify_sequence():
    assert classify_sequence('') == 'vide'
    assert classify_sequence('ACDEFGHIKX') == 'lettres_ambigues'
    assert classify_sequence('ACDEFGHIKB') == 'lettres_ambigues'
    assert classify_sequence('ACDEFGHIKU') == 'lettres_rares'
    assert classify_sequence('ACDEFGHIKO') == 'lettres_rares'
    assert classify_sequence('ACDEF12345') == 'autres_caracteres'
    assert classify_sequence('ACDEF-GHIK') == 'autres_caracteres'
    assert classify_sequence('ACDEFGHIK') == 'trop_courte'
    assert classify_sequence('A' * 51) == 'trop_longue'
    assert classify_sequence('ACDEFGHIKL') is None
    assert classify_sequence('A' * 50) is None

def test_filter_dataframe_all_causes():
    data = pd.DataFrame({
        'sequence': [np.nan, None, '', '   ', 'ACDEFGHIKX', 'ACDEFGHIKU', 'ACDEF12345', 'ACDEFGHIK', 'A' * 51, 'acdefghikl', '  ACDEFGHIKLMNP  '],
        'source': ['test'] * 11,
        'id': [f'id_{i}' for i in range(11)],
    })
    df_kept, rejections = filter_dataframe(data, seq_col='sequence', min_len=10, max_len=50)
    assert len(df_kept) == 2
    assert list(df_kept['sequence']) == ['ACDEFGHIKL', 'ACDEFGHIKLMNP']
    assert rejections['manquante'] == 2
    assert rejections['vide'] == 2
    assert rejections['lettres_ambigues'] == 1
    assert rejections['lettres_rares'] == 1
    assert rejections['autres_caracteres'] == 1
    assert rejections['trop_courte'] == 1
    assert rejections['trop_longue'] == 1
    assert len(df_kept) + sum(rejections.values()) == len(data)
    assert all(CANONICAL_REGEX.match(s) for s in df_kept['sequence'])

def test_standardize_and_pipeline(tmp_path):
    raw_dir = tmp_path / 'raw'
    proc_dir = tmp_path / 'processed'
    raw_dir.mkdir()
    for prefix in ['uniprot_reviewed_', 'amp_database_', 'bactibase_bacteriocins_']:
        (raw_dir / f'{prefix}20261017.tsv').write_text('Entry\tSequence\tLength\nP1\tACDEFGHIKL\t10\nP2\tACDEFX\t6\n', encoding='utf-8')
    df_rej = run_filtering_pipeline(raw_dir=str(raw_dir), processed_dir=str(proc_dir))
    assert len(df_rej) == 3
    assert (df_rej['total_brut'] == 2).all()
    assert (df_rej['conservees'] == 1).all()
    assert (df_rej['total_rejetees'] == 1).all()
    assert os.path.exists(proc_dir / 'uniprot.csv')
    assert os.path.exists(proc_dir / 'rejection_table_t14.csv')
