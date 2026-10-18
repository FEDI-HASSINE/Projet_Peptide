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


def test_classifie_premiere_cause_gagne():
    """T15 — règle de première cause : l'ordre manquante→…→trop_longue gagne."""
    # Ambigue ET trop longue (51 X) -> ambigue, PAS trop_longue.
    assert classify_sequence('X' * 51) == 'lettres_ambigues'
    # Rare ET trop courte (9 U) -> rare, PAS trop_courte.
    assert classify_sequence('U' * 9) == 'lettres_rares'
    # Ambigue ET trop courte (XBZJ, 4 résidus) -> ambigue.
    assert classify_sequence('XBZJ') == 'lettres_ambigues'
    # Autres caractères ET trop courte -> autres_caracteres.
    assert classify_sequence('123') == 'autres_caracteres'
    # Rares ET trop longue (51 O) -> rares, PAS trop_longue.
    assert classify_sequence('O' * 51) == 'lettres_rares'


def test_toutes_lettres_ambigues_et_rares():
    """T15 — toutes les lettres ambigues (XBZJ) et rares (UO)."""
    # 'ACDEFGHIK' = 9 résidus ; +1 lettre = 10 => longueur valide,
    # seule la lettre spéciale cause le rejet.
    for lettre in 'XBZJ':
        assert classify_sequence(lettre + 'ACDEFGHIK') == 'lettres_ambigues'
    for lettre in 'UO':
        assert classify_sequence(lettre + 'ACDEFGHIK') == 'lettres_rares'


def test_bordures_longueur():
    """T15 — bornes exactes : 9 rejetée, 10 et 50 comprises, 51 rejetée."""
    assert classify_sequence('A' * 9) == 'trop_courte'
    assert classify_sequence('A' * 10) is None
    assert classify_sequence('A' * 50) is None
    assert classify_sequence('A' * 51) == 'trop_longue'


def test_rejet_par_sequence_et_invariant():
    """T15 — chaque séquence comptée dans UNE seule cause ; invariant."""
    data = pd.DataFrame({
        'sequence': [np.nan, '', '   ', 'X' + 'A' * 9, 'U' + 'A' * 9,
                     'ACDEF1', 'A' * 9, 'A' * 51,
                     'acdefghikl', '  ACDEFGHIKLM  '],
        'source': ['test'] * 10,
        'id': [f'id_{i}' for i in range(10)],
    })
    df_kept, rejections = filter_dataframe(data, seq_col='sequence',
                                           min_len=10, max_len=50)
    # Mapping cause par cause, chaque non-valide dans sa PREMIERE cause.
    attendu = {
        'manquante': 1,
        'vide': 2,               # '' et '   ' (vide après strip)
        'lettres_ambigues': 1,   # X+9A
        'lettres_rares': 1,      # U+9A
        'autres_caracteres': 1,  # ACDEF1 (autres AVANT trop_courte)
        'trop_courte': 1,        # 9 résidus
        'trop_longue': 1,        # 51 résidus
    }
    assert rejections == attendu
    assert sum(rejections.values()) == 8
    assert len(df_kept) == 2
    # Invariant du spec : conservés + rejetés = total brut.
    assert len(df_kept) + sum(rejections.values()) == len(data)
    # Les 2 valides conservées, normalisées upper/strip.
    assert sorted(df_kept['sequence']) == ['ACDEFGHIKL', 'ACDEFGHIKLM']


def test_nan_est_manquante_pas_autres():
    """T15 — piège du spec : dropna() AVANT astype(str)."""
    data = pd.DataFrame({
        'sequence': [np.nan, 'ACDEFGHIKL'],
        'source': ['test', 'test'],
        'id': ['a', 'b'],
    })
    df_kept, rejections = filter_dataframe(data, seq_col='sequence',
                                           min_len=10, max_len=50)
    assert rejections['manquante'] == 1
    # 'NAN' ne doit PAS passer par la regex dans autres_caracteres.
    assert rejections['autres_caracteres'] == 0
    assert rejections['vide'] == 0
    assert len(df_kept) == 1
