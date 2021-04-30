import os

import numpy as np
import pandas as pd

from pytrapment import entrapment

# fixtures is used to store files used for the tests only
fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_generate_proteins():
    fasta_file = os.path.join(fixtures_loc, "eight_sequences.fasta")
    proteins_df = entrapment.fasta2dataframe(fasta_file)
    assert proteins_df.shape == (8, 2)


def test_compute_composition_df():
    fasta_file = os.path.join(fixtures_loc, "mock_seq.fasta")
    proteins_df = entrapment.fasta2dataframe(fasta_file)
    mock_comp = entrapment.compute_composition_df(proteins_df)
    assert bool(np.all(mock_comp == 1)) is True
    assert mock_comp.shape[1] == 20


def test_get_nearest_neighbor_proteins():
    fasta_host = os.path.join(fixtures_loc, "one_sequence_host.fasta")
    fasta_trap = os.path.join(fixtures_loc, "two_sequences_trap.fasta")
    final_fasta = entrapment.get_nearest_neighbor_proteins(fasta_host, fasta_trap)

    assert np.all(final_fasta.index == ["mock1_host_R", "mock1_unique_Ks"])
    assert final_fasta.shape[0] == 2


def test_digest_protein_df():
    mock_df = pd.DataFrame()
    mock_df["protein"] = ["A"]
    mock_df["sequence"] = ["ELVISKLIVESR"]

    peptide_df = entrapment.digest_protein_df(mock_df)

    assert np.all(sorted(peptide_df["sequence"].tolist()) == ["ELVISK", "LIVESR"])


def test_filter_trap_fasta():
    host_df = pd.DataFrame()
    host_df["protein"] = ["A", "B"]
    host_df["sequence"] = ["ELVISKLIVESR", "PEPTIDERLEPTIDEK"]
    host_df = host_df.set_index("protein")

    df_peptides_host = entrapment.digest_protein_df(host_df)

    trap_df = pd.DataFrame()
    trap_df["protein"] = ["C", "D"]
    trap_df["sequence"] = ["ELVISKSVEN", "RANDEMPEPTIDE"]
    trap_df = trap_df.set_index("protein")
    df_comp_trap = entrapment.compute_composition_df(trap_df)
    df_peptides_trap = entrapment.digest_protein_df(trap_df)

    # perform the filtering
    df_comp_trap, df_prot_trap = \
        entrapment.filter_trap_fasta(trap_df, df_comp_trap, df_peptides_trap, df_peptides_host)

    assert df_comp_trap.shape[0] == 1
    assert df_prot_trap.index[0] == "D"
