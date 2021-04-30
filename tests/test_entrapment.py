import os

import numpy as np

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
    fasta_host = os.path.join(fixtures_loc, "eight_sequences.fasta")
    fasta_trap = os.path.join(fixtures_loc, "UP000321567.fasta")
    entrapment.get_nearest_neighbor_proteins(fasta_host, fasta_trap)
    assert False
