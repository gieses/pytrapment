import os

from pytrapment import entrapment, qc

# fixtures is used to store files used for the tests only
fixtures_loc = os.path.join(os.path.dirname(__file__), 'fixtures')


def test_compute_peptide_features():
    # poor test for completeness not accurate features
    fasta_file = os.path.join(fixtures_loc, "eight_sequences.fasta")
    proteins_df = entrapment.fasta2dataframe(fasta_file)

    peptides_df = entrapment.digest_protein_df(proteins_df)
    features_df = qc.compute_sequence_features(peptides_df)

    assert features_df.shape[1] == 10


def test_qc_peptides(tmpdir):
    fasta_host = os.path.join(fixtures_loc, "eight_sequences.fasta")
    fasta_trap = os.path.join(fixtures_loc, "multiple_sequences.fasta")

    final_df = entrapment.get_nearest_neighbor_proteins(fasta_host, fasta_trap)

    host_peptides = entrapment.digest_protein_df(final_df[final_df["db_type"] == "host"])
    trap_peptides = entrapment.digest_protein_df(final_df[final_df["db_type"] == "trap"])

    features_df_host = qc.compute_sequence_features(host_peptides)
    features_df_host["Type"] = "host"

    features_df_trap = qc.compute_sequence_features(trap_peptides)
    features_df_trap["Type"] = "trap"

    p = tmpdir.mkdir("pytrament_test")
    qc.qc_peptides(features_df_host, features_df_trap, p)
    assert True
