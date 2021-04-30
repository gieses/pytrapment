"""Module to perform QC on the xiRT performance."""
import numpy as np
import pandas as pd
from pyteomics import fasta
from pyteomics import parser
from scipy.spatial import distance


def compute_composition_df(seq_df):
    """
    Compute the composition matrix for all proteins.

    Args:
        seq_df:  df, dataframe with sequences

    Returns:
        df, with the composition of the proteins
    """
    # get composition table
    df_seq_comp = pd.DataFrame(
        list(seq_df["sequence"].apply(parser.amino_acid_composition).values)) * 1.0

    # add column with 0s for amino acids that didnt occur in the protein fasta file
    for i in parser.std_amino_acids:
        if i not in df_seq_comp.columns:
            df_seq_comp[i] = 0

    df_seq_comp = df_seq_comp.fillna(0.0)
    df_seq_comp.index = seq_df.index
    return df_seq_comp


def get_nearest_neighbor_proteins(fasta_host, fasta_trap):
    """
    Retrieve the nearest neighbors for all proteins in the host fasta.

    Args:
        fasta_host:
        fasta_trap:

    Returns:
        df, dataframe with proteins for host and entrapment database.
    """
    # make nearest neighbor thing
    # ger protein table
    df_prot_host = fasta2dataframe(fasta_host)
    df_comp_host = compute_composition_df(df_prot_host)
    df_peptides_host = digest_protein_df(df_prot_host)

    df_prot_trap = fasta2dataframe(fasta_trap)
    df_comp_trap = compute_composition_df(df_prot_trap)
    df_peptides_trap = digest_protein_df(df_prot_trap)

    # perform the filtering
    df_comp_trap, df_prot_trap = filter_trap_fasta(df_prot_trap, df_comp_trap,
                                                   df_peptides_trap, df_peptides_host)

    # get best human protein matching by euclidean distance
    neighbor = []
    distances = np.zeros(df_comp_host.shape[0])
    for ii, row in enumerate(df_comp_host.iterrows()):
        # compute the distance from the query (current row) to all other proteins in the
        # trap database
        ary = distance.cdist(df_comp_trap, pd.DataFrame(row[1]).transpose(), metric='euclidean')
        # retrieve minimal disance entry here and use the index as neighbor and include
        # it to the fasta df later
        neighbor.append(df_comp_trap[ary == ary.min()].index.values[0])
        distances[ii] = ary.min()
    # print identifier for id mapping
    # neighbors = [i.split("|")[1] for i in np.ravel(neighbor)]
    fasta_df_entrapment = df_prot_trap.loc[neighbor]
    # store seed-neighbor pairs
    fasta_df_entrapment["host_seed"] = df_comp_host.index

    final_fasta_df = pd.concat([df_prot_host, fasta_df_entrapment])
    final_fasta_df["db_type"] = ["host"] * len(df_prot_host) + ["trap"] * len(fasta_df_entrapment)
    return final_fasta_df


def filter_trap_fasta(df_prot_trap, df_comp_trap, df_peptides_trap, df_peptides_host):
    """
    Remove proteins with peptides that also occur in the host database.

    Args:
        df_comp_trap: df, protein entries from the entrapment database (composition)
        df_peptides_host: df, peptides from the host fasta
        df_peptides_trap: df, peptides from the entrapment fasta
        df_prot_trap: df, proteins from the entrapment database

    Returns:
        (df_comp_trap, df_prot_trap), returns a tuple of valid (unique) trapment ids.
    """
    # make sure I/L witht he same mass doesnt mess with overlapping peptides
    df_peptides_host["sequence"] = df_peptides_host["sequence"].str.replace("I", "L")
    df_peptides_trap["sequence"] = df_peptides_trap["sequence"].str.replace("I", "L")

    df_peptides_host = df_peptides_host.set_index("sequence")
    df_peptides_trap = df_peptides_trap.set_index("sequence")
    df_joined = df_peptides_trap.join(df_peptides_host, rsuffix="_host", lsuffix="_trap",
                                      how="left")
    blacklist_proteins_trap = df_joined.dropna(subset=["protein_host"])["protein_trap"].unique()
    # drop proteins from dfs
    df_prot_trap = df_prot_trap.drop(blacklist_proteins_trap)
    df_comp_trap = df_comp_trap.drop(blacklist_proteins_trap)
    return df_comp_trap, df_prot_trap


def fasta2dataframe(FASTA):
    """
    Convert the entries in the FASTA file to a dataframe with ID, sequence and Type as column.

    Parameters
    FASTA : str
        Location of the FASTA file..

    Returns
        dataframe
    """
    # store proteins and ids here
    unique_proteins = []
    unique_ids = []

    with open(FASTA, mode='rt') as ffile:
        for description, sequence in fasta.FASTA(ffile):
            unique_proteins.append(sequence)
            unique_ids.append(description)

    # convert to dataframe
    df = pd.DataFrame(unique_proteins)
    df.columns = ["sequence"]
    df["Type"] = "Protein"
    df.index = unique_ids

    # make sure ordering is lost
    df = df.sample(frac=1, random_state=42)
    return df


def generete_peptides_proteins(FASTA, limit=-1):
    """
    Generate a sets of peptides and proteins for the two species in the FASTA.

    Parameters
    ----------
    FASTA : str
        Location of the FASTA file..

    Returns
    -------
    None.


    unique_peptides_HS = set()
    unique_peptides_EC = set()

    unique_proteins_HS = set()
    unique_proteins_EC = set()
    ec = 0
    hs = 0

    with open(FASTA, mode='rt') as ffile:
        for description, sequence in fasta.FASTA(ffile):
            if limit != -1:
                if hs == limit:
                    print("Found enough proteins!)")
                    break

            # digest
            new_peptides = parser.cleave(sequence, 'trypsin')

            unique_peptides_EC.update(new_peptides)
            unique_proteins_EC.update([sequence])

            unique_peptides_HS.update(new_peptides)
            unique_proteins_HS.update([sequence])
            hs += 1

            else:
                print(description, sequence)

    # create id for df
    lengths = [len(unique_peptides_EC), len(unique_proteins_EC),
               len(unique_peptides_HS), len(unique_proteins_HS)]
    ids = ["Peptide", "Protein", "Peptide", "Protein"]
    id_column = np.repeat(ids, lengths)
    os_column = np.repeat(["E. coli", "E. coli", "H. sapiens", "H. sapiens"],
                          lengths)
    sequences = np.concatenate([list(unique_peptides_EC), list(unique_proteins_EC),
                                list(unique_peptides_HS), list(unique_proteins_HS)])
    df = pd.DataFrame({"Type": id_column, "OS": os_column,
                       "sequence": sequences})
    return df
    """
    return 1


def digest_protein_df(df_fasta, rule="trypsin", min_length=6):
    """
    Digest a dataframe of proteins into a dataframe with unique  peptides.

    Args:
        df_fasta: df, dataframe with protein, sequence columns
        rule: str, pyteomics string identifier
        min_length: int, minimal length for a peptide

    Returns:
        peptide_df with <protein_name:peptide> entries.
    """
    # create all peptide sequences first
    cleaved = df_fasta["sequence"].apply(parser.cleave, args=(rule,)).explode()
    return cleaved[cleaved.apply(len) >= min_length].rename_axis("protein").reset_index()
