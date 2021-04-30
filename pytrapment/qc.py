"""Module to perform QC on the xiRT performance."""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pyteomics import electrochem

sns.set(context="notebook", style="white", palette="deep", font_scale=1)


def compute_features(peptides_df):
    """
    Compute features to evaluate database.

    Parameters:
        peptides_df : df
            dataframe with peptides.

    Returns:
        None.
    """
    # features
    seq = "sequence"
    peptides_df["length"] = peptides_df["sequence"].apply(len)

    # amino acid counts
    peptides_df["KR"] = peptides_df[seq].str.count("K") + peptides_df[seq].str.count("R")
    peptides_df["aromatic"] = peptides_df[seq].str.count("F") + peptides_df[seq].str.count("W") + \
        peptides_df[seq].str.count("Y")
    peptides_df["acids"] = peptides_df[seq].str.count("D") + peptides_df[seq].str.count("E")
    peptides_df["aliphatic"] = peptides_df[seq].str.count("A") + peptides_df[seq].str.count("I") + \
        peptides_df[seq].str.count("L") + peptides_df[seq].str.count("M") + \
        peptides_df[seq].str.count("V")
    peptides_df["HGP"] = peptides_df[seq].str.count("G") + peptides_df[seq].str.count("P") + \
        peptides_df[seq].str.count("H")

    # sequence properties
    bp_seqs = peptides_df["sequence"].values
    peptides_df["isoelectric_point"] = [electrochem.pi(x) for x in bp_seqs]
    peptides_df["gravy"] = [electrochem.gravy(x) for x in bp_seqs]

    # transform shape
    df_melt = peptides_df.melt(id_vars=["Type", "OS", "sequence"],
                               value_vars=["length", "KR", "aromatic", "acids", "isoelectric_point",
                                           "gravy"])
    df_melt["log_value"] = np.log2(df_melt["value"])

    df_peptides = df_melt[df_melt["Type"] == "Peptide"]
    df_proteins = df_melt[df_melt["Type"] == "Protein"]

    # peptides
    g = sns.FacetGrid(df_peptides, col='variable', sharey=False, col_wrap=3)
    g.map_dataframe(sns.boxplot, x='OS', y='value', color="darkgrey")
    g.map_dataframe(sns.boxplot, x='OS', y='value', color="darkgrey")
    g.axes[4].set(xlabel="Peptide Properties")
    g.axes[0].set(ylabel="counts / value")
    g.axes[3].set(ylabel="counts / value")
    plt.tight_layout()

    return df_peptides, df_proteins
