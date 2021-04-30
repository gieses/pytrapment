"""Module to perform QC on the xiRT performance."""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pyteomics import electrochem

sns.set(context="notebook", style="white", palette="deep", font_scale=1)


def compute_features(df):
    """
    Compute features to evaluate database

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # features
    df["length"] = df["sequence"].apply(len)

    # amino acid counts
    df["KR"] = df["sequence"].str.count("K") + df["sequence"].str.count("R")
    df["aromatic"] = df["sequence"].str.count("F") + \
                     df["sequence"].str.count("W") + df["sequence"].str.count("Y")
    df["acids"] = df["sequence"].str.count("D") + df["sequence"].str.count("E")
    df["aliphatic"] = df["sequence"].str.count("A") + df["sequence"].str.count("I") + \
                      df["sequence"].str.count("L") + df["sequence"].str.count("M") + df[
                          "sequence"].str.count("V")
    df["HGP"] = df["sequence"].str.count("G") + df["sequence"].str.count("P") + \
                df["sequence"].str.count("H")

    # sequence properties
    bp_seqs = df["sequence"].values
    df["isoelectric_point"] = [electrochem.pi(x) for x in bp_seqs]
    df["gravy"] = [electrochem.gravy(x) for x in bp_seqs]

    # transform shape
    df_melt = df.melt(id_vars=["Type", "OS", "sequence"],
                      value_vars=["length", "KR", "aromatic", "acids", "isoelectric_point",
                                  "gravy"])
    df_melt["log_value"] = np.log2(df_melt["value"])

    df_peptides = df_melt[df_melt["Type"] == "Peptide"]
    df_proteins = df_melt[df_melt["Type"] == "Protein"]
    return df_peptides, df_proteins

    # peptides
    g = sns.FacetGrid(df_peptides, col='variable', sharey=False, col_wrap=3)
    g.map_dataframe(custom_boxplot, x='OS', y='value', color="darkgrey")
    g.map_dataframe(custom_boxplot2, x='OS', y='value', color="darkgrey")
    g.axes[4].set(xlabel="Peptide Properties")
    g.axes[0].set(ylabel="counts / value")
    g.axes[3].set(ylabel="counts / value")
    plt.tight_layout()
