"""Module to perform QC on the xiRT performance."""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pyteomics import electrochem

sns.set(context="notebook", style="white", palette="deep", font_scale=1)


def compute_sequence_features(peptides_df):
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
    peptides_df["isoelectric_point"] = [electrochem.pI(x) for x in peptides_df["sequence"].values]
    peptides_df["gravy"] = [electrochem.gravy(x) for x in peptides_df["sequence"].values]
    return peptides_df


def qc_peptides(features_df_host, features_df_trap, out_dir):
    """
    Plot a boxplot for the different peptide features.

    Args:
        features_df_host: df, host peptides with features
        features_df_trap: df, trap peptides with features
        out_dir: str, path to store qc data

    Returns:
        None
    """
    # make single df
    features_df = pd.concat([features_df_host, features_df_trap])

    # melt for easier plotting
    features_melt = features_df.melt(id_vars=["Type", "sequence"],
                                     value_vars=['length', 'KR', 'aromatic', 'acids', 'aliphatic',
                                                 'HGP', 'isoelectric_point', 'gravy'])

    # some log-values
    features_melt["log_value"] = np.log2(features_melt["value"])

    # peptides
    g = sns.FacetGrid(features_melt, col='variable', sharey=False, col_wrap=4)
    g.map_dataframe(sns.boxplot, x='Type', y='value', color="darkgrey")
    g.map_dataframe(sns.boxplot, x='Type', y='value', color="darkgrey")
    g.axes[4].set(xlabel="Peptide Properties")
    g.axes[7].set(xlabel="Peptide Properties")
    g.axes[5].set(xlabel="Peptide Properties")
    g.axes[6].set(xlabel="Peptide Properties")
    g.axes[0].set(ylabel="counts / value")
    g.axes[4].set(ylabel="counts / value")
    plt.tight_layout()

    g.savefig(os.path.join(out_dir, "qc_plot.png"))
