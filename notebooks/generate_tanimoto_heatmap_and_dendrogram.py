#!/usr/bin/env python

from __future__ import annotations

import click
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage


def tanimoto_calc(inc1: str, inc2: str) -> float:
    """Calculate the Tanimoto score between two SMILES strings."""
    mol1 = MolFromSmiles(inc1)
    mol2 = MolFromSmiles(inc2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    score = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)
    return score


@click.command()
@click.option(
    "--metadata-file",
    default="miadbviz/data/skeletons-smiles.xlsx",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the metadata Excel file containing SMILES and Skeletons columns.",
)
@click.option(
    "--smiles-column",
    type=str,
    default="SMILES",
    show_default=True,
    help="Name of the column in the Excel file containing SMILES strings.",
)
@click.option(
    "--skeleton-column",
    type=str,
    default="Skeletons",
    show_default=True,
    help="Name of the column in the Excel file containing skeleton names.",
)
@click.option(
    "--method",
    type=click.Choice(["single", "complete", "average", "ward"]),
    default="single",
    show_default=True,
    help="Clustering method to use in the linkage step.",
)
@click.option(
    "--dendrogram-figsize",
    type=(int, int),
    default=(10, 15),
    show_default=True,
    help="Figure size for the dendrogram plot (width, height).",
)
@click.option(
    "--output-dendrogram",
    type=click.Path(writable=True),
    default=None,
    help="Path to save the dendrogram plot image (optional).",
)
def generate_tanimoto_heatmap_and_dendrogram(
    metadata_file: str,
    smiles_column: str,
    skeleton_column: str,
    method: str,
    dendrogram_figsize: tuple[int, int],
    output_dendrogram: str,
):
    """Generate a Tanimoto similarity heatmap and a dendrogram from SMILES in an Excel file."""
    # Read the metadata Excel file
    df = pd.read_excel(metadata_file)
    print(df.columns)  # Print column names for debugging

    sm = df[smiles_column].values.tolist()  # SMILES list
    skeletons = df[skeleton_column].values.tolist()  # Skeleton list

    # Create the similarity matrix
    scores = np.zeros((len(sm), len(sm)))
    for i in range(len(sm)):
        for j in range(len(sm)):
            score = tanimoto_calc(sm[i], sm[j])  # Calculate Tanimoto similarity
            scores[i][j] = score
            scores[j][i] = score  # Symmetric matrix

    # Perform hierarchical clustering
    linkage_matrix = linkage(scores, method=method, metric="euclidean")

    # Prepare data for heatmap and dendrogram
    s = [[skeletons[i]] + list(scores[i]) for i in range(len(scores))]
    df_scores = pd.DataFrame(s, columns=[skeleton_column] + skeletons)

    print(df_scores)

    # Plot the dendrogram
    plt.figure(figsize=dendrogram_figsize)
    dendrogram(
        linkage_matrix,
        orientation="top",
        labels=skeletons,
        distance_sort="descending",
        show_leaf_counts=False,
        leaf_font_size=9,
    )

    # Save or show dendrogram
    if output_dendrogram:
        plt.savefig(output_dendrogram)
        print(f"Dendrogram saved to {output_dendrogram}")
    else:
        plt.show()


if __name__ == "__main__":
    generate_tanimoto_heatmap_and_dendrogram()
