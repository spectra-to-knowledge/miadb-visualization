#!/usr/bin/env python

from __future__ import annotations

import os
from typing import List

import click
import numpy as np
import pandas as pd
import seaborn as sns
from matchms import calculate_scores
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosine
from matplotlib import pyplot as plt


@click.command()
@click.option(
    "--path-root",
    default="src/miadbviz/data",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Directory path where the MGF file is located.",
)
@click.option(
    "--filename",
    default="MIADB-monomers.mgf",
    type=str,
    required=True,
    help="Name of the MGF file to process.",
)
@click.option(
    "--tolerance",
    type=float,
    default=0.005,
    show_default=True,
    help="Tolerance for the ModifiedCosine similarity measure.",
)
@click.option(
    "--min-score-threshold",
    type=float,
    default=0.0,
    show_default=True,
    help="Minimum score threshold to filter similarity scores.",
)
@click.option(
    "--max-spectra",
    type=int,
    default=321,
    show_default=True,
    help="Maximum number of spectra to include in the heatmap.",
)
@click.option(
    "--output-size",
    type=(int, int),
    default=(50, 50),
    show_default=True,
    help="Figure size for the heatmap (width, height).",
)
@click.option(
    "--cmap",
    type=str,
    default="viridis",
    show_default=True,
    help="Colormap for the heatmap.",
)
def generate_heatmap(
    path_root: str,
    filename: str,
    tolerance: float,
    min_score_threshold: float,
    max_spectra: int,
    output_size: tuple[int, int],
    cmap: str,
):
    """Calculate and visualize a spectral similarity heatmap using Modified Cosine similarity."""
    # Define file path
    file_path = os.path.join(path_root, filename)

    # Load spectra from the MGF file
    spectra: List = list(load_from_mgf(file_path))

    # Define the similarity measure
    similarity_measure = ModifiedCosine(tolerance=tolerance)

    # Calculate similarity scores
    scores = calculate_scores(spectra, spectra, similarity_measure, is_symmetric=True)
    scores_array = scores.scores.to_array()

    # Limit the number of spectra to max_spectra
    m = min(max_spectra, len(spectra))
    spectra = spectra[:m]
    scores_array = scores_array[:m, :m]

    # Convert to a regular NumPy array of floats
    scores_array = scores_array["ModifiedCosine_score"]
    
    # Apply minimum score threshold filter
    filtered_scores_array = np.where(scores_array >= min_score_threshold, scores_array, 0)

    # Extract metadata for heatmap labels
    titles = [s.metadata.get("title", f"Spectrum {i}") for i, s in enumerate(spectra)]
    skeletons = [s.metadata.get("skeleton", f"Skeleton {i}") for i, s in enumerate(spectra)]

    # Create DataFrame for the heatmap
    df = pd.DataFrame(filtered_scores_array, columns=titles, index=skeletons)

    # Plot the heatmap
    sns.set(rc={"figure.figsize": output_size})
    sns.heatmap(df, cmap=cmap, linewidths=0.5, linecolor="white")
    plt.title("MIADB Modified Cosine Heatmap")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    generate_heatmap()
