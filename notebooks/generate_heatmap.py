from __future__ import annotations

import os

import click
import gensim
import numpy as np
import pandas as pd
import seaborn as sns

from matchms import calculate_scores
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosine
from matplotlib import pyplot as plt
from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model
from spec2vec import Spec2Vec


def get_similarity_measure(measure_type: str, model_path: str | None = None) -> object:
    """Return the appropriate similarity measure based on the user's choice."""
    if measure_type == "modified_cosine":
        return ModifiedCosine()
    elif measure_type == "ms2deepscore":
        if not model_path:
            raise ValueError("Model path is required for MS2DeepScore")
        model = load_model(model_path)
        return MS2DeepScore(model)
    elif measure_type == "spec2vec":
        if not model_path:
            raise ValueError("Model path is required for Spec2Vec")
        model = gensim.models.Word2Vec.load(model_path)
        return Spec2Vec(
            model=model,
            intensity_weighting_power=0.5,
            allowed_missing_percentage=30.0,
        )
    else:
        raise ValueError("Unsupported similarity measure")


@click.command()
@click.option(
    "--path-root",
    default="miadbviz/data",
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
    "--measure-type",
    type=click.Choice(["modified_cosine", "ms2deepscore", "spec2vec"]),
    default="modified_cosine",
    show_default=True,
    help="Similarity measure to use for calculating similarity scores.",
)
@click.option(
    "--model-path",
    type=click.Path(exists=True),
    default=None,
    help="Path to the pre-trained model (required for MS2DeepScore and Spec2Vec).",
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
    measure_type: str,
    model_path: str,
    tolerance: float,
    min_score_threshold: float,
    max_spectra: int,
    output_size: tuple[int, int],
    cmap: str,
):
    """Calculate and visualize a spectral similarity heatmap using various similarity measures."""
    # Define file path
    file_path = os.path.join(path_root, filename)

    # Load spectra from the MGF file
    spectra: list = list(load_from_mgf(file_path))

    # Get the appropriate similarity measure
    similarity_measure = get_similarity_measure(measure_type, model_path)

    # Calculate similarity scores
    scores = calculate_scores(spectra, spectra, similarity_measure, is_symmetric=True)
    scores_array = scores.scores.to_array()

    # Limit the number of spectra to max_spectra
    m = min(max_spectra, len(spectra))
    spectra = spectra[:m]
    scores_array = scores_array[:m, :m]

    # If the similarity measure is ModifiedCosine, extract the score field
    if measure_type == "modified_cosine":
        scores_array = scores_array["ModifiedCosine_score"]

    # Apply minimum score threshold filter
    filtered_scores_array = np.where(
        scores_array >= min_score_threshold,
        scores_array,
        0,
    )

    # Extract metadata for heatmap labels
    titles = [s.metadata.get("title", f"Spectrum {i}") for i, s in enumerate(spectra)]
    skeletons = [
        s.metadata.get("skeleton", f"Skeleton {i}") for i, s in enumerate(spectra)
    ]

    # Create DataFrame for the heatmap visualization
    df = pd.DataFrame(filtered_scores_array, columns=titles, index=skeletons)

    # Plot the heatmap
    sns.set(rc={"figure.figsize": output_size})
    sns.heatmap(df, cmap=cmap, linewidths=0.5, linecolor="white")
    plt.title(f"MIADB {measure_type.capitalize()} Heatmap")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    generate_heatmap()
