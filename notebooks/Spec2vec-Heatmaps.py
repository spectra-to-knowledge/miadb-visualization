from __future__ import annotations

import os

import gensim
import pandas as pd
import seaborn as sns
from matchms import calculate_scores
from matchms.importing import load_from_mgf
from matplotlib import pyplot as plt
from spec2vec import Spec2Vec

# Define paths for the spectra file and the pretrained model (update to your file paths)
path_root = "path/to/your/data/spectra_file.mgf"
spectra = list(load_from_mgf(path_root))

# Load pretrained Spec2Vec model
path_model = "path/to/your/models"
filename_model = "trained_spec2vec_model.model"
filename = os.path.join(path_model, filename_model)
model = gensim.models.Word2Vec.load(filename)

# Check spectra metadata and print if any intensity is above threshold
for s in spectra:
    if any(s.intensities > 1):
        print(s.metadata)

# Define Spec2Vec similarity measure with intensity weighting and allowed missing data threshold
spec2vec_similarity = Spec2Vec(
    model=model, intensity_weighting_power=0.5, allowed_missing_percentage=30.0
)

# Calculate similarity scores between spectra (symmetric scoring)
scores = calculate_scores(spectra, spectra, spec2vec_similarity, is_symmetric=True)
scores_array = scores.scores.to_array()

# Extract metadata for heatmap labels
titles = [s.metadata["title"] for s in spectra]  # Retrieve titles for spectra
skeletons = [s.metadata["skeleton"] for s in spectra]  # Retrieve skeleton info for spectra

# Create DataFrame for the heatmap visualization
df = pd.DataFrame(scores_array, columns=skeletons, index=skeletons)

# Generate and display a heatmap
sns.set(rc={"figure.figsize": (100, 100)})
sns.heatmap(df, cmap="viridis")
plt.title("MIADB Spec2Vec Heatmap")  # Set your title here
plt.show()
