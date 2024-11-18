from __future__ import annotations

import os

import pandas as pd
import seaborn as sns
from matchms import calculate_scores
from matchms.importing import load_from_mgf
from matplotlib import pyplot as plt
from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model

# Define paths for the spectra file and pretrained model (to be specified by the user)
path_root = "path/to/your/data"
file = os.path.join(path_root, "spectra_file.mgf")
path_root_model = "path/to/your/models"
filename_model = "trained_model.hdf5"  # You can access and download the MS2DeepScore model via https://github.com/matchms/ms2deepscore (created by Florian Huber)

# Load spectra from the .mgf file
spectra = list(load_from_mgf(file))

# Load the pretrained MS2DeepScore model
model = load_model(os.path.join(path_root_model, filename_model))
similarity_measure = MS2DeepScore(model)

# Calculate similarity scores between spectra (symmetric scoring)
scores = calculate_scores(spectra, spectra, similarity_measure, is_symmetric=True)
scores_array = scores.scores.to_array()

# Create a DataFrame for the heatmap visualization
titles = [s.metadata["title"] for s in spectra]  # Retrieve titles for spectra
df = pd.DataFrame(scores_array, columns=titles, index=titles)

# Generate and display a heatmap
sns.set(rc={"figure.figsize": (150, 150)})
sns.heatmap(df, cmap="viridis")
plt.title("MIADB MS2DeepScore Heatmap")
plt.show()
