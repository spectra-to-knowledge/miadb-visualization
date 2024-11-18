# Libraries and modules import
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import seaborn as sns
from matchms import calculate_scores
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosine
from matplotlib import pyplot as plt

# Setting your path
path_root = "directory_path"  # input the directory path
file = os.path.join(path_root, "filename")  # input the filename

# Creating a spectra list
spectra = list(load_from_mgf(file))

# Calculating the similarity score
similarity_measure = ModifiedCosine(
    tolerance=0.005
)  # you're free to select another similarity score and to set the tolerance (don't forget to import the module before)
scores = calculate_scores(spectra, spectra, similarity_measure, is_symmetric=True)
scores_array = scores.scores.to_array()

m = 321  # Modify the value with the number of spectra that you want to plot

scores_array = np.zeros((m, m), dtype=float)  # Creating an array of size m by m
for sp_i in range(len(spectra)):
    for sp_j in range(sp_i, len(spectra)):
        sc = float(similarity_measure.pair(spectra[sp_i], spectra[sp_j])["score"])
        scores_array[sp_i][sp_j] = sc
        scores_array[sp_j][sp_i] = sc

""" If you want to apply a minimum score filter : 
spectra=spectra[:m]
n= len(spectra)
scores_array=np.zeros((m,m),dtype=float)
for sp_i in range(len(spectra)):
    for sp_j in range(sp_i,len(spectra)):
        sc=float(similarity_measure.pair(spectra[sp_i],spectra[sp_j])['score'])
        if sc <= 0.9:
            sc = 0
        scores_array[sp_i][sp_j] = sc
        scores_array[sp_j][sp_i] = sc """

# Axes definition
# You're free to input the metadata you want to visualize in your heatmap
l = [s.metadata["title"] for s in spectra]
sk = [s.metadata["skeleton"] for s in spectra]
df = pd.DataFrame(scores_array, columns=l, index=sk)

# Heatmap plotting
sns.set(rc={"figure.figsize": (100, 100)})  # You can modulate the output size
sns.heatmap(df, cmap="viridis")  # You can change the viridis value to set another color palette
plt.title("MIADB Modified Cosine Heatmap")  # Set your title here
plt.show()
