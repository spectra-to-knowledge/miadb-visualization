# Libraries import
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from rdkit import DataStructs
from rdkit.Chem import AllChem, MolFromSmiles
from scipy.cluster.hierarchy import dendrogram, linkage


# Tanimoto scores calculation function
def tanimoto_calc(inc1, inc2):
    mol1 = MolFromSmiles(inc1)  # Conversion in inc1 and 2 (smiles) into molecular objects
    mol2 = MolFromSmiles(inc2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(
        mol1, 2, nBits=2048
    )  # Fingerprint generation using the Morgan algorithm (size 2)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
    s = round(
        DataStructs.TanimotoSimilarity(fp1, fp2), 3
    )  # Tanimoto scoring between the two fingerprints
    return s


df = pd.read_excel("metadata_file.xlsx")  # Replace with your excel metadata file's path

print(df.columns)  # We used 'Skeletons' and 'SMILES' columns for that particular example

sm = df["SMILES"].values.tolist()  # Creation of a SMILES list named sm

scores = np.zeros(
    (len(sm), len(sm))
)  # Creation of a square matrix (n x n); n being your number of skeletons (or smiles)
for i in range(len(sm)):
    for j in range(len(sm)):
        t = tanimoto_calc(sm[i], sm[j])  # Tanimoto scoring between SMILES i and SMILES j
        scores[i][j] = t
        scores[j][i] = t
l = linkage(scores, method="single", metric="euclidean")

sq = df[
    "Skeletons"
].values.tolist()  # Creation of the sq list containing all the skeletons contained in the excel file

s = []
for i in range(len(scores)):
    s.append(
        [sq[i]] + list(scores[i])
    )  # Every skeleton is being paired into a list with its corresponding Tanimoto scores (with all other skeletons + itself)

DF = pd.DataFrame(
    s, columns=["skeletons"] + sq
)  # Creation of a DataFrame containing the skeletons name and their corresponding Tanimoto score one to each other

print(DF)

labelList = sm = df["Skeletons"].values.tolist()  # Creation of the labels
plt.figure(figsize=(10, 15))
dendrogram(
    l,
    orientation="top",
    labels=sm,
    distance_sort="descending",
    show_leaf_counts=False,
    leaf_font_size=9,
)  # Creation of the dendogram
plt.show()
