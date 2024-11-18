===================
MIADB-visualization
===================

    This repository provides tools and scripts designed for chemists and researchers who are new to Python scripting, specifically in the context of mass spectrometry (MS/MS) and molecular structure analysis.

The repository leverages the **MatchMS** Python library to calculate various MS/MS similarity scores, allowing for a more effective comparison of mass spectrometry data.
It also includes functionalities for visualizing these similarity scores using **Seaborn** heatmaps. 

Additionally, the repository contains a script for calculating **Tanimoto similarity scores** between multiple molecular structures and visualizing structural similarities as a tree.
Another script is included to visualize the ion distribution within complex datasets, enhancing the understanding of mass spectrometry results.

This repository was used during the **MIADB update** to explore the chemical diversity of **Monoterpene Indole Alkaloids (MIAs)**.
The tools and methods in this repository are applicable to a wide range of chemical and biological research involving mass spectrometry and molecular similarity analysis.

Features
--------

- **MS/MS Similarity Calculation**: 
  - Use of the **MatchMS** library to calculate various MS/MS similarity scores. 
  - Reference: [Huber et al., (2020). *matchms - processing and similarity evaluation of mass spectrometry data*. Journal of Open Source Software, 5(52), 2411.](https://doi.org/10.21105/joss.02411)
  - [MatchMS documentation](https://matchms.readthedocs.io/en/latest/)

- **Heatmap Visualization**: 
  - Visualize MS/MS similarity scores as heatmaps using the **Seaborn** library.

- **Tanimoto Structural Similarity**: 
  - Compute **Tanimoto similarity scores** between molecular structures.
  - Visualize the structural similarities as a **tree diagram**.

- **Ion Distribution Visualization**:
  - A script to visualize the repartition of ions in complex mass spectrometry datasets.


💪 Getting Started
------------------

⚠️ WORK IN PROGRESS

.. code-block:: sh

    if command -v poetry &> /dev/null; then echo "Poetry is already installed."; poetry --version; else echo "Poetry is not installed. Installing Poetry..."; curl -sSL https://install.python-poetry.org | python3 -; fi


.. code-block:: sh

    poetry install
    poetry run python3 notebooks/TODO.py

⚠️ WORK IN PROGRESS
