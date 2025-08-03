===================
MIADB-visualization
===================

    This repository provides tools and scripts designed for chemists and researchers who are new to Python scripting, specifically in the context of mass spectrometry (MS/MS) and molecular structure analysis.

The repository leverages `matchms <https://matchms.readthedocs.io/en/latest/>`_ to calculate various MS/MS similarity scores, allowing for a more effective comparison of mass spectrometry data.
It also includes functionalities for visualizing these similarity scores using `seaborn <https://github.com/mwaskom/seaborn>`_ heatmaps.

Additionally, the repository contains a script for calculating **Tanimoto coefficients** between multiple molecular structures and visualizing structural similarities as a tree.
Another script is included to visualize the ion distribution within complex datasets, enhancing the understanding of mass spectrometry results.

This repository was used during the **MIADB update** to explore the chemical diversity of **Monoterpene Indole Alkaloids (MIAs)**.
The tools and methods in this repository are applicable to a wide range of chemical and biological research involving mass spectrometry and molecular similarity analysis.

Features
--------

* **MS/MS Similarity Calculation**:
    * Use of the **matchms** library to calculate various MS/MS similarity scores.
        * `Reference <https://doi.org/10.21105/joss.02411>`_
        * `Documentation <https://matchms.readthedocs.io/en/latest/>`_

* **Heatmap Visualization**:
    * Visualize MS/MS similarity scores as heatmaps using the **Seaborn** library.

* **Tanimoto Structural Similarity**:
    * Compute **Tanimoto similarity scores** between molecular structures.
    * Visualize the structural similarities as a **tree diagram**.

* **Ion Distribution Visualization**:
    * A script to visualize the repartition of ions in complex mass spectrometry datasets.


💪 Getting Started
------------------

We use `uv <https://docs.astral.sh/uv/>`_ to manage dependencies.

.. code-block:: sh

    pip install uv


.. code-block:: sh

    uv sync


Then, you should be able to run the different notebooks:

.. code-block:: sh

    uv run python3 notebooks/generate_heatmap.py
    # uv run python3 notebooks/generate_heatmap.py --measure-type spec2vec --model-path ../../Downloads/spec2vec_AllPositive_ratio05_filtered_201101_iter_15.model
    # uv run python3 notebooks/generate_heatmap.py --measure-type ms2deepscore --model-path ../../Downloads/ms2deepscore_model.pt


⚠️ For the Spec2Vec and MS2DeepScore metrics to work, you'll need to download the models first as stated in their respective documentations.

.. code-block:: sh

    uv run python3 notebooks/generate_tanimoto_heatmap_and_dendrogram.py


.. code-block:: sh

    uv run python3 notebooks/generate_pie_chart.py
    uv run python3 notebooks/generate_pie_chart.py --file-path src/miadbviz/data/ajmalicine-spirooxindoles-corynantheane-spirooxindoles-combined-repartition.xlsx --sheet-name ajmalicine-spiro-corynanthean-s
    uv run python3 notebooks/generate_pie_chart.py --file-path src/miadbviz/data/corynanthean-spirooxindoles-best-queries-repartition.xlsx --sheet-name corynanthean-spiro-best-queries


In case you need help:

.. code-block:: sh

    uv run python3 notebooks/generate_heatmap.py --help
