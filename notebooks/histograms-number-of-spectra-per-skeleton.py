import matplotlib.pyplot as plt
from matchms.importing import load_from_mgf
from collections import Counter
import numpy as np

# Define the path to the MGF file containing annotated spectra
# Replace "path_to_file.mgf" with the actual file path
mgf_path = "path_to_file.mgf"

# Load spectra from the MGF file
spectra = load_from_mgf(mgf_path)

# Extract the "skeleton" metadata from each spectrum
# If the metadata key "skeleton" is missing, it will return None by default
skeletons = [s.metadata.get("skeleton") for s in spectra]

# Count the number of spectra associated with each unique skeleton
skeleton_counts = Counter(skeletons)

# Prepare data for the histogram
# skeleton_names contains unique skeleton identifiers
# skeleton_values contains the corresponding counts of spectra for each skeleton
skeleton_names = list(skeleton_counts.keys())
skeleton_values = list(skeleton_counts.values())

# Generate a color gradient for the bars in the histogram
colors = plt.cm.viridis(np.linspace(0, 1, len(skeleton_names)))

# Create the vertical bar chart
plt.figure(figsize=(10, 6))
bars = plt.bar(skeleton_names, skeleton_values, color=colors, width=0.7)

# Adjust the y-axis limits to ensure all annotations fit within the plot
plt.ylim(0, max(skeleton_values) * 1.2)

# Annotate the bars with their respective counts, oriented vertically
for bar, value in zip(bars, skeleton_values):
    plt.text(
        bar.get_x() + bar.get_width() / 2,  # Horizontal position
        bar.get_height() + 0.5,  # Vertical position slightly above the bar
        str(value),  # Text to display
        ha="center",  # Horizontal alignment
        va="bottom",  # Vertical alignment
        rotation=90,  # Rotate text vertically
        fontsize=10,  # Font size for the annotation
    )

# Configure axis labels and tick rotation for clarity
plt.xticks(rotation=90, ha="center")  # Rotate x-axis labels vertically
plt.xlabel("Skeleton")  # Label for the x-axis
plt.ylabel("Number of spectra")  # Label for the y-axis
plt.title("Number of spectra per skeleton")  # Title of the plot
plt.tight_layout()  # Adjust layout for better appearance

# Display the histogram
plt.show()
