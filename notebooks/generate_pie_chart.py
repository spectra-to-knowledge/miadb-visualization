#!/usr/bin/env python

import click
import matplotlib.pyplot as plt
import pandas as pd


@click.command()
@click.option(
    "--file-path",
    default="miadbviz/data/ajmalicine-spirooxindoles-best-queries-repartition.xlsx",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the Excel file containing ion intensity data.",
)
@click.option(
    "--sheet-name",
    default="ajmalicine-spiro-best-queries-r",
    type=str,
    required=True,
    help="Name of the sheet in the Excel file to be used.",
)
@click.option(
    "--output-figure",
    type=click.Path(writable=True),
    default=None,
    help="Path to save the pie chart figure (optional).",
)
def generate_pie_chart(file_path: str, sheet_name: str, output_figure: str):
    """Generate a pie chart of ion distribution by plant genus based on intensity."""

    # Read the data from the Excel file
    df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)

    # Extract plant names from the first row
    plant_names = df.iloc[0].tolist()

    # Extract data from the second row onward
    data = df.iloc[1:].reset_index(drop=True)
    data.columns = plant_names  # Assign plant names as column headers

    # Group duplicate plant names by averaging their intensities
    data = data.groupby(by=data.columns, axis=1).mean()

    # Print the columns that have been averaged (duplicates) and the total column count after averaging
    print("Modified columns (duplicates averaged):")
    print([col for col in plant_names if plant_names.count(col) > 1])
    print("Total number of columns after averaging duplicates:", len(data.columns))

    # Transpose the data so plants become rows
    data_transposed = data.T
    data_transposed.columns = [f"Ion_{i + 1}" for i in range(len(data_transposed.columns))]  # Rename columns

    # Reset index so plant names become a column
    data_transposed = data_transposed.reset_index()
    data_transposed.columns = ["Plant"] + list(data_transposed.columns[1:])

    # Convert the data to long format for each plant
    data_long = pd.melt(data_transposed, id_vars=["Plant"], var_name="Ion", value_name="Intensity")

    # Extract genus from plant names
    data_long["Genus"] = data_long["Plant"].apply(lambda x: x.split()[0] if isinstance(x, str) else "Unknown")

    # Remove rows with zero intensity
    data_long = data_long[data_long["Intensity"] > 0]

    # Group by genus and calculate total intensity and count of ions
    grouped_data = (
        data_long.groupby("Genus")
        .agg(total_intensity=("Intensity", "sum"), total_count=("Intensity", "count"))
        .reset_index()
    )

    # Exclude genera with zero total intensity
    grouped_data = grouped_data[grouped_data["total_intensity"] > 0]

    # Combine genera with less than 1% of total intensity into an "Other" category
    total_intensity_sum = grouped_data["total_intensity"].sum()
    grouped_data["Percentage"] = (grouped_data["total_intensity"] / total_intensity_sum) * 100
    small_genres = grouped_data[grouped_data["Percentage"] < 1]

    # Create an "Other" category if necessary
    if not small_genres.empty:
        other_intensity = small_genres["total_intensity"].sum()
        other_count = small_genres["total_count"].sum()
        other_genres = ", ".join(small_genres["Genus"])

        grouped_data = grouped_data[grouped_data["Percentage"] >= 1]
        other_data = pd.DataFrame({
            "Genus": ["Other"],
            "total_intensity": [other_intensity],
            "total_count": [other_count],
            "Percentage": [small_genres["Percentage"].sum()],
        })
        grouped_data = pd.concat([grouped_data, other_data], ignore_index=True)

    # Create a pie chart based on the total intensity for each genus
    plt.figure(figsize=(12, 8))
    plt.pie(
        grouped_data["total_intensity"],
        labels=grouped_data["Genus"],
        autopct="%1.1f%%",
        startangle=140,
        pctdistance=0.85,  # Distance of percentage labels from center
        labeldistance=1.05,  # Distance of labels from center
        textprops={"fontsize": 10},  # Label font size
        wedgeprops={"edgecolor": "black"},  # Adds borders to the slices for clarity
    )
    plt.title("Ion distribution by plant genus (based on intensity)")

    # Add a legend detailing the genera included in "Other"
    if "Other" in grouped_data["Genus"].values:
        plt.figtext(
            0.5,
            -0.1,
            f"Genera included in 'Other': {other_genres}",
            ha="center",
            fontsize=10,
            bbox=dict(facecolor="white", alpha=0.5),
        )

    plt.tight_layout()

    # Save or show the plot
    if output_figure:
        plt.savefig(output_figure)
        print(f"Pie chart saved to {output_figure}")
    else:
        plt.show()


if __name__ == "__main__":
    generate_pie_chart()
