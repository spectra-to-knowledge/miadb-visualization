import pandas as pd
import matplotlib.pyplot as plt

# Load data from an Excel file
file_path = "your_file_path.xlsx"  # Placeholder for the file path (you can use the .tsv output from your molecular network and delete manually the unuseful columns to only retain the biological sources and their corresponding ion intensities)
sheet_name = "your_sheet_name"      # Placeholder for the sheet name

# Load data from the Excel file
df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)

# Extract plant names from the first row
plant_names = df.iloc[0].tolist()

# Data starts from the second row
data = df.iloc[1:].reset_index(drop=True)
data.columns = plant_names  # Set column names using plant names

# Transpose data to have plants as columns
data_transposed = data.T
data_transposed.columns = [f'Ion_{i + 1}' for i in range(len(data_transposed.columns))]  # Rename columns

# Reset index to make plant names a column
data_transposed = data_transposed.reset_index()
data_transposed.columns = ['Plant'] + list(data_transposed.columns[1:])

# Convert data to long format for each plant
data_long = pd.melt(data_transposed, id_vars=['Plant'], var_name='Ion', value_name='Intensity')

# Extract genus from plant names
data_long['Genus'] = data_long['Plant'].apply(lambda x: x.split()[0] if isinstance(x, str) else 'Unknown')

# Remove rows with zero intensity
data_long = data_long[data_long['Intensity'] > 0]

# Group data by genus and calculate necessary statistics
grouped_data = data_long.groupby('Genus').agg(
    total_intensity=('Intensity', 'sum'),
    total_count=('Intensity', 'count')
).reset_index()

# Exclude genera with zero total intensity
grouped_data = grouped_data[grouped_data['total_intensity'] > 0]

# Group genera with less than 1% intensity into an "Others" category
total_intensity_sum = grouped_data['total_intensity'].sum()
grouped_data['Percentage'] = (grouped_data['total_intensity'] / total_intensity_sum) * 100
small_genres = grouped_data[grouped_data['Percentage'] < 1]

# Create an "Others" category if needed
if not small_genres.empty:
    other_intensity = small_genres['total_intensity'].sum()
    other_count = small_genres['total_count'].sum()
    other_genres = ', '.join(small_genres['Genus'])

    grouped_data = grouped_data[grouped_data['Percentage'] >= 1]
    other_data = pd.DataFrame({
        'Genus': ['Others'],
        'total_intensity': [other_intensity],
        'total_count': [other_count],
        'Percentage': [small_genres['Percentage'].sum()]
    })
    grouped_data = pd.concat([grouped_data, other_data], ignore_index=True)

# Plot a pie chart based on the sum of intensities for each genus
plt.figure(figsize=(12, 8))
plt.pie(
    grouped_data['total_intensity'],
    labels=grouped_data['Genus'],
    autopct='%1.1f%%',
    startangle=140,
    pctdistance=0.85,
    labeldistance=1.05,
    textprops={'fontsize': 10},
    wedgeprops={'edgecolor': 'black'}
)
plt.title('Distribution of ions by plant genus (based on intensity)')

# Add a legend for genera grouped under "Others"
if 'Others' in grouped_data['Genus'].values:
    plt.figtext(
        0.5, -0.1,
        f"Genera included in 'Others': {other_genres}",
        ha='center',
        fontsize=10,
        bbox=dict(facecolor='white', alpha=0.5)
    )

plt.tight_layout()
plt.show()
