import os

import pandas as pd
from rdkit import Chem


def convert_smiles_to_smarts(input_file, output_file, sheet_name="Sheet1"):
    """
    Converts chemical structures from SMILES format to SMARTS format and saves them in an Excel file.

    Parameters:
    - input_file: str, path to the input Excel file containing a 'SMILES' column.
    - output_file: str, path to the output Excel file where the SMARTS will be added.
    - sheet_name: str, name of the Excel sheet to read from (default is 'Sheet1').
    """

    # Read the Excel file containing SMILES
    excel_data = pd.read_excel(input_file, sheet_name=sheet_name)

    # Get the list of SMILES from the Excel file
    smiles_list = list(excel_data["SMILES"])

    smarts_list = []
    # Convert SMILES to SMARTS
    for smiles in smiles_list:
        print(smiles)
        mol = Chem.MolFromSmiles(smiles)
        smarts = (
            Chem.MolToSmarts(mol, isomericSmiles=True) if mol else None
        )  # Check if mol is valid
        smarts_list.append(smarts)

    # Add the SMARTS column to the existing DataFrame
    excel_data["SMARTS"] = smarts_list

    # Save the updated Excel file with the new 'SMARTS' column
    excel_data.to_excel(output_file, index=False)

    # Print information for verification
    print(f"Number of SMARTS generated: {len(smarts_list)}")
    print(excel_data.head())  # Display the first few rows for verification


if __name__ == "__main__":
    # Define file paths (you can modify these paths or use relative paths)
    input_file_path = input("Enter the path to the input Excel file: ")
    output_file_path = input("Enter the path to save the output Excel file: ")

    # Run the conversion
    convert_smiles_to_smarts(input_file=input_file_path, output_file=output_file_path)
