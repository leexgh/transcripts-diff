import pandas as pd
import os

def process_grch37():
    # Define the directories using relative paths
    input_dir = os.path.join('input', 'grch37')
    output_dir = os.path.join('output', 'grch37')

    # Read the comparison file
    df = pd.read_csv(os.path.join(output_dir, 'TEMP_compare_result_all_enst.txt'), sep='\t')

    # Function to check if all values in a row are the same (ignoring empty strings)
    def are_values_different(row, columns):
        # Get non-empty values
        values = [str(row[col]).strip() for col in columns if str(row[col]).strip()]
        if not values:  # If all values are empty
            return False
        # Check if all values are the same
        return not all(v == values[0] for v in values)

    # Columns to compare (excluding 'ensembl-enst', 'msk-Impact-enst' for now)
    columns_to_compare = ['oncokb', 'mskcc', "MSK-Impact_enst_id"]

    # Find rows where values differ
    different_rows = df[df.apply(lambda row: are_values_different(row, columns_to_compare), axis=1)]

    # Sort by Hugo Symbol
    non_empty_mask = different_rows[columns_to_compare].applymap(
        lambda val: bool(pd.notna(val) and str(val).strip())
    )
    different_rows['non_empty_count'] = non_empty_mask.sum(axis=1)

    # Sort by count of non-empty columns (descending), then by Hugo Symbol
    different_rows = different_rows.sort_values(by=['identical_protein_sequence', 'non_empty_count', 'Hugo Symbol', ], ascending=[True, False, True])
    different_rows.drop(columns=['non_empty_count', 'ensembl_protein_sequence', 'oncokb_protein_sequence', 'mskcc_protein_sequence', 'MSK-Impact_enst_id_protein_sequence', 'oncokb_vs_MSK-Impact', 'oncokb_vs_mskcc'], inplace=True)
    # Save results
    output_path = os.path.join(output_dir, 'enst_differences.txt')
    different_rows.to_csv(output_path, sep='\t', index=False)

    # Print summary
    print(f"\nGRCh37 Results:")
    print(f"Total number of genes analyzed: {len(df)}")
    print(f"Number of genes with different transcript values: {len(different_rows)}")
    print(f"\nResults have been saved to: {output_path}")

    # Print some examples
    print("\nExample differences:")
    for _, row in different_rows.head().iterrows():
        print(f"\nGene: {row['Hugo Symbol']}")
        for col in columns_to_compare:
            print(f"{col}: {row[col]}")

    # Print instructions for including msk-Impact-enst
    print("\nTo include msk-Impact-enst in the comparison:")
    print("1. Add 'msk-Impact-enst' to the columns_to_compare list")
    print("2. The script will automatically include it in the comparison")

def process_grch38():
    # Define the directories using relative paths
    input_dir = os.path.join('input', 'grch38')
    output_dir = os.path.join('output', 'grch38')

    # Read the comparison file
    df = pd.read_csv(os.path.join(output_dir, 'compare_result_grch38.txt'), sep='\t')

    # Function to check if all values in a row are the same (ignoring empty strings)
    def are_values_different(row, columns):
        # Get non-empty values
        values = [str(row[col]).strip() for col in columns if str(row[col]).strip()]
        if not values:  # If all values are empty
            return False
        # Check if all values are the same
        return not all(v == values[0] for v in values)

    # Columns to compare
    columns_to_compare = ['oncokb', 'mskcc']

    # Find rows where values differ
    different_rows = df[df.apply(lambda row: are_values_different(row, columns_to_compare), axis=1)]

    # Sort by Hugo Symbol
    non_empty_mask = different_rows[columns_to_compare].applymap(
        lambda val: bool(pd.notna(val) and str(val).strip())
    )
    different_rows['non_empty_count'] = non_empty_mask.sum(axis=1)

    # Sort by count of non-empty columns (descending), then by Hugo Symbol
    different_rows = different_rows.sort_values(by=['non_empty_count', 'Hugo Symbol'], ascending=[False, True])
    different_rows.drop(columns=['non_empty_count'], inplace=True)

    # Save results
    output_path = os.path.join(output_dir, 'enst_differences.txt')
    different_rows.to_csv(output_path, sep='\t', index=False)

    # Print summary
    print(f"\nGRCh38 Results:")
    print(f"Total number of genes analyzed: {len(df)}")
    print(f"Number of genes with different transcript values: {len(different_rows)}")
    print(f"\nResults have been saved to: {output_path}")

    # Print some examples
    print("\nExample differences:")
    for _, row in different_rows.head().iterrows():
        print(f"\nGene: {row['Hugo Symbol']}")
        for col in columns_to_compare:
            print(f"{col}: {row[col]}")

if __name__ == "__main__":
    process_grch37()
    # process_grch38() 