import pandas as pd
import os

def process_grch37():
    # Define the directories using relative paths
    input_dir = os.path.join('input', 'grch37')
    output_dir = os.path.join('output', 'grch37')

    # Read the input files
    def read_file(file_path, index_col, sep='\t'):
        try:
            df = pd.read_csv(file_path, sep=sep, index_col=index_col, low_memory=False)
            df.index = df.index.astype(str)
            return df
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            return None

    # Read all files with full paths
    ensembl_df = read_file(os.path.join(input_dir, 'ensembl_biomart_canonical_transcripts_per_hgnc.txt'), 'hgnc_symbol')
    oncokb_df = read_file(os.path.join(input_dir, 'oncokb_curated_genes.txt'), 'Hugo Symbol')
    mskcc_df = read_file(os.path.join(input_dir, 'isoform_overrides_at_mskcc.txt'), 'gene_name')
    msk_impact_df = read_file(os.path.join(input_dir, 'Iv7_dmp_isoform_merged_overrides.txt'), 'gene_name')

    # Get unique gene symbols from the three sources
    oncokb_genes = set(oncokb_df.index)
    mskcc_genes = set(mskcc_df.index)
    msk_impact_genes = set(msk_impact_df.index)

    # Combine all unique genes
    all_genes = sorted(list(oncokb_genes.union(mskcc_genes).union(msk_impact_genes)), key=str)

    # Create a new DataFrame with the desired structure
    result_df = pd.DataFrame(index=all_genes, columns=[
        'Hugo Symbol',
        'ensembl',
        'oncokb',
        'mskcc',
        'MSK-Impact'
    ])

    # Fill in the Hugo Symbol column
    result_df['Hugo Symbol'] = result_df.index

    # Handle duplicates by grouping and joining values with commas
    def get_duplicate_handled_dict(df, col_name):
        if col_name not in df.columns:
            print(f"Warning: Column {col_name} not found in DataFrame")
            return {}
        # Group by index and join values with comma
        return df.groupby(df.index)[col_name].agg(lambda x: ','.join(str(i) for i in x if pd.notna(i))).to_dict()

    # Create dictionaries with duplicate handling
    ensembl_map = get_duplicate_handled_dict(ensembl_df, 'ensembl_canonical_transcript')
    oncokb_map = get_duplicate_handled_dict(oncokb_df, 'GRCh37 Isoform')
    mskcc_map = get_duplicate_handled_dict(mskcc_df, '#enst_id')
    msk_impact_map = get_duplicate_handled_dict(msk_impact_df, '#isoform_override')

    # Fill in the transcript columns
    result_df['ensembl'] = result_df.index.map(lambda x: ensembl_map.get(x, ''))
    result_df['oncokb'] = result_df.index.map(lambda x: oncokb_map.get(x, ''))
    result_df['mskcc'] = result_df.index.map(lambda x: mskcc_map.get(x, ''))
    result_df['MSK-Impact'] = result_df.index.map(lambda x: msk_impact_map.get(x, ''))

    # Save the result in the output directory
    output_path = os.path.join(output_dir, 'compare_result_grch37.txt')
    # Save without index
    result_df.to_csv(output_path, sep='\t', index=False)

    print(f"\nGRCh37 Results:")
    print(f"Total number of unique genes: {len(all_genes)}")
    print(f"Results have been saved to {output_path}")

    # Print some statistics
    print("\nStatistics:")
    print(f"Number of genes from Oncokb: {len(oncokb_genes)}")
    print(f"Number of genes from MSKCC: {len(mskcc_genes)}")
    print(f"Number of genes from MSK-Impact: {len(msk_impact_genes)}")
    print(f"Number of genes from Ensembl: {len(ensembl_df.index)}")

    # Print example of duplicate handling
    print("\nExample of duplicate handling:")
    duplicate_example = result_df[result_df['mskcc'].str.contains(',')].head(1)
    if not duplicate_example.empty:
        print("\nExample of a gene with multiple transcripts:")
        print(duplicate_example.to_string())

def process_grch38():
    # Define the directories using relative paths
    input_dir = os.path.join('input', 'grch38')
    output_dir = os.path.join('output', 'grch38')

    # Read the input files
    def read_file(file_path, index_col, sep='\t'):
        try:
            df = pd.read_csv(file_path, sep=sep, index_col=index_col, low_memory=False)
            df.index = df.index.astype(str)
            return df
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            return None

    # Read all files with full paths
    ensembl_df = read_file(os.path.join(input_dir, 'ensembl_biomart_canonical_transcripts_per_hgnc.txt'), 'hgnc_symbol')
    oncokb_df = read_file(os.path.join(input_dir, 'oncokb_curated_genes.txt'), 'Hugo Symbol')
    mskcc_df = read_file(os.path.join(input_dir, 'isoform_overrides_at_mskcc_grch38.txt'), 'gene_name')
    mane_df = read_file(os.path.join(input_dir, 'isoform_overrides_mane_grch38.txt'), 'gene_name')

    # Get unique gene symbols from Oncokb and MSKCC only
    oncokb_genes = set(oncokb_df.index)
    mskcc_genes = set(mskcc_df.index)

    # Combine unique genes from Oncokb and MSKCC only
    all_genes = sorted(list(oncokb_genes.union(mskcc_genes)), key=str)

    # Create a new DataFrame with the desired structure
    result_df = pd.DataFrame(index=all_genes, columns=[
        'Hugo Symbol',
        'ensembl',
        'oncokb',
        'mskcc',
        'mane'
    ])

    # Fill in the Hugo Symbol column
    result_df['Hugo Symbol'] = result_df.index

    # Handle duplicates by grouping and joining values with commas
    def get_duplicate_handled_dict(df, col_name):
        if col_name not in df.columns:
            print(f"Warning: Column {col_name} not found in DataFrame")
            return {}
        # Group by index and join values with comma
        return df.groupby(df.index)[col_name].agg(lambda x: ','.join(str(i) for i in x if pd.notna(i))).to_dict()

    # Create dictionaries with duplicate handling
    ensembl_map = get_duplicate_handled_dict(ensembl_df, 'ensembl_canonical_transcript')
    oncokb_map = get_duplicate_handled_dict(oncokb_df, 'GRCh38 Isoform')
    mskcc_map = get_duplicate_handled_dict(mskcc_df, '#enst_id')
    mane_map = get_duplicate_handled_dict(mane_df, '#enst_id')

    # Fill in the transcript columns
    result_df['ensembl'] = result_df.index.map(lambda x: ensembl_map.get(x, ''))
    result_df['oncokb'] = result_df.index.map(lambda x: oncokb_map.get(x, ''))
    result_df['mskcc'] = result_df.index.map(lambda x: mskcc_map.get(x, ''))
    result_df['mane'] = result_df.index.map(lambda x: mane_map.get(x, ''))

    # Save the result in the output directory
    output_path = os.path.join(output_dir, 'compare_result_grch38.txt')
    # Save without index
    result_df.to_csv(output_path, sep='\t', index=False)

    print(f"\nGRCh38 Results:")
    print(f"Total number of unique genes: {len(all_genes)}")
    print(f"Results have been saved to {output_path}")

    # Print some statistics
    print("\nStatistics:")
    print(f"Number of genes from Oncokb: {len(oncokb_genes)}")
    print(f"Number of genes from MSKCC: {len(mskcc_genes)}")
    print(f"Number of genes from MANE: {len(mane_df.index)}")
    print(f"Number of genes from Ensembl: {len(ensembl_df.index)}")

    # Print example of duplicate handling
    print("\nExample of duplicate handling:")
    duplicate_example = result_df[result_df['mskcc'].str.contains(',')].head(1)
    if not duplicate_example.empty:
        print("\nExample of a gene with multiple transcripts:")
        print(duplicate_example.to_string())

if __name__ == "__main__":
    process_grch37()
    process_grch38() 