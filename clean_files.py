import os
import pandas as pd
from rdkit import Chem

# ----------------------- Configuration -----------------------
INPUT_DIR = 'ABC-Ts-v3/csv'  # Set this to the directory containing the csv files
MASTER_FILE = os.path.join(INPUT_DIR, 'abc_transporters_master.csv')
CLEANED_DIR = os.path.join(INPUT_DIR, 'cleaned')
os.makedirs(CLEANED_DIR, exist_ok=True)

# ----------------------- Helper Functions -----------------------
def is_valid_smiles(smiles):
    if not smiles or pd.isna(smiles):
        return False
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

# ----------------------- Main Script -----------------------
def main():
    # Load master data for mapping pref_name
    master_df = pd.read_csv(MASTER_FILE)
    # Create a lookup for pref_name by target_chembl_id
    target_name_map = dict(zip(master_df['target_chembl_id'], master_df['pref_name']))
    
    # We will accumulate stats in a dictionary keyed by target_chembl_id
    stats = {}
    # stats[target_chembl_id] = {
    #    'pref_name': str,
    #    'inhibitors_count': int,
    #    'substrates_count': int,
    #    'cleaned_inhibitors_count': int,
    #    'cleaned_substrates_count': int
    # }

    # Iterate through all CSV files (except the master file) in INPUT_DIR
    for filename in os.listdir(INPUT_DIR):
        if not filename.endswith('.csv'):
            continue
        if filename == 'abc_transporters_master.csv':
            continue

        file_path = os.path.join(INPUT_DIR, filename)
        # Determine if inhibitors or substrates file
        if '_inhibitors.csv' in filename:
            file_type = 'inhibitors'
        elif '_substrates.csv' in filename:
            file_type = 'substrates'
        else:
            # Unknown file type, skip
            continue
        
        # Extract target_chembl_id from filename
        # Format assumed: <target_chembl_id>_<type>.csv
        target_chembl_id = filename.split('_')[0]

        df = pd.read_csv(file_path)

        # Count raw rows
        raw_count = len(df)

        # Validate SMILES
        df['valid_smiles'] = df['smiles'].apply(is_valid_smiles)

        # For inhibitors: pchembl_value should not be None
        if file_type == 'inhibitors':
            cleaned_df = df[(df['valid_smiles'] == True) & (df['pchembl_value'].notna()) & (df['pchembl_value'] != '')]
        else:
            # For substrates: only SMILES validity mentioned
            cleaned_df = df[df['valid_smiles'] == True]

        cleaned_count = len(cleaned_df)

        # Save cleaned file
        cleaned_path = os.path.join(CLEANED_DIR, filename)
        cleaned_df.drop(columns=['valid_smiles'], inplace=True)
        cleaned_df.to_csv(cleaned_path, index=False)

        # Update stats
        if target_chembl_id not in stats:
            stats[target_chembl_id] = {
                'pref_name': target_name_map.get(target_chembl_id, ''),
                'inhibitors_count': 0,
                'substrates_count': 0,
                'cleaned_inhibitors_count': 0,
                'cleaned_substrates_count': 0
            }

        if file_type == 'inhibitors':
            stats[target_chembl_id]['inhibitors_count'] = raw_count
            stats[target_chembl_id]['cleaned_inhibitors_count'] = cleaned_count
        else:
            stats[target_chembl_id]['substrates_count'] = raw_count
            stats[target_chembl_id]['cleaned_substrates_count'] = cleaned_count

    # Produce overall summary
    records = []
    for tid, val in stats.items():
        records.append({
            'pref_name': val['pref_name'],
            'target_chembl_id': tid,
            'number_of_inhibitors': val['inhibitors_count'],
            'number_of_substrates': val['substrates_count'],
            'number_of_cleaned_inhibitors': val['cleaned_inhibitors_count'],
            'number_of_cleaned_substrates': val['cleaned_substrates_count']
        })

    summary_df = pd.DataFrame(records)
    summary_path = os.path.join(CLEANED_DIR, 'summary.csv')
    summary_df.to_csv(summary_path, index=False)

    print("Cleaning complete. Summary file saved at:", summary_path)

if __name__ == '__main__':
    main()
