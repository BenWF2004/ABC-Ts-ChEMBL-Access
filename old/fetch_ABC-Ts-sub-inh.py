import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from typing import List, Dict, Any
from chembl_webresource_client.new_client import new_client

# -------------------------- Configuration --------------------------

SPECIES = ['Homo sapiens']                      # Target species (Homo sapiens)
INHIBITOR_ACTIVITY_TYPES = ['IC50', 'Ki', 'Kd', 'IC90', 'EC50']
SUBSTRATE_ACTIVITY_TYPES = ['Km', 'Vmax', 'Kcat']

OUTPUT_DIR = 'chembl_output'
CSV_DIR = os.path.join(OUTPUT_DIR, 'csv')
PLOTS_DIR = os.path.join(OUTPUT_DIR, 'plots')
SUMMARY_FILE = os.path.join(OUTPUT_DIR, 'summary.txt')
LOG_FILE = os.path.join(OUTPUT_DIR, 'chembl_api.log')

# Create output directories if they don't exist
os.makedirs(CSV_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

# Setup logging
logging.basicConfig(
    filename=LOG_FILE,
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# -------------------------- Helper Functions --------------------------

def extract_family(pref_name: str) -> str:
    match = re.match(r'ABC([A-Z])\d+', pref_name)
    if match:
        return f"ABC-{match.group(1)}"
    else:
        return "Unknown"

def fetch_target_info(target_ids: List[str], species: List[str]) -> pd.DataFrame:
    logging.info(f"Fetching target information for target IDs: {target_ids} and species: {species}")
    target = new_client.target
    targets_data = []

    for target_id in target_ids:
        try:
            res = target.get(target_id)
            if res is None:
                logging.warning(f"Target with ChEMBL ID {target_id} not found.")
                continue
            if res['organism'] in species:
                family = extract_family(res['pref_name'])
                target_info = {
                    'target_chembl_id': res['target_chembl_id'],
                    'pref_name': res['pref_name'],
                    'family': family,
                    'organism': res['organism'],
                    'target_type': res['target_type']
                }
                targets_data.append(target_info)
        except Exception as e:
            logging.error(f"Error fetching target information for {target_id}: {e}")
    
    df_targets = pd.DataFrame(targets_data)
    if not df_targets.empty:
        df_targets.drop_duplicates(inplace=True)
    logging.info(f"Retrieved {len(df_targets)} target entries.")
    return df_targets

def fetch_activity_data(target_chembl_id: str, activity_types: List[str], activity_category: str) -> pd.DataFrame:
    logging.info(f"Fetching {activity_category} activities for target: {target_chembl_id}")
    activity = new_client.activity
    activities = []

    try:
        query = activity.filter(target_chembl_id=target_chembl_id, standard_type__in=activity_types)
        for act in query:
            activity_info = {
                'molecule_chembl_id': act.get('molecule_chembl_id', ''),
                'activity_type': act.get('standard_type', ''),
                'standard_value': act.get('standard_value', None),
                'standard_units': act.get('standard_units', ''),
                'confidence_score': act.get('confidence_score', None)
            }
            activities.append(activity_info)
    except Exception as e:
        logging.error(f"Error fetching {activity_category} activities for target {target_chembl_id}: {e}")
    
    df_activities = pd.DataFrame(activities)
    logging.info(f"Retrieved {len(df_activities)} {activity_category} activities for target {target_chembl_id}.")
    return df_activities

def process_and_save_data(df_targets: pd.DataFrame) -> Dict[str, Dict[str, pd.DataFrame]]:
    logging.info("Processing and saving activity data.")
    target_data = {}

    for _, row in df_targets.iterrows():
        target_id = row['target_chembl_id']
        target_name = row['pref_name']

        try:
            df_inhibitors = fetch_activity_data(target_id, INHIBITOR_ACTIVITY_TYPES, 'inhibitor')
            df_substrates = fetch_activity_data(target_id, SUBSTRATE_ACTIVITY_TYPES, 'substrate')

            target_data[target_name] = {
                'inhibitors': df_inhibitors,
                'substrates': df_substrates
            }

            if not df_inhibitors.empty:
                csv_path = os.path.join(CSV_DIR, f"{target_name.replace(' ', '_')}_inhibitors.csv")
                df_inhibitors.to_csv(csv_path, index=False)
                logging.info(f"Saved inhibitors data to {csv_path}.")

            if not df_substrates.empty:
                csv_path = os.path.join(CSV_DIR, f"{target_name.replace(' ', '_')}_substrates.csv")
                df_substrates.to_csv(csv_path, index=False)
                logging.info(f"Saved substrates data to {csv_path}.")
        except Exception as e:
            logging.error(f"Failed to process data for target {target_name}: {e}")

    return target_data

def get_chembl_ids_from_csv(csv_path: str) -> List[str]:
    try:
        df = pd.read_csv(csv_path)
        if 'target_chembl_id' not in df.columns:
            logging.error(f"CSV file {csv_path} does not contain a 'target_chembl_id' column.")
            return []
        return df['target_chembl_id'].dropna().unique().tolist()
    except Exception as e:
        logging.error(f"Error reading CSV file {csv_path}: {e}")
        return []

# -------------------------- Main Execution Flow --------------------------

def main():
    try:
        print("\nHow would you like to provide ChEMBL IDs?")
        print("1. Enter ChEMBL IDs manually (comma-separated)")
        print("2. Use a CSV file with a 'target_chembl_id' column")
        choice = input("Enter 1 or 2: ").strip()

        if choice == '1':
            user_input = input("Enter ChEMBL IDs separated by commas: ").strip()
            if user_input:
                target_ids = [tid.strip() for tid in user_input.split(',') if tid.strip()]
            else:
                print("No ChEMBL IDs entered.")
                return

        elif choice == '2':
            csv_path = input("Enter the path to the CSV file: ").strip()
            if not os.path.isfile(csv_path):
                print(f"File {csv_path} not found.")
                return
            target_ids = get_chembl_ids_from_csv(csv_path)
            if not target_ids:
                print("No valid ChEMBL IDs found in the CSV file.")
                return

        else:
            print("Invalid choice. Please enter 1 or 2.")
            return

        df_targets = fetch_target_info(target_ids, SPECIES)
        if df_targets.empty:
            print("No target information found. Please check the ChEMBL IDs.")
            return

        print("\nFetched Targets:")
        print(df_targets[['target_chembl_id', 'pref_name']])
        
        target_data = process_and_save_data(df_targets)
        if target_data:
            print("Data processing complete. Check the CSV directory for outputs.")
        else:
            print("No data found for the specified targets.")
    except Exception as e:
        logging.error(f"An error occurred during the main execution: {e}")
        print("An error occurred. Please check the log file for details.")

if __name__ == '__main__':
    main()
