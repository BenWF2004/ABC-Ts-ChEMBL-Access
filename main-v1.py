import os
import sys
import time
import logging
import pandas as pd
from chembl_webresource_client.new_client import new_client

# -------------------------- Configuration --------------------------
SPECIES = 'Homo sapiens'
INHIBITOR_ACTIVITY_TYPES = ['IC50', 'Ki', 'Kd', 'IC90', 'EC50']
SUBSTRATE_ACTIVITY_TYPES = ['Km', 'Vmax', 'Kcat']

OUTPUT_DIR = 'ABC-Ts'
CSV_DIR = os.path.join(OUTPUT_DIR, 'csv')
LOG_FILE = os.path.join(OUTPUT_DIR, 'chembl_api.log')

# Create output directories if they don't exist
os.makedirs(CSV_DIR, exist_ok=True)

# Setup logging
logging.basicConfig(
    filename=LOG_FILE,
    filemode='w',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logging.getLogger().addHandler(console_handler)

# -------------------------- Helper Functions --------------------------

def fetch_all_human_targets():
    """Fetch all human targets from ChEMBL."""
    logging.info("Fetching all human targets from ChEMBL...")
    target = new_client.target
    results = target.filter(organism=SPECIES)
    targets = list(results)
    logging.info(f"Total human targets retrieved: {len(targets)}")
    return targets

def target_has_abc_keyword(target_data):
    """
    Check if the target is an ABC transporter by:
    - Checking 'ABC' in pref_name
    - Checking 'ABC' in description
    - Checking 'ABC' in component synonyms
    """
    name = target_data.get('pref_name', '') or ''
    desc = target_data.get('description', '') or ''

    if 'ABC' in name or 'ABC' in desc:
        return True

    # Check component synonyms
    components = target_data.get('target_components', [])
    for comp in components:
        synonyms = comp.get('target_component_synonyms', [])
        for syn in synonyms:
            syn_name = syn.get('component_synonym', '')
            if 'ABC' in syn_name:
                return True

    return False

def filter_abc_transporters(all_targets):
    """Filter targets to include only ABC transporters."""
    logging.info("Filtering for ABC transporters among all human targets...")
    abc_transporters = [t for t in all_targets if target_has_abc_keyword(t)]
    logging.info(f"Total ABC transporters found: {len(abc_transporters)}")
    return abc_transporters

def save_master_abc_list(abc_transporters):
    """Save a master CSV of ABC transporters with relevant details."""
    if not abc_transporters:
        logging.info("No ABC transporters to save.")
        return None

    # Extract relevant fields
    records = []
    for t in abc_transporters:
        row = {
            'target_chembl_id': t.get('target_chembl_id', ''),
            'pref_name': t.get('pref_name', ''),
            'organism': t.get('organism', ''),
            'target_type': t.get('target_type', ''),
            'synonyms': '; '.join([
                syn.get('component_synonym', '')
                for comp in t.get('target_components', [])
                for syn in comp.get('target_component_synonyms', [])
            ])
        }
        records.append(row)

    df = pd.DataFrame(records)
    output_path = os.path.join(CSV_DIR, 'abc_transporters_master.csv')
    df.to_csv(output_path, index=False)
    logging.info(f"Master ABC transporters list saved to {output_path}")
    return df

def fetch_activity_data(target_chembl_id, activity_types):
    """
    Fetch activity data for a given target and a list of activity types.
    Returns a DataFrame of activities.
    """
    activity = new_client.activity
    activities = []
    try:
        query = activity.filter(target_chembl_id=target_chembl_id, standard_type__in=activity_types)
        for act in query:
            activity_info = {
                'activity_id': act.get('activity_id', ''),
                'assay_chembl_id': act.get('assay_chembl_id', ''),
                'molecule_chembl_id': act.get('molecule_chembl_id', ''),
                'standard_type': act.get('standard_type', ''),
                'standard_value': act.get('standard_value', ''),
                'standard_units': act.get('standard_units', ''),
                'standard_relation': act.get('standard_relation', ''),
                'published_value': act.get('published_value', ''),
                'published_units': act.get('published_units', ''),
                'published_relation': act.get('published_relation', ''),
                'pchembl_value': act.get('pchembl_value', ''),
                'activity_comment': act.get('activity_comment', ''),
                'data_validity_comment': act.get('data_validity_comment', ''),
                'confidence_score': act.get('confidence_score', ''),
                'ligand_efficiency': act.get('ligand_efficiency', {}),
                'document_chembl_id': act.get('document_chembl_id', ''),
                'record_id': act.get('record_id', ''),
                'target_chembl_id': act.get('target_chembl_id', ''),
            }
            activities.append(activity_info)
    except Exception as e:
        logging.error(f"Error fetching activity data for {target_chembl_id}: {e}")

    df = pd.DataFrame(activities)
    return df

def process_abc_transporters(abc_df):
    """
    For each ABC transporter, fetch inhibitor and substrate activities, 
    and save them as separate CSV files.
    """
    if abc_df is None or abc_df.empty:
        logging.info("No ABC transporters to process for activity data.")
        return

    total_targets = len(abc_df)
    start_time = time.time()

    for i, row in abc_df.iterrows():
        target_chembl_id = row['target_chembl_id']
        pref_name = row['pref_name']

        completed = i
        remaining = total_targets - completed
        elapsed = time.time() - start_time
        avg_time_per_target = elapsed / completed if completed > 0 else 0
        eta = avg_time_per_target * remaining if avg_time_per_target > 0 else 0
        eta_str = time.strftime("%H:%M:%S", time.gmtime(eta))

        logging.info(f"Processing {pref_name} ({target_chembl_id}) "
                     f"[{completed+1}/{total_targets}] - ETA: {eta_str}")

        # Fetch inhibitors
        df_inhibitors = fetch_activity_data(target_chembl_id, INHIBITOR_ACTIVITY_TYPES)
        if not df_inhibitors.empty:
            inhibitors_path = os.path.join(CSV_DIR, f"{target_chembl_id}_inhibitors.csv")
            df_inhibitors.to_csv(inhibitors_path, index=False)
            logging.info(f"Inhibitors for {pref_name} saved to {inhibitors_path}")

        # Fetch substrates
        df_substrates = fetch_activity_data(target_chembl_id, SUBSTRATE_ACTIVITY_TYPES)
        if not df_substrates.empty:
            substrates_path = os.path.join(CSV_DIR, f"{target_chembl_id}_substrates.csv")
            df_substrates.to_csv(substrates_path, index=False)
            logging.info(f"Substrates for {pref_name} saved to {substrates_path}")

    logging.info("Processing of ABC transporters completed.")


# -------------------------- Main Execution Flow --------------------------

def main():
    try:
        print("Starting to fetch and process ABC transporters. Please wait...")
        all_human_targets = fetch_all_human_targets()
        abc_transporters = filter_abc_transporters(all_human_targets)
        abc_df = save_master_abc_list(abc_transporters)

        if abc_df is not None and not abc_df.empty:
            process_abc_transporters(abc_df)
            print("All tasks completed successfully. Check the log file and CSV directory for results.")
        else:
            print("No ABC transporters found.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        print("An error occurred. Please check the log file for details.")

if __name__ == '__main__':
    main()
