import os
import logging
from typing import List, Dict, Any
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from chembl_webresource_client.new_client import new_client

# -------------------------- Configuration --------------------------

# Define target ChEMBL IDs directly to avoid synonym issues
# These correspond to ABCB1 and ABCG2 respectively
INITIAL_TARGETS = ['CHEMBL2091', 'CHEMBL1930']  # ABCB1 and ABCG2 respectively
SPECIES = ['Homo sapiens']                      # Start with human; can add 'Mus musculus', etc.

# Define activity types for inhibitors and substrates
INHIBITOR_ACTIVITY_TYPES = ['IC50', 'Ki', 'Kd', 'IC90', 'EC50']
SUBSTRATE_ACTIVITY_TYPES = ['Km', 'Vmax', 'Kcat']

# Output directories
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
    """
    Extracts the family from the preferred name of the ABC transporter.
    For example, 'ABCB1' -> 'ABC-B'.
    """
    match = re.match(r'ABC([A-Z])\d+', pref_name)
    if match:
        return f"ABC-{match.group(1)}"
    else:
        return "Unknown"

def fetch_target_info(target_ids: List[str], species: List[str]) -> pd.DataFrame:
    """
    Fetch target information from ChEMBL based on target IDs and species.
    """
    logging.info(f"Fetching target information for target IDs: {target_ids} and species: {species}")
    target = new_client.target
    targets_data = []

    try:
        for target_id in target_ids:
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
        df_targets = pd.DataFrame(targets_data)
        df_targets.drop_duplicates(inplace=True)
        logging.info(f"Retrieved {len(df_targets)} target entries.")
        return df_targets
    except Exception as e:
        logging.error(f"Error fetching target information: {e}")
        raise

def fetch_all_abc_ts(species: List[str]) -> List[Dict[str, Any]]:
    """
    Fetch all ATP-Binding Cassette Transporters (ABC-Ts) from ChEMBL.
    """
    logging.info(f"Fetching all ABC-Ts for species: {species}")
    target = new_client.target
    abc_ts = []

    try:
        # Define search criteria for ABC-Ts
        # Filter targets where target_synonym contains 'ABC'
        query = target.filter(target_synonym__icontains='ABC')
        for res in query:
            if res['organism'] in species:
                family = extract_family(res['pref_name'])
                abc_t_info = {
                    'target_chembl_id': res['target_chembl_id'],
                    'pref_name': res['pref_name'],
                    'family': family,
                    'organism': res['organism']
                }
                abc_ts.append(abc_t_info)
        logging.info(f"Retrieved {len(abc_ts)} ABC-Ts.")
        return abc_ts
    except Exception as e:
        logging.error(f"Error fetching ABC-Ts: {e}")
        return abc_ts

def fetch_activity_data(target_chembl_id: str, activity_types: List[str], activity_category: str) -> pd.DataFrame:
    """
    Fetch activity data (inhibitors or substrates) for a given target from the activity endpoint.
    """
    logging.info(f"Fetching {activity_category} activities for target: {target_chembl_id}")
    activity = new_client.activity
    activities = []

    try:
        # Fetch activities based on provided activity types
        query = activity.filter(target_chembl_id=target_chembl_id, standard_type__in=activity_types)
        for act in query:
            molecule_chembl_id = act.get('molecule_chembl_id', '')
            molecule_record = new_client.molecule.get(molecule_chembl_id)

            if molecule_record is None:
                logging.warning(f"Molecule with ChEMBL ID {molecule_chembl_id} not found.")
                continue

            max_phase = molecule_record.get('max_phase')
            if isinstance(max_phase, (int, float)):
                if max_phase > 0:
                    source = 'Approved Drug'
                else:
                    source = 'Experimental'
            else:
                source = 'Experimental'  # Default if max_phase is None or invalid

            activity_info = {
                'molecule_chembl_id': molecule_chembl_id,
                'pref_name': molecule_record.get('pref_name', ''),
                'canonical_smiles': molecule_record.get('molecule_structures', {}).get('canonical_smiles', ''),
                'activity_type': act.get('standard_type', ''),
                'standard_value': act.get('standard_value', None),
                'standard_units': act.get('standard_units', ''),
                'confidence_score': act.get('confidence_score', None),
                'max_phase': max_phase,
                'source': source,
                'target': molecule_record.get('target', {}).get('target_chembl_id', ''),
                'organism': molecule_record.get('organism', '')
            }
            activities.append(activity_info)

        df_activities = pd.DataFrame(activities)
        logging.info(f"Retrieved {len(df_activities)} {activity_category} activities for target {target_chembl_id}.")
        return df_activities
    except Exception as e:
        logging.error(f"Error fetching {activity_category} activities for target {target_chembl_id}: {e}")
        return pd.DataFrame()

def process_and_save_data(df_targets: pd.DataFrame, output_dir: str) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Process activity data (inhibitors and substrates) for each target and save to CSV files.
    Returns a nested dictionary with target names as keys and their respective activity DataFrames.
    """
    logging.info("Processing and saving activity data.")
    target_data = {}

    for _, row in df_targets.iterrows():
        target_id = row['target_chembl_id']
        target_name = row['pref_name']
        family = row['family']
        organism = row['organism']

        try:
            # Fetch inhibitors
            df_inhibitors = fetch_activity_data(target_id, INHIBITOR_ACTIVITY_TYPES, 'inhibitor')
            # Fetch substrates
            df_substrates = fetch_activity_data(target_id, SUBSTRATE_ACTIVITY_TYPES, 'substrate')

            target_data[target_name] = {
                'inhibitors': df_inhibitors,
                'substrates': df_substrates
            }

            # Save inhibitors to CSV
            if not df_inhibitors.empty:
                df_inhibitors['family'] = family
                df_inhibitors['organism'] = organism
                # Handle missing standard values
                df_inhibitors['standard_value'] = pd.to_numeric(df_inhibitors['standard_value'], errors='coerce')

                # Define CSV filename
                csv_inhibitors = f"{target_name.replace(' ', '_')}_inhibitors.csv"
                csv_inhibitors_path = os.path.join(output_dir, csv_inhibitors)
                df_inhibitors.to_csv(csv_inhibitors_path, index=False)
                logging.info(f"Saved inhibitors data to {csv_inhibitors_path}.")
            else:
                logging.info(f"No inhibitors data found for target {target_name}.")

            # Save substrates to CSV
            if not df_substrates.empty:
                df_substrates['family'] = family
                df_substrates['organism'] = organism
                # Handle missing standard values
                df_substrates['standard_value'] = pd.to_numeric(df_substrates['standard_value'], errors='coerce')

                # Define CSV filename
                csv_substrates = f"{target_name.replace(' ', '_')}_substrates.csv"
                csv_substrates_path = os.path.join(output_dir, csv_substrates)
                df_substrates.to_csv(csv_substrates_path, index=False)
                logging.info(f"Saved substrates data to {csv_substrates_path}.")
            else:
                logging.info(f"No substrates data found for target {target_name}.")

        except Exception as e:
            logging.error(f"Failed to process data for target {target_name}: {e}")

    return target_data

def generate_summary(target_data: Dict[str, Dict[str, pd.DataFrame]], summary_file: str):
    """
    Generate a summary text file with key statistics for inhibitors and substrates.
    """
    logging.info("Generating summary file.")
    with open(summary_file, 'w') as f:
        for target, activities in target_data.items():
            f.write(f"=== {target} ===\n")
            
            # Inhibitors Summary
            df_inhibitors = activities['inhibitors']
            if not df_inhibitors.empty:
                total_inhibitors = len(df_inhibitors)
                f.write(f"-- Inhibitors --\n")
                f.write(f"Total Inhibitors: {total_inhibitors}\n")

                # Confidence scores distribution
                confidence_scores = df_inhibitors['confidence_score'].dropna()
                if not confidence_scores.empty:
                    f.write(f"Confidence Score - Mean: {confidence_scores.mean():.2f}, "
                            f"Median: {confidence_scores.median():.2f}, "
                            f"Std Dev: {confidence_scores.std():.2f}\n")
                else:
                    f.write("Confidence Score: No data available.\n")

                # Phase information
                phases = df_inhibitors['max_phase'].dropna().unique()
                phases_str = ', '.join(map(str, phases)) if len(phases) > 0 else 'No data'
                f.write(f"Max Phase: {phases_str}\n")
            else:
                f.write("-- Inhibitors --\n")
                f.write("No inhibitor data available.\n")

            # Substrates Summary
            df_substrates = activities['substrates']
            if not df_substrates.empty:
                total_substrates = len(df_substrates)
                f.write(f"-- Substrates --\n")
                f.write(f"Total Substrates: {total_substrates}\n")

                # Confidence scores distribution
                confidence_scores = df_substrates['confidence_score'].dropna()
                if not confidence_scores.empty:
                    f.write(f"Confidence Score - Mean: {confidence_scores.mean():.2f}, "
                            f"Median: {confidence_scores.median():.2f}, "
                            f"Std Dev: {confidence_scores.std():.2f}\n")
                else:
                    f.write("Confidence Score: No data available.\n")

                # Phase information
                phases = df_substrates['max_phase'].dropna().unique()
                phases_str = ', '.join(map(str, phases)) if len(phases) > 0 else 'No data'
                f.write(f"Max Phase: {phases_str}\n")
            else:
                f.write("-- Substrates --\n")
                f.write("No substrates data available.\n")

            f.write("\n")
    logging.info(f"Summary file created at {summary_file}.")

def create_plots(target_data: Dict[str, Dict[str, pd.DataFrame]], plots_dir: str):
    """
    Create and save plots for inhibitor and substrate activities and confidence scores.
    """
    logging.info("Creating plots.")
    for target, activities in target_data.items():
        # Inhibitors Plots
        df_inhibitors = activities['inhibitors']
        if not df_inhibitors.empty:
            # Inhibitor Activity Distribution
            if df_inhibitors['standard_value'].notna().any():
                plt.figure(figsize=(10, 6))
                sns.histplot(data=df_inhibitors, x='standard_value', hue='activity_type', kde=True, bins=30)
                plt.title(f'Inhibitor Activity Distribution for {target}')
                plt.xlabel('Activity Value')
                plt.ylabel('Frequency')
                plt.legend(title='Activity Type')
                affinity_plot_path = os.path.join(plots_dir, f"{target.replace(' ', '_')}_inhibitor_activity.png")
                plt.savefig(affinity_plot_path)
                plt.close()
                logging.info(f"Saved inhibitor activity plot to {affinity_plot_path}.")
            else:
                logging.info(f"No inhibitor activity data available for {target}; skipping activity plot.")

            # Inhibitor Confidence Score Distribution
            confidence_scores = df_inhibitors['confidence_score'].dropna()
            if not confidence_scores.empty:
                plt.figure(figsize=(10, 6))
                sns.histplot(data=df_inhibitors, x='confidence_score', kde=True, bins=30, color='blue')
                plt.title(f'Inhibitor Confidence Score Distribution for {target}')
                plt.xlabel('Confidence Score')
                plt.ylabel('Frequency')
                confidence_plot_path = os.path.join(plots_dir, f"{target.replace(' ', '_')}_inhibitor_confidence_score.png")
                plt.savefig(confidence_plot_path)
                plt.close()
                logging.info(f"Saved inhibitor confidence score plot to {confidence_plot_path}.")
            else:
                logging.info(f"No inhibitor confidence score data available for {target}; skipping confidence score plot.")
        else:
            logging.info(f"No inhibitors data available for {target}; skipping inhibitor plots.")

        # Substrates Plots
        df_substrates = activities['substrates']
        if not df_substrates.empty:
            # Substrate Activity Distribution
            if df_substrates['standard_value'].notna().any():
                plt.figure(figsize=(10, 6))
                sns.histplot(data=df_substrates, x='standard_value', hue='activity_type', kde=True, bins=30)
                plt.title(f'Substrate Activity Distribution for {target}')
                plt.xlabel('Activity Value')
                plt.ylabel('Frequency')
                plt.legend(title='Activity Type')
                affinity_plot_path = os.path.join(plots_dir, f"{target.replace(' ', '_')}_substrate_activity.png")
                plt.savefig(affinity_plot_path)
                plt.close()
                logging.info(f"Saved substrate activity plot to {affinity_plot_path}.")
            else:
                logging.info(f"No substrate activity data available for {target}; skipping activity plot.")

            # Substrate Confidence Score Distribution
            confidence_scores = df_substrates['confidence_score'].dropna()
            if not confidence_scores.empty:
                plt.figure(figsize=(10, 6))
                sns.histplot(data=df_substrates, x='confidence_score', kde=True, bins=30, color='green')
                plt.title(f'Substrate Confidence Score Distribution for {target}')
                plt.xlabel('Confidence Score')
                plt.ylabel('Frequency')
                confidence_plot_path = os.path.join(plots_dir, f"{target.replace(' ', '_')}_substrate_confidence_score.png")
                plt.savefig(confidence_plot_path)
                plt.close()
                logging.info(f"Saved substrate confidence score plot to {confidence_plot_path}.")
            else:
                logging.info(f"No substrate confidence score data available for {target}; skipping confidence score plot.")
        else:
            logging.info(f"No substrates data available for {target}; skipping substrate plots.")

def list_all_abc_ts(species: List[str]) -> List[Dict[str, Any]]:
    """
    Fetch and return a list of all ABC-Ts from ChEMBL for the specified species.
    """
    abc_ts = fetch_all_abc_ts(species)
    return abc_ts

def display_abc_ts(abc_ts: List[Dict[str, Any]]):
    """
    Display the list of ABC-Ts in the terminal as an array.
    """
    if not abc_ts:
        print("No ABC-Ts found.")
        return

    # Convert list of dicts to list of tuples for better readability
    abc_ts_list = [(t['target_chembl_id'], t['pref_name'], t['family'], t['organism']) for t in abc_ts]
    print("\nList of ATP-Binding Cassette Transporters (ABC-Ts):")
    for abc in abc_ts_list:
        print(abc)

def save_abc_ts_to_csv(abc_ts: List[Dict[str, Any]], output_dir: str):
    """
    Save the list of ABC-Ts into a CSV file, sorted by family and preferred name.
    """
    if not abc_ts:
        print("No ABC-Ts to save.")
        logging.info("No ABC-Ts available to save to CSV.")
        return

    df_abc_ts = pd.DataFrame(abc_ts)
    # Sort by family and then by preferred name
    df_abc_ts_sorted = df_abc_ts.sort_values(by=['family', 'pref_name'])
    # Define CSV path
    csv_path = os.path.join(output_dir, 'abc_ts_sorted.csv')
    # Save to CSV
    df_abc_ts_sorted.to_csv(csv_path, index=False)
    logging.info(f"Saved ABC-Ts to CSV at {csv_path}.")
    print(f"\nABC-Ts have been saved to {csv_path}.")

# -------------------------- Main Execution Flow --------------------------

def main():
    try:
        # Interactive Menu
        print("Choose an option:")
        print("1. List all ATP-Binding Cassette Transporters (ABC-Ts) from ChEMBL")
        print("2. Fetch inhibitor and substrate data for ABC-Ts")
        choice = input("Enter 'list' to list ABC-Ts or 'fetch' to retrieve inhibitor and substrate data: ").strip().lower()

        if choice == 'list' or choice == '1':
            # Option 1: List all ABC-Ts
            abc_ts = list_all_abc_ts(SPECIES)
            display_abc_ts(abc_ts)
            save_abc_ts_to_csv(abc_ts, CSV_DIR)
            logging.info("Displayed and saved all ABC-Ts to CSV.")
            print("\nProcess completed. Check the log file for more details if needed.")

        elif choice == 'fetch' or choice == '2':
            # Option 2: Fetch inhibitor and substrate data for ABC-Ts
            # Prompt user to input additional ChEMBL IDs if desired
            user_input = input("Enter ChEMBL IDs separated by commas (or press Enter to use default ABC-Ts): ").strip()
            if user_input:
                user_targets = [tid.strip() for tid in user_input.split(',') if tid.strip()]
                if not user_targets:
                    print("No valid ChEMBL IDs entered. Exiting.")
                    logging.warning("User entered invalid or empty ChEMBL IDs.")
                    return
                df_targets = fetch_target_info(user_targets, SPECIES)
            else:
                df_targets = fetch_target_info(INITIAL_TARGETS, SPECIES)

            if df_targets.empty:
                logging.warning("No target information retrieved. Exiting program.")
                print("No target information found. Please check the target names and species.")
                return

            # Display fetched targets
            print("\nFetched Targets:")
            print(df_targets[['target_chembl_id', 'pref_name', 'family', 'organism']])

            # Save target metadata to CSV
            targets_csv = os.path.join(CSV_DIR, 'target_metadata.csv')
            df_targets.to_csv(targets_csv, index=False)
            logging.info(f"Saved target metadata to {targets_csv}.")

            # Step 2 & 3: Fetch Inhibitor and Substrate Activity Data and Save to CSV
            target_data = process_and_save_data(df_targets, CSV_DIR)
            if not target_data:
                logging.warning("No activity data retrieved for any targets.")
                print("No activity data found for the specified targets.")
                return

            # Step 4: Generate Summary
            generate_summary(target_data, SUMMARY_FILE)

            # Step 5: Create Plots
            create_plots(target_data, PLOTS_DIR)

            logging.info("Program completed successfully.")
            print(f"\nData retrieval and processing complete. Check the '{OUTPUT_DIR}' directory for outputs.")

        else:
            print("Invalid choice. Please run the script again and enter 'list' or 'fetch'.")
            logging.warning(f"User entered invalid choice: {choice}")

    except Exception as e:
        logging.error(f"An error occurred in the main execution: {e}")
        print("An error occurred. Please check the log file for details.")

if __name__ == '__main__':
    main()
