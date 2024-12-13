import os
import logging
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Configuration
SPECIES = 'Homo sapiens'
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

def fetch_human_targets():
    """Fetch all human targets from ChEMBL."""
    target = new_client.target
    human_targets = target.filter(organism=SPECIES)
    return human_targets

def filter_abc_transporters(targets):
    """Filter targets to include only ABC transporters."""
    abc_transporters = []
    for t in targets:
        if 'ABC' in t.get('pref_name', '') or 'ABC' in t.get('description', ''):
            abc_transporters.append(t)
    return abc_transporters

def save_to_csv(data, output_path):
    """Save data to CSV."""
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)

def main():
    try:
        logging.info("Fetching human targets from ChEMBL...")
        human_targets = fetch_human_targets()
        logging.info(f"Total human targets retrieved: {len(human_targets)}")

        logging.info("Filtering for ABC transporters...")
        abc_transporters = filter_abc_transporters(human_targets)
        logging.info(f"Total ABC transporters found: {len(abc_transporters)}")

        if abc_transporters:
            output_path = os.path.join(CSV_DIR, 'abc_transporters.csv')
            save_to_csv(abc_transporters, output_path)
            print(f"ABC transporters have been saved to {output_path}.")
        else:
            print("No ABC transporters found.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        print("An error occurred. Please check the log file for details.")

if __name__ == '__main__':
    main()
