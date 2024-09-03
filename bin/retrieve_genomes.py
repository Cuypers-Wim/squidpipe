#!/usr/bin/env python3

"""
retrieve_genomes.py

This script processes a CSV file containing taxon IDs and sample names, attempts to download the most recent RefSeq reference genome for those taxon IDs.
If no RefSeq reference can be found, the script will try to obtain the most complete Genbank .fna file. If no genome is found for a given taxon ID, it will attempt to download genomes for
its parent taxon IDs.

Usage:
    python retrieve_genomes.py -c <path_to_csv_file> -o <output_directory>

Setup:
    Ensure NCBI's 'datasets' software is installed and accessible in the environment.
"""

import csv
import os
import subprocess
import zipfile
import shutil
import logging
import json
import argparse

def setup_logging():
    """
    Set up logging to file and console.
    """
    logging.basicConfig(filename='retrieve_genomes.log', level=logging.INFO,
                        format='%(asctime)s %(levelname)s:%(message)s')

def log(message):
    """
    Log a message to both console and file.

    Parameters:
    message (str): The message to log.
    """
    print(message)
    logging.info(message)

def download_genome(taxon_id, filename, output_dir, assembly_source='refseq', assembly_version='latest', assembly_level=None):
    """
    Attempt to download the genome for a given taxon ID.

    Parameters:
    taxon_id (int): The taxon ID for which to download the genome.
    filename (str): The base filename for the downloaded genome.
    output_dir (str): The directory where the downloaded files will be stored.
    assembly_source (str): The source of the assembly (default 'refseq').
    assembly_version (str): The version of the assembly (default 'latest').
    assembly_level (str): The assembly level (default None).

    Returns:
    bool: True if the download was successful, False otherwise.
    """
    try:
        command = [
            'datasets', 'download', 'genome', 'taxon', str(taxon_id),
            '--filename', os.path.join(output_dir, f"{filename}.zip")
        ]
        if assembly_level:
            command.extend(['--assembly-level', assembly_level])
            log(f"Attempting to download the most complete GenBank genome for {filename} or a parent taxon")
        else:
            command.extend(['--assembly-source', assembly_source, '--assembly-version', assembly_version])
            log(f"Attempting to download the latest RefSeq genome for {filename} or a parent taxon")

        subprocess.run(command, check=True)
        return True
    except subprocess.CalledProcessError as e:
        log(f"Error downloading genome for taxon ID {taxon_id}: {e}")
        return False

def get_parent_taxon_ids(taxon_id):
    """
    Get the parent taxon IDs for a given taxon ID.

    Parameters:
    taxon_id (int): The taxon ID for which to get parent taxon IDs.

    Returns:
    list: A list of parent taxon IDs in reversed order.
    """
    result = subprocess.run([
        'datasets', 'summary', 'taxonomy', 'taxon', str(taxon_id)
    ], capture_output=True, check=True, text=True)
    summary = result.stdout

    # Parse the JSON response
    summary_json = json.loads(summary)

    # Navigate to the parents key
    try:
        parents = summary_json['reports'][0]['taxonomy']['parents']
        return parents[::-1]  # Return the list of parent taxon IDs in reversed order
    except (KeyError, IndexError):
        return []

def process_taxon(taxon_id, name, output_dir):
    """
    Process a given taxon ID to download and extract the genome.

    Parameters:
    taxon_id (int): The taxon ID to process.
    name (str): The name of the sample.
    output_dir (str): The directory where the processed files will be stored.

    Returns:
    bool: True if processing was successful, False otherwise.
    """
    if download_genome(taxon_id, name, output_dir):
        if extract_and_move_largest_fna(name, output_dir):
            return True

    # Try downloading with --assembly-level chromosome,complete
    if download_genome(taxon_id, name, output_dir, assembly_level='chromosome,complete'):
        if extract_and_move_largest_fna(name, output_dir):
            return True

    return False

def extract_and_move_largest_fna(name, output_dir):
    """
    Extract and move the largest .fna file from the downloaded zip file.

    Parameters:
    name (str): The name of the sample.
    output_dir (str): The directory where the extracted files will be stored.

    Returns:
    bool: True if extraction and moving was successful, False otherwise.
    """
    zip_path = os.path.join(output_dir, f"{name}.zip")
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(output_dir)
    fna_files = [
        os.path.join(root, file)
        for root, _, files in os.walk(os.path.join(output_dir, 'ncbi_dataset', 'data'))
        for file in files if file.endswith('.fna')
    ]
    if not fna_files:
        log(f"No .fna file found for {name}")
        return False

    # Move the largest or only .fna file
    largest_fna = max(fna_files, key=os.path.getsize)
    shutil.move(largest_fna, os.path.join(output_dir, f"{name}.fna"))
    return True

def main(csv_file, output_dir):
    """
    Main function to process the CSV file and download genomes.

    Parameters:
    csv_file (str): The path to the CSV file containing taxon IDs and sample names.
    output_dir (str): The directory where the processed files will be stored.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Open the CSV file and process each row
    with open(csv_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            taxon_id, name = row[2], row[1]
            try:
                log(f"Processing {taxon_id} ({name})")
                if not process_taxon(taxon_id, name, output_dir):
                    parent_taxon_ids = get_parent_taxon_ids(taxon_id)
                    for parent_id in parent_taxon_ids:
                        log(f"Attempting to download genome for parent taxon ID {parent_id}")
                        if process_taxon(parent_id, name, output_dir):
                            break
                    else:
                        log(f"No .fna file found for any parent taxon ID of {taxon_id} ({name})")
            except Exception as e:
                log(f"Error processing {taxon_id} ({name}): {e}")
            finally:
                zip_path = os.path.join(output_dir, f"{name}.zip")
                if os.path.exists(zip_path):
                    os.remove(zip_path)
                dataset_path = os.path.join(output_dir, 'ncbi_dataset')
                if os.path.exists(dataset_path):
                    shutil.rmtree(dataset_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a CSV file to download genomes based on taxon IDs.")
    parser.add_argument('-c', '--csv_file', required=True, help="Path to the CSV file containing taxon IDs and sample names.")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory where the downloaded and processed files will be stored.")
    args = parser.parse_args()

    setup_logging()
    log("Starting genome retrieval script")
    log(f"Current working directory: {os.getcwd()}")

    main(args.csv_file, args.output_dir)
    log("Genome retrieval script completed")
