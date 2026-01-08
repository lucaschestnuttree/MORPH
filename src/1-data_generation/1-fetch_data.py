import os
import sys
import gzip
import shutil
import requests


def preprocess_file(file_path: str) -> None:
    """Remove the first three lines from a file."""
    temp_file_path = file_path + ".tmp"
    with open(file_path, 'r') as src_file, open(temp_file_path, 'w') as temp_file:
        # Skip the first three lines
        for _ in range(3):
            try:
                next(src_file)
            except StopIteration:
                # File has fewer than 3 lines; just stop
                return
        # Write the rest of the file
        for line in src_file:
            temp_file.write(line)
    os.replace(temp_file_path, file_path)


def missense_db_file_exists(file_name: str, missense_db_folder: str) -> str:
    """
    Checks if a Missense DB file exists in the specified folder and its unzipped state.
    
    Returns:
        'unzipped'  -> TSV (or uncompressed file) exists.
        'zipped'    -> Only .gz file exists.
        'missing'   -> Neither exists.
    """
    file_name_zip = file_name if file_name.endswith('.gz') else f"{file_name}.gz"
    file_name_uncompressed = file_name.replace('.gz', '')
    zip_path = os.path.join(missense_db_folder, file_name_zip)
    uncompressed_path = os.path.join(missense_db_folder, file_name_uncompressed)
    
    if os.path.isfile(uncompressed_path):
        return 'unzipped'
    elif os.path.isfile(zip_path):
        return 'zipped'
    return 'missing'


def fetch_missense_db(record_id: str, missense_db_dir: str) -> None:
    """Fetch all files from the given Zenodo record, save, unzip, and preprocess them."""
    url = f"https://zenodo.org/api/records/{record_id}"
    os.makedirs(missense_db_dir, exist_ok=True)

    try:
        response = requests.get(url)
        response.raise_for_status()
        record_metadata = response.json()
    except requests.HTTPError as e:
        print(f"[ERROR] Fetching Missense DB metadata from Zenodo: {e}")
        return

    # Zenodo API may structure files differently; handle common patterns
    files = record_metadata.get('files') or record_metadata.get('files', {}).get('entries', [])
    if isinstance(files, dict):
        files = files.values()

    for file_info in files:
        file_name = file_info['key']
        file_url = file_info['links']['self']

        status = missense_db_file_exists(file_name, missense_db_dir)
        zip_file_path = os.path.join(missense_db_dir, file_name)

        if status == 'unzipped':
            print(f"[SKIP] {file_name} already exists (unzipped).")
            continue
        elif status == 'zipped':
            print(f"[INFO] {file_name} exists but is still zipped.")
        else:
            print(f"[DOWNLOAD] {file_name}...")
            try:
                file_response = requests.get(file_url)
                file_response.raise_for_status()
                with open(zip_file_path, 'wb') as file:
                    file.write(file_response.content)
            except requests.HTTPError as e:
                print(f"[ERROR] Failed to download {file_name}: {e}")
                continue

        # Unzip .gz files
        if file_name.endswith('.gz'):
            unzipped_file_path = zip_file_path[:-3]  # strip .gz
            try:
                with gzip.open(zip_file_path, 'rb') as f_in:
                    with open(unzipped_file_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                print(f"[UNZIP] {file_name} -> {unzipped_file_path}")
                preprocess_file(unzipped_file_path)

                os.remove(zip_file_path)
                print(f"[CLEANUP] Removed {zip_file_path}")

            except Exception as e:
                print(f"[ERROR] Unzipping {file_name}: {e}")
        else:
            print(f"[PREPROCESS] {file_name}")
            preprocess_file(zip_file_path)


if __name__ == "__main__":
    # Usage:
    #   python script.py BASE_PATH [MISSENSE_DIR] [RECORD_ID]
    #
    # Defaults:
    #   MISSENSE_DIR -> data/1-raw/missense_db
    #   RECORD_ID    -> 10813168

    if len(sys.argv) < 2:
        print("Usage: python script.py BASE_PATH [MISSENSE_DIR] [RECORD_ID]")
        print("  BASE_PATH    : base directory (script will chdir into it)")
        print("  MISSENSE_DIR : directory for Missense DB files (default: data/1-raw/missense_db)")
        print("  RECORD_ID    : Zenodo record ID (default: 10813168)")
        sys.exit(1)

    base_path = sys.argv[1]
    missense_db_folder_path = sys.argv[2] if len(sys.argv) >= 3 else "data/1-raw/missense_db"
    record_id = sys.argv[3] if len(sys.argv) >= 4 else "10813168"

    os.chdir(base_path)
    print(f"[INFO] Working directory: {os.getcwd()}")
    print(f"[INFO] Missense DB folder: {missense_db_folder_path}")
    print(f"[INFO] Zenodo record ID: {record_id}")

    print("[RUN] Fetching Missense DB files...")
    fetch_missense_db(record_id=record_id, missense_db_dir=missense_db_folder_path)
    print("[DONE] Missense DB fetching completed.")
