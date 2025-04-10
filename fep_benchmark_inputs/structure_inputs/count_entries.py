import os
import glob
import pandas as pd

if __name__ == "__main__":
    subset_dict = {}  # key: subset name, value: number of entries
    all_dirs = [i for i in glob.glob('*') if os.path.isdir(i)]
    for dir_path in all_dirs:
        n_entries = 0
        all_csv = [i for i in glob.glob(dir_path + '/*.csv') if 'subset_metadata.csv' not in i]
        for csv_path in all_csv:
            # print(csv_path)
            df = pd.read_csv(csv_path)
            try:
                n_entries += len(set(df['Ligand 1']).union(set(df['Ligand 2'])))
            except KeyError:
                n_entries += len(set(df['Ligand1']).union(set(df['Ligand2'])))
        subset_dict[dir_path] = n_entries
    
    for dataset, n_entries in subset_dict.items():
        print(f"{dataset}: {n_entries}")
    print(f"\nTotal number of entries: {sum(subset_dict.values())}\n")

    # Get the number from subset_metadata.csv
    n_total = 0
    for dir_path in all_dirs:
        subset_metadata = pd.read_csv(dir_path + '/subset_metadata.csv')
        n = subset_metadata['Number of nodes'].sum()
        print(f"{dir_path}: {n}")
        n_total += n
    
    print(f"\nTotal number of entries: {n_total}")

    # Get the number from benchmark_metadata.csv
    df_benchmark = pd.read_csv('../benchmark_metadata.csv')
    n_total = df_benchmark['Number of nodes'].sum()
    print(f"\nTotal number of entries: {n_total}")

    # Check overlap with PDBbind/HiQBind
    fep_entries = df_benchmark['Reference PDB'].unique()
    fep_entries = set([i.lower() for i in fep_entries])

    pdbbind_entries = pd.read_csv('/home/bioc1870/Data/PDBbind_Opt/pdbbind_opt.csv')['PDBID'].unique()
    hiqbind_entries = pd.read_csv('/home/bioc1870/Data/HiQBind/hiqbind_metadata.csv')['PDBID'].unique()

    fep_pdbbind_overlap = fep_entries.intersection(set(pdbbind_entries))
    fep_hiqbind_overlap = fep_entries.intersection(set(hiqbind_entries))
    print(f"Overlap with PDBbind: {len(fep_pdbbind_overlap)} (out of {len(fep_entries)})")
    print(f"Overlap with HiQBind: {len(fep_hiqbind_overlap)} (out of {len(fep_entries)})")

    pdbbind_hiqbind_overlap = set(pdbbind_entries).intersection(set(hiqbind_entries))
    print(f"Overlap between PDBbind and HiQBind: {len(pdbbind_hiqbind_overlap)}")
