import os
import sys
import glob
import natsort
import numpy as np
import pandas as pd

class Logger:
    def __init__(self, logfile):
        self.terminal = sys.stdout
        self.log = open(logfile, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()  # Ensure the message is written immediately

    def flush(self):
        pass

def split_sdf(sdf_file, ref_dir, output_dir):
    """
    Split the SDF file into multiple files, each containing one ligand.

    Parameters
    ----------
    sdf_file : str
        Path to the input SDF file.
    ref_dir: str
        The directory where the index CSV files are located.
    output_dir : str
        Directory where the output files will be saved.

    Returns
    -------
    data_dict: dict
        A dictionary containing keys as the paths to the output files and values as the delta_G values.
    """
    data_dict = {}
    protein_base = os.path.basename(os.path.dirname(sdf_file))
    prefix = os.path.basename(sdf_file).replace("_ligands.sdf", "")
    
    # Find the reference CSV file
    if 'bace_ciordia_retro_ligands.sdf' in sdf_file:
        # The exception that cannot be handled by the code below ...
        ref_csv = os.path.join(ref_dir, protein_base, f"retrospective_custom_core_extra_pkacorr_out.csv")
    else:
        ref_csv = os.path.join(ref_dir, protein_base, f"{prefix}_out.csv")
        if not os.path.exists(ref_csv):
            parts = prefix.split("_")
            for i in range(len(parts), 0, -1):
                new_prefix = "_".join(parts[:i])
                ref_csv = glob.glob(os.path.join(ref_dir, protein_base, f"{new_prefix}*.csv"))
                if len(ref_csv) > 0:
                    break

            if len(ref_csv) == 0:
                all_csv = glob.glob(os.path.join(ref_dir, protein_base, "*.csv"))
                for csv_file in all_csv:
                    csv_parts = os.path.basename(csv_file).split(".csv")[0].split("_")
                    # If there is intersection between parts and csv_parts, then we have a match
                    if len(set(parts).intersection(csv_parts)) > 0:
                        ref_csv = [csv_file]
                        break
                if len(ref_csv) == 0:
                    print(f"Note: No reference CSV file found for {sdf_file}. Skipping ...")
                    return 
            if len(ref_csv) > 1:
                raise ValueError(f"Multiple reference CSV files found for {sdf_file}: {ref_csv}")
            ref_csv = ref_csv[0]

    df = pd.read_csv(ref_csv)
    if 'Ligand name' in df.columns:
        ligand_name_to_dg = dict(zip(df['Ligand name'], df['Exp. dG (kcal/mol)']))
    elif 'Ligand' in df.columns:
        ligand_name_to_dg = dict(zip(df['Ligand'], df['Exp. dG (kcal/mol)']))
    else:
        raise KeyError(f"Expected columns 'Ligand name' or 'Ligand' not found in {ref_csv}.")

    # Make sure all keys are strings (some may be integers)
    ligand_name_to_dg = {str(k): v for k, v in ligand_name_to_dg.items()}

    with open(sdf_file, "r") as f:
        content = f.read()
    molecules = content.strip().split("$$$$\n")
    molecules = [mol.strip() for mol in molecules if mol.strip()]  # Remove any possible empty entries

    for i, molecule in enumerate(molecules):
        lines = molecule.strip().split('\n')

        # Check if the ligand has a binding affinity measurement
        ligand_name = lines[0].strip()
        if ligand_name not in ligand_name_to_dg:
            if 'hsp90_frag_2rings' in sdf_file:
                print()
                print(ligand_name)
                print(type(ligand_name))
                print(ligand_name_to_dg)
                print(ligand_name in ligand_name_to_dg)
                
            print(f"Note: Ligand {ligand_name} in {sdf_file} does not have a binding affinity measurement in {ref_csv}. Skipping ...")
            continue

        delta_G = ligand_name_to_dg[ligand_name]
        suffix = lines[0].replace(' ', '_').replace(',', '').replace('/', '_').replace('(', '').replace(')', '')
        output_name = os.path.join(output_dir, f"{protein_base}_{prefix}_{suffix}.sdf")

        with open(output_name, "w") as f:
            f.write(molecule)
            if not molecule.endswith('$$$$'):
                f.write('\n$$$$\n')
        
        # print('Created file:', output_name)
        data_dict[output_name] = delta_G

    return data_dict

def extract_delta_G(sdf_file):
    """
    Given a processed SDF file, extract the delta_G value from the SDF file.

    Parameters
    ----------
    sdf_file : str
        Path to the input SDF file.

    Returns
    -------
    delta_G : float
        The free energy change extracted from the SDF file.
    """
    f = open(sdf_file, 'r')
    lines = f.readlines()
    f.close()

    delta_G = None
    for i, line in enumerate(lines):
        if line == "> <r_exp_dg>\n":
            delta_G = float(lines[i + 1].strip())
            break
    
    if delta_G is None:
        raise ValueError(f"Delta G not found in {sdf_file}.")
    
    return delta_G

def delta_G_to_pK(delta_G):
    """
    Convert delta_G (in kcal/mol) to pK value.

    Parameters
    ----------
    delta_G : float
        The free energy change in kcal/mol.
    
    Returns
    -------
    pK : float
        The pK value calculated from the delta_G.
    """
    R = 1.987e-3  # kcal/(mol*K)
    T = 297  # Kelvin
    pK = -1 / np.log(10) * (delta_G / (R * T))
    return pK

if __name__ == "__main__":
    sys.stdout = Logger('process_benchmark.log')
    sys.stderr = Logger('process_benchmark.log')
    data_dir = os.path.join('.', "fep_benchmark_inputs", "structure_inputs")
    ref_dir = os.path.join('.', "21_4_results", "ligand_predictions")
    protein_dirs = [d for d in natsort.natsorted(glob.glob(os.path.join(data_dir, "*"))) if os.path.isdir(d)]

    data = []
    for protein_dir in protein_dirs:
        split_sdf_dir = os.path.join(protein_dir, "split_sdf")
        if not os.path.exists(split_sdf_dir):
            os.makedirs(split_sdf_dir)
        sdf_files = glob.glob(os.path.join(protein_dir, "*.sdf"))
        for sdf_file in sdf_files:
            data_dict = split_sdf(sdf_file, ref_dir, split_sdf_dir)
            if data_dict is None:
                continue
            processed_sdfs = data_dict.keys()
            for processed_sdf in processed_sdfs:
                system_id = f'{os.path.basename(processed_sdf).replace(".sdf", "")}'
                protein_path = sdf_file.replace('ligands.sdf', 'protein.pdb')
                delta_G = data_dict[processed_sdf]
                pK = delta_G_to_pK(delta_G)

                data.append({
                    'system_id': system_id,
                    'pK': pK,
                    'delta_G': delta_G,
                    'protein_path': protein_path,
                    'sdf_path': processed_sdf,
                })

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(data_dir, "processed_benchmark.csv"), index=False)
    print(f"Data saved to {os.path.join(data_dir, 'processed_benchmark.csv')}")
