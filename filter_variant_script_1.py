import pandas as pd

file_path = r"C:\Users\pz24886\Metabolic_Burden_Project\bpsA_param_sweep_2.pkl"
runs = pd.read_pickle(file_path)

# Columns you absolutely need
core_columns = [
    'generation', 'time', 'cellMass', 'rnaMass', 'proteinMass',
    'growth', 'ATP[c]', 'countsToMolar', 'GLN[c]', 'NG-BpsA-MONOMER[c]',
    'strain_id', 'repeat', 'IND[c]'
]

# Identify tRNA columns
tRNA_cols = [col for col in runs.columns if 'tRNA' in col]
charged_tRNAs = [col for col in tRNA_cols if col.startswith('charged-')]
uncharged_tRNAs = [col for col in tRNA_cols if not col.startswith('charged-')]

# Keep only the filtered columns
filtered_runs = runs[core_columns + charged_tRNAs + uncharged_tRNAs]

# Save filtered dataframe to pickle (smaller file)
filtered_runs.to_pickle("filtered_runs.pkl")
