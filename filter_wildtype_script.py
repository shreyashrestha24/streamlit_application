import numpy as np
import pandas as pd

# Load full dataset
baseline_df = pd.read_pickle(r"C:\Users\pz24886\Metabolic_Burden_Project\wildtype_1.pkl")
baseline_df = baseline_df.rename(columns={'gene_ko': 'simulation_identifier'})

# Columns you want to keep
columns_of_interest = [
    'simulation_identifier', 'time', 'countsToMolar', 'generation', 
    'growth', 'ATP[c]', 'GTP[c]', 'GLN[c]', 'cellMass', 'rnaMass', 'proteinMass', 'GLNS-MONOMER[c]', 'RELA-MONOMER[c]', 'RPOH-MONOMER[c]', 'RPOS-MONOMER[c]', 'ADP[c]'
]

# Get tRNA columns
tRNA_cols = [col for col in baseline_df.columns if 'tRNA' in col]
charged_tRNAs = [col for col in tRNA_cols if col.startswith('charged-')]
uncharged_tRNAs = [col for col in tRNA_cols if not col.startswith('charged-')]

# Filter to keep necessary columns
filtered_baseline = baseline_df[columns_of_interest + charged_tRNAs + uncharged_tRNAs]

# Compute charging_ratio per row (array of values)
def compute_trna_charging_ratio_per_row(row):
    charged_arrs = [np.array(row[col]).flatten() for col in charged_tRNAs]
    uncharged_arrs = [np.array(row[col]).flatten() for col in uncharged_tRNAs]

    charged_sum = np.sum(charged_arrs, axis=0)
    uncharged_sum = np.sum(uncharged_arrs, axis=0)
    total_trna = charged_sum + uncharged_sum

    ratio = np.divide(charged_sum, total_trna, out=np.zeros_like(charged_sum, dtype=float), where=total_trna != 0)
    return ratio

filtered_baseline['charging_ratio'] = filtered_baseline.apply(compute_trna_charging_ratio_per_row, axis=1)

print(filtered_baseline.columns.tolist())

# Save filtered baseline with charging_ratio arrays included
filtered_baseline.to_pickle("filtered_baseline_with_trna.pkl")

