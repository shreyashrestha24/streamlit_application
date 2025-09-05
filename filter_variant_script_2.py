import numpy as np
import pandas as pd

# Load the filtered DataFrame from part1
file_path = "C:\\Users\\pz24886\\Metabolic_Burden_Project\\NEWWW_variant.pkl"
df = pd.read_pickle(file_path)

# Identify all tRNA-related columns
tRNA_cols = [col for col in df.columns if 'tRNA' in col]
charged_tRNAs = [col for col in tRNA_cols if col.startswith('charged-')]
uncharged_tRNAs = [col for col in tRNA_cols if not col.startswith('charged-')]

# Function to compute charging ratio per row
def compute_trna_charging_ratio_per_row(row):
    charged_arrs = [np.array(row[col]).flatten() for col in charged_tRNAs]
    uncharged_arrs = [np.array(row[col]).flatten() for col in uncharged_tRNAs]

    charged_sum = np.sum(charged_arrs, axis=0)
    uncharged_sum = np.sum(uncharged_arrs, axis=0)
    total_trna = charged_sum + uncharged_sum

    ratio = np.divide(
        charged_sum,
        total_trna,
        out=np.zeros_like(charged_sum, dtype=float),
        where=total_trna != 0
    )
    return ratio

# Add charging_ratio column
df['charging_ratio'] = df.apply(compute_trna_charging_ratio_per_row, axis=1)

# Save updated DataFrame
df.to_pickle("NEWWW_variant_with_charging.pkl")

print("Saved DataFrame with charging_ratio to NEWWW_variant_with_charging.pkl")


