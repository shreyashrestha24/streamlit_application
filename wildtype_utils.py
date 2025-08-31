import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import trapz

def load_data():
    data_file = st.file_uploader("Upload Baseline data file (.pkl)", type=["pkl"])

    if data_file is not None:
        try:
            df = pd.read_pickle(data_file)
            st.success("Baseline Data loaded successfully.")
            return df
        except Exception as e:
            st.error(f"Failed to load data: {e}")
            return None
    else:
        st.warning("Please upload a data file.")
        return None


def melt_baseline_timeseries_columns(df, columns, time_column='time'):
    records = []

    for _, row in df.iterrows():
        times_arr = np.array(row[time_column]).flatten()
        ctm_arr = np.array(row['countsToMolar']).flatten()

        for column_name in columns:
            values_arr = np.array(row[column_name]).flatten()

            if len(times_arr) != len(values_arr):
                continue

            for t, v, ctm in zip(times_arr, values_arr, ctm_arr):
                if pd.isna(v):
                    continue

                # Use countsToMolar if meaningful for variable
                if column_name in ['countsToMolar', 'charging_ratio', 'growth']:
                    concentration = v
                else:
                    concentration = v * ctm

                records.append({
                    'simulation_identifier': row['simulation_identifier'],
                    'generation': row['generation'],
                    'time': t,
                    'variable': column_name,
                    'value': v,
                    'countsToMolar': ctm,
                    'concentration': concentration
                })

    return pd.DataFrame(records)

def plot_max_generation_distribution(long_df: pd.DataFrame, title: str = "Max Generation per Simulation Run"):
    """Plot bar chart of max generation reached per repeat"""
    if long_df.empty or 'simulation_identifier' not in long_df.columns or 'generation' not in long_df.columns:
        st.warning("long_df is empty or missing required columns.")
        return

    max_gens = long_df.groupby('simulation_identifier')['generation'].max().reset_index()

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(data=max_gens, x='simulation_identifier', y='generation', ax=ax)
    ax.set_xlabel("Simulation Identifier")
    ax.set_ylabel("Max Generation")
    ax.set_title(title)
    plt.xticks(rotation=90)
    plt.tight_layout()

    st.pyplot(fig)

def plot_generation_distribution(long_df: pd.DataFrame, title: str = "Distribution of Max Generations Across Simulations"):
    """Plot frequency of simulations by their max generation."""
    if long_df.empty or 'simulation_identifier' not in long_df.columns or 'generation' not in long_df.columns:
        st.warning("long_df is empty or missing required columns.")
        return

    # Get max generation per sim
    max_gens = long_df.groupby('simulation_identifier')['generation'].max()

    # Count how many sims hit each max generation
    gen_counts = max_gens.value_counts().sort_index().reset_index()
    gen_counts.columns = ['Max Generation', 'Simulation Count']

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.barplot(data=gen_counts, x='Max Generation', y='Simulation Count', color='skyblue', ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Max Generation")
    ax.set_ylabel("Number of Simulations")
    plt.tight_layout()

    st.pyplot(fig)

def mean_max_generation(long_df, sim_id_col='simulation_identifier', generation_col='generation'):
    # Check if required columns exist
    if long_df.empty or sim_id_col not in long_df.columns or generation_col not in long_df.columns:
        raise ValueError(f"DataFrame missing required columns: '{sim_id_col}' and/or '{generation_col}'")

    # Get max generation per simulation
    max_gens = long_df.groupby(sim_id_col)[generation_col].max()

    # Calculate mean of max generations
    mean_max_gen = max_gens.mean()

    return mean_max_gen

def compute_full_auc(long_df, variable_name):
    """ calculates AUC for the specified variable across simulations"""
    # Filter for variable
    var_df = long_df[long_df['variable'] == variable_name]

    auc_results = []
    for sim_id, group in var_df.groupby('simulation_identifier'):
        # sort by time for integration
        group_sorted = group.sort_values('time')

        times = group_sorted['time'].values
        concentrations = group_sorted['concentration'].values

        if len(times) > 1:
            auc = np.trapz(y=concentrations, x=times)
        else:
            auc = np.nan

        auc_results.append({'simulation_identifier': sim_id, 'AUC': auc})

    return pd.DataFrame(auc_results)


def compute_mean_auc(long_df, metabolite):
    """ Computes the average AUC for a specific metabolite across all repats"""
    aucs = []

    for sim_id, group in long_df[long_df['variable'] == metabolite].groupby('simulation_identifier'):
        group = group.sort_values('time')
        x = group['time'].values
        y = group['concentration'].values

        if len(x) < 2:
            continue

        sim_auc = trapz(y, x)
        aucs.append(sim_auc)

    mean_auc = np.mean(aucs) if aucs else 0
    return mean_auc, aucs

def average_compound_concentration(df, compound):
    # compute average concentration for any compound
    filtered = df[df['variable'] == compound]
    avg_df = filtered.groupby('simulation_identifier')['concentration'].mean().reset_index(name='mean_concentration')
    return avg_df

def overall_mean_compound_concentration(df, compound):
    filtered = df[df['variable'] == compound]
    return filtered['concentration'].mean()

#new function
def summary_compound_concentrations(df, compound):
    """ computes summary concentration statistics (min, max, mean) per repeat for a given compound."""
    # Filter for the chosen compound
    filtered = df[df['variable'] == compound]

    # per repeat summary
    per_sim_df = (
        filtered.groupby('simulation_identifier')['concentration']
        .agg(['mean', 'min', 'max'])
        .reset_index()
        .rename(columns={'mean': 'mean_concentration',
                         'min': 'min_concentration',
                         'max': 'max_concentration'})
    )

    # Overall summary
    overall_stats = {
        'compound': compound,
        'overall_mean_concentration': filtered['concentration'].mean(),
        'overall_min_concentration': filtered['concentration'].min(),
        'overall_max_concentration': filtered['concentration'].max()
    }

    return per_sim_df, overall_stats


def average_concentration_time_series_binned(df, variable, bin_size=60, group_by='time', separate_repeats=False):
    """ Bin time points and average concentrations for a single compound"""
    df = df[df['variable'] == variable].copy()
    if df.empty:
        return pd.DataFrame()

    if group_by == 'time':
        bins = np.arange(df['time'].min(), df['time'].max() + bin_size, bin_size)
        df['time_bin'] = pd.cut(df['time'], bins=bins, right=False)
        if separate_repeats:
            grouped = (
                df.groupby(['simulation_identifier', 'time_bin'])['concentration']
                .mean()
                .reset_index(name='average_concentration')
            )
        else:
            grouped = (
                df.groupby('time_bin')['concentration']
                .mean()
                .reset_index(name='average_concentration')
            )
        grouped['time_mid'] = grouped['time_bin'].apply(lambda x: x.left + bin_size / 2)
        if separate_repeats:
            return grouped[['time_mid', 'average_concentration', 'simulation_identifier']]
        else:
            return grouped[['time_mid', 'average_concentration']]

    elif group_by == 'generation':
        if separate_repeats:
            grouped = df.groupby(['simulation_identifier', 'generation'])['concentration'].mean().reset_index(name='average_concentration')
            return grouped[['generation', 'average_concentration', 'simulation_identifier']]
        else:
            grouped = df.groupby('generation')['concentration'].mean().reset_index(name='average_concentration')
            return grouped[['generation', 'average_concentration']]
    else:
        raise ValueError("group_by must be 'time' or 'generation'")


def plot_baseline_heatmap(long_df, compounds):
    """ Plots heatmap of AUCs per simulation and compound using existing compute_full_auc"""
    auc_dfs = []
    for compound in compounds:
        auc_df = compute_full_auc(long_df, compound)
        auc_df = auc_df.rename(columns={'AUC': compound})
        auc_dfs.append(auc_df.set_index('simulation_identifier'))

    if not auc_dfs:
        st.warning("No compounds selected or no data available.")
        return

    combined_auc_df = pd.concat(auc_dfs, axis=1)

    combined_auc_df = combined_auc_df.fillna(0)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(combined_auc_df, cmap='viridis', annot=True, fmt=".2f", cbar_kws={'label': 'AUC'}, ax=ax)
    ax.set_title("Baseline AUC Heatmap per Simulation and Compound")
    ax.set_xlabel("Compound")
    ax.set_ylabel("Simulation Identifier")
    plt.xticks(rotation=45)
    plt.tight_layout()

    st.pyplot(fig)

def plot_avg_concentration_heatmap(df, compound, figsize=(12, 8), cmap="viridis"):
    """ Plots a heatmap of average concentrations per simulation across generations for a given compound."""
    comp_df = df[df['variable'] == compound]

    if comp_df.empty:
        return None

    heatmap_data = (
        comp_df
        .groupby(['simulation_identifier', 'generation'])['concentration']
        .mean()
        .unstack(level='generation')
    )

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(heatmap_data, cmap=cmap, ax=ax)
    ax.set_title(f"{compound} by Simulation and Generation")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Simulation ID")
    return fig
