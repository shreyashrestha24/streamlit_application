import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st

# Loading required files into App
def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_data(data_path):
    return pd.read_pickle(data_path)

#Melting time-series data into long-format
def melt_timeseries_column(df, column_name, time_column='time', metabolites=None):
    """
    Melt time series column into long format
    """
    records = []
    metabolites = metabolites or []

    for _, row in df.iterrows():

        times_arr = np.array(row[time_column]).flatten() if np.iterable(row[time_column]) else np.array([row[time_column]])
        values_arr = np.array(row[column_name]).flatten() if np.iterable(row[column_name]) else np.array([row[column_name]])

        ctm = row.get('countsToMolar', 1)
        ctm_arr = np.array(ctm).flatten() if np.iterable(ctm) else np.full(len(values_arr), ctm)

        # Skip rows where lengths don’t match
        if len(times_arr) != len(values_arr) or len(values_arr) != len(ctm_arr):
            continue

        for t, v, c in zip(times_arr, values_arr, ctm_arr):
            if v is None or (isinstance(v, (int, float, np.number)) and np.isnan(v)):
                continue

            concentration = v * c if column_name in metabolites else v
            unit = "mM" if column_name in metabolites else "unitless"

            records.append({
                'strain_id': row['strain_id'],
                'generation': row['generation'],
                'time': t,
                'compound': column_name,
                'value': v,
                'countsToMolar': c,
                'concentration': concentration,
                'unit': unit,
                'gene_expression_factor': row.get('gene_expression_factor', np.nan),
                'translation_efficiency': row.get('translation_efficiency', np.nan),
                'repeat': row['repeat']
            })

    return pd.DataFrame(records)

def melt_all_compounds(df, columns, metabolites=None):
    dfs = [melt_timeseries_column(df, col, metabolites=metabolites) for col in columns]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# Viability Analysis
def compute_lifespan_by_repeat(df):
    """computes the max generation per strain, repeat, GE, and TE"""
    lifespans = (
        df.groupby(['strain_id', 'repeat', 'gene_expression_factor', 'translation_efficiency'])
        ['generation']
        .max()
        .reset_index()
        .rename(columns={'generation': 'max_generation'})
    )
    return lifespans


def compute_avg_lifespan(lifespans_df):
    """computes average max generation per variant (strain_id) across repeats"""
    avg_lifespan = (
        lifespans_df.groupby('strain_id', as_index=False)
        .agg(
            avg_max_generation=('max_generation', 'mean'),
            gene_expression_factor=('gene_expression_factor', 'first'),
            translation_efficiency=('translation_efficiency', 'first')
        )
    )
    return avg_lifespan


def plot_avg_lifespan_per_variant(avg_lifespan):
    """Barplot for average max generation per variant"""
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.barplot(ax=ax, data=avg_lifespan, x='strain_id', y='avg_max_generation', palette='viridis')
    ax.set_title('Average Max Generation per Variant')
    ax.set_xlabel('Variant (strain_id)')
    ax.set_ylabel('Average Max Generation')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.tight_layout()
    st.pyplot(fig)


# Summary Statistics

def calculate_compound_metrics(long_df, compound, average_repeats=False):
    """ calculates key metrics (max, mean, min, final concentration) for a selected compound"""
    # Filter for selected compound
    comp_df = long_df[long_df['compound'] == compound].copy()

    # Group by strain and repeat
    grouped = comp_df.groupby(['strain_id', 'repeat']) if not average_repeats else comp_df.groupby('strain_id')

    metrics_records = []
    for name, group in grouped:
        max_conc = group['concentration'].max()
        mean_conc = group['concentration'].mean()
        min_conc = group['concentration'].min()
        # final cooncentration at last time point
        final_time = group['time'].max()
        final_conc = group[group['time'] == final_time]['concentration'].mean()

        record = {
            'max_conc': max_conc,
            'mean_conc': mean_conc,
            'min_conc': min_conc,
            'final_conc': final_conc
        }

        # grouping info
        if not average_repeats:
            record['strain_id'], record['repeat'] = name
        else:
            record['strain_id'] = name

        record['gene_expression_factor'] = group['gene_expression_factor'].iloc[0] if 'gene_expression_factor' in group else np.nan
        record['translation_efficiency'] = group['translation_efficiency'].iloc[0] if 'translation_efficiency' in group else np.nan

        metrics_records.append(record)

    return pd.DataFrame(metrics_records)

def calculate_compound_auc(long_df, selected_compound, average_repeats=False):
    """calculates area under the curve (AUC) for the selected compound"""

    # Filter for the selected compound
    compound_data = long_df[long_df['compound'] == selected_compound].copy()

    # Ensure sorted by strain, repeat, and time
    compound_data = compound_data.sort_values(by=['strain_id', 'repeat', 'time'])

    auc_records = []
    grouped = compound_data.groupby(['strain_id', 'repeat'])

    for (strain_id, repeat), group in grouped:
        if len(group) >= 2:
            auc = np.trapz(group['concentration'], group['time'])
        else:
            auc = np.nan

        auc_records.append({
            'strain_id': strain_id,
            'repeat': repeat,
            f'{selected_compound}_auc': auc,
            'gene_expression_factor': group['gene_expression_factor'].iloc[0]
                                      if 'gene_expression_factor' in group else np.nan,
            'translation_efficiency': group['translation_efficiency'].iloc[0]
                                      if 'translation_efficiency' in group else np.nan
        })

    auc_df = pd.DataFrame(auc_records)

    if average_repeats:
        # Aggregate mean over repeats
        auc_df = (
            auc_df.groupby('strain_id', as_index=False)
                  .agg({
                      f'{selected_compound}_auc': 'mean',
                      'gene_expression_factor': 'first',
                      'translation_efficiency': 'first'
                  })
        )

    return auc_df


# Time-Series plots

def plot_metric_timeseries(long_df, strain_id, metric, average_repeats=False, bin_size=60):
    """Plot a metric over time for a given strain"""

    # Filter for strain and metric
    if metric in long_df['compound'].unique():
        df_plot = long_df[(long_df['strain_id'] == strain_id) & (long_df['compound'] == metric)].copy()
        y_col = 'concentration'
    else:
        df_plot = long_df[long_df['strain_id'] == strain_id].copy()
        if metric not in df_plot.columns:
            st.warning(f"Metric '{metric}' not found in dataframe.")
            st.write("Available columns:", df_plot.columns.tolist())
            return
        y_col = metric

    if df_plot.empty:
        st.warning(f"No data found for {strain_id} and metric '{metric}'")
        return

    df_plot = df_plot.sort_values('time')
    df_plot['y_value'] = df_plot[y_col]

    fig, ax = plt.subplots(figsize=(10, 5))

    if average_repeats:
        # Bin time and average
        df_plot['time_bin'] = (df_plot['time'] // bin_size) * bin_size
        grouped = df_plot.groupby('time_bin')['y_value'].mean().reset_index()
        sns.lineplot(data=grouped, x='time_bin', y='y_value', ax=ax, color='blue')
    else:
        # Plot each repeat separately
        if 'repeat' in df_plot.columns:
            df_plot['repeat'] = df_plot['repeat'].astype(str)
            sns.lineplot(data=df_plot, x='time', y='y_value', hue='repeat', palette='tab10', ax=ax)
        else:
            sns.lineplot(data=df_plot, x='time', y='y_value', ax=ax, color='blue')

    ax.set_title(f"{metric} over Time for {strain_id}")
    ax.set_xlabel("Time")
    ax.set_ylabel(metric)
    ax.grid(True)
    plt.xticks(rotation=45)
    plt.tight_layout()

    st.pyplot(fig)






