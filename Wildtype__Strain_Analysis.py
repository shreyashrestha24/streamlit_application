import json
import streamlit as st
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from wildtype_utils import (
    melt_baseline_timeseries_columns,
    plot_max_generation_distribution,
    plot_generation_distribution,
    average_compound_concentration,
    overall_mean_compound_concentration,
    average_concentration_time_series_binned,
    plot_baseline_heatmap,
    plot_avg_concentration_heatmap,
    mean_max_generation,
    summary_compound_concentrations
)

st.title("Wildtype Data Analysis")

tabs = st.tabs([
    "Lifespan Analysis",
    "Heatmap Analysis",
    "Wildtype Summary Statistics",
    "Time-Series Plots"
])

st.sidebar.header("Upload Baseline Data ☺")

data_file = st.sidebar.file_uploader("Upload Pickle Data (.pkl)", type=["pkl"])
config_file = st.sidebar.file_uploader("Upload YAML Config", type=["yaml", "yml"])

# load data file into session_state only if not already loaded
if 'df_baseline' not in st.session_state and data_file is not None:
    try:
        df = pd.read_pickle(data_file)
        st.session_state['df_baseline'] = df
        st.success(" Data loaded and saved to session_state.")
    except Exception as e:
        st.error(f"Error reading data file: {e}")

# load config file into session_state only if not already loaded
if 'config_baseline' not in st.session_state and config_file is not None:
    try:
        config = yaml.safe_load(config_file)
        st.session_state['config_baseline'] = config
    except Exception as e:
        st.error(f"Failed to load config: {e}")

# proceed only if both data and config are loaded
if 'df_baseline' in st.session_state and 'config_baseline' in st.session_state:
    df = st.session_state['df_baseline']
    config = st.session_state['config_baseline']

    st.subheader("Columns present in Baseline Data")

    if st.checkbox("Show DataFrame Columns"):
        st.write(df.columns.tolist())

    # get columns of interest from config
    all_columns = config.get('wildtype_columns_of_interest', [])

    # using only these columns
    df = df[all_columns]

    # define columns that are NOT compounds
    non_compound_cols = {'simulation_identifier', 'time', 'countsToMolar', 'generation'}

    # filter to keep only compound columns
    compound_columns = [col for col in all_columns if col not in non_compound_cols]

    # melting data
    long_df = melt_baseline_timeseries_columns(df, compound_columns)


    with tabs[0]:
        st.subheader("Lineage Survival of Baseline Data")
        # max generation summary
        max_gens = long_df.groupby('simulation_identifier')['generation'].max().reset_index()
        max_gens.rename(columns={'simulation_identifier': 'Simulation ID', 'generation': 'Max Generation'}, inplace=True)

        st.subheader("Simulation Lifespan (Max Generation) per Run")
        st.dataframe(max_gens)

        st.subheader("Frequency of Max Generations")
        plot_generation_distribution(long_df)

        st.subheader("Max Generation Distribution")
        plot_max_generation_distribution(long_df)


        mean_gen = mean_max_generation(long_df)
        st.write(f"Average max generation: {mean_gen:.2f}")

    with tabs[1]:
        st.subheader("Heatmap Analysis")

        compound = st.selectbox("Select compound:", sorted(long_df['variable'].unique()))
        fig = plot_avg_concentration_heatmap(long_df, compound)

        if fig:
            st.pyplot(fig)
        else:
            st.warning("No data available for the selected compound.")

    with tabs[2]:
        st.subheader("Summary Statistics: Wildtype Baseline Values")

        if not compound_columns:
            st.warning("No compound columns found in config")
        else:
            # let user pick the compound
            compound = st.selectbox("Select compound for analysis", compound_columns, index=0)

            # concentration summaries
            per_sim_df, overall_stats = summary_compound_concentrations(long_df, compound)

            # select metrics to display
            metric_choices = st.multiselect(
                "Select summary metrics to display",
                ["Mean", "Min", "Max"],
                default=["Mean", "Min", "Max"]
            )

            # display only selected metrics
            if "Mean" in metric_choices:
                st.metric(f"Overall Mean {compound} Concentration",
                          f"{overall_stats['overall_mean_concentration']:.4f} mM")

            if "Min" in metric_choices:
                st.metric(f"Overall Min {compound} Concentration",
                          f"{overall_stats['overall_min_concentration']:.4f} mM")

            if "Max" in metric_choices:
                st.metric(f"Overall Max {compound} Concentration",
                          f"{overall_stats['overall_max_concentration']:.4f} mM")

            # JSON export for selected stats
            export_stats = {}
            if "Mean" in metric_choices:
                export_stats["overall_mean_concentration"] = round(overall_stats['overall_mean_concentration'], 4)
            if "Min" in metric_choices:
                export_stats["overall_min_concentration"] = round(overall_stats['overall_min_concentration'], 4)
            if "Max" in metric_choices:
                export_stats["overall_max_concentration"] = round(overall_stats['overall_max_concentration'], 4)

            json_data = json.dumps(export_stats, indent=2)

            st.download_button(
                label=f"Download {compound} Concentration Summary (JSON)",
                data=json_data,
                file_name=f"{compound}_concentration_summary.json",
                mime='application/json'
            )

            # per simulation checkbox
            show_table = st.checkbox(f"Show Per-Simulation Concentration Table for {compound}")
            if show_table:
                column_map = {
                    "Mean": "mean_concentration",
                    "Min": "min_concentration",
                    "Max": "max_concentration"
                }
                selected_cols = ["simulation_identifier"] + [
                    column_map[m] for m in metric_choices if m in column_map
                ]
                st.dataframe(per_sim_df[selected_cols])


    with tabs[3]:
        compound = st.selectbox("Select compound to visualize", sorted(long_df['variable'].unique()))

        group_by = st.radio("Group data by:", options=['time', 'generation'], index=0)

        # Option for plotting repeats separately
        plot_repeats_separately = st.checkbox("Plot repeats separately", value=False)

        # Bin size slider
        bin_size = st.slider("Bin size", min_value=10, max_value=300, step=10, value=60)

        if plot_repeats_separately:
            # Keep repeats separate
            binned_df = average_concentration_time_series_binned(long_df, compound, bin_size, group_by=group_by,
                                                                 separate_repeats=True)
        else:
            # Average over repeats
            binned_df = average_concentration_time_series_binned(long_df, compound, bin_size, group_by=group_by,
                                                                 separate_repeats=False)

        # Determine x-axis column and label
        x_col = 'time_mid' if group_by == 'time' else 'generation'
        x_label = "Time (s)" if group_by == 'time' else "Generation"

        st.subheader(f"{compound} Over {x_label} (Binned every {bin_size} {x_label.lower()})")

        fig, ax = plt.subplots(figsize=(10, 5))

        if plot_repeats_separately:
            # Plot each repeat with a different color
            sns.lineplot(data=binned_df, x=x_col, y='average_concentration', hue='simulation_identifier', ax=ax, palette='tab10')
        else:
            # Plot averaged values
            sns.lineplot(data=binned_df, x=x_col, y='average_concentration', ax=ax, color='blue')

        ax.set_xlabel(x_label)
        ax.set_ylabel(f"{compound}")
        ax.grid(True)

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
        plt.tight_layout()

        st.pyplot(fig)

else:
    st.warning("Please upload both a data file and a config file to proceed.")
    st.stop()
