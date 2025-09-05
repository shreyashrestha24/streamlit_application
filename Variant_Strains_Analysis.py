import streamlit as st
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from variant_utils import (
    load_config,
    load_data,
melt_timeseries_column,
melt_all_compounds,
compute_avg_lifespan,
compute_lifespan_by_repeat,
plot_avg_lifespan_per_variant,
calculate_compound_auc,
calculate_compound_metrics,
plot_metric_timeseries)

st.sidebar.header("Upload Variant Data ☺")

# data & config file upload
data_file = st.sidebar.file_uploader("Please Upload Pickle Data (.pkl)", type=["pkl"])
config_file = st.sidebar.file_uploader("Please Upload YAML Config", type=["yaml", "yml"])

# Load variant data
if 'df_variant' not in st.session_state and data_file is not None:
    try:
        df = load_data(data_file)
        st.session_state['df_variant'] = df
        st.success("Variant data loaded and saved to session_state.")
    except Exception as e:
        st.error(f"Error reading variant data file: {e}")

# Load variant config
if 'config_variant' not in st.session_state and config_file is not None:
    try:
        config = yaml.safe_load(config_file)
        st.session_state['config_variant'] = config
        st.success("Variant config loaded and saved to session_state.")
    except Exception as e:
        st.error(f"Failed to load variant config: {e}")

# only proceed if both are loaded
if 'df_variant' in st.session_state and 'config_variant' in st.session_state:
    df = st.session_state['df_variant']
    config = st.session_state['config_variant']

    # user selects compound in sidebar
    metabolites = config.get('metabolites', [])
    non_metabolites = config.get('non_metabolites', [])
    all_compounds = metabolites + non_metabolites

    st.sidebar.subheader("Select Compound for All Analyses")

    default_compound = st.session_state.get('selected_compound', None)

    selected_compound = st.sidebar.selectbox(
        "Select Compound",
        options=all_compounds,
        index=all_compounds.index(default_compound) if default_compound in all_compounds else 0
    )

    st.session_state['selected_compound'] = selected_compound

    # melt selected compounds
    long_df = melt_all_compounds(
        df,
        [selected_compound],
        metabolites=metabolites
    )

    print(long_df['compound'].unique())

    if selected_compound in long_df['compound'].unique():
        st.success(f"Successfully melted data for **{selected_compound}** "
                   f"({len(long_df)} rows).")
    else:
        st.error(f"No data found for {selected_compound}. Please recheck config or raw dataframe.")

    tabs = st.tabs([
            "Lineage Lifespan Analysis",
            "Summary Statistics (bar plots)",
            "Time-Series Plots"
        ])

    with tabs[0]:
        # Viability analysis
        lifespans_df = compute_lifespan_by_repeat(long_df)

        # Average lifespan per variant strain
        avg_lifespan = compute_avg_lifespan(lifespans_df)
        st.write("### Average Lifespan per Variant", avg_lifespan)

        # plot
        st.subheader("Average Lifespan per Variant")
        plot_avg_lifespan_per_variant(avg_lifespan)


        with tabs[1]:
            selected_compound = st.session_state['selected_compound']

            # Let user choose what metric to view
            options = ['min', 'mean', 'max', 'final', 'AUC']
            selected_options = st.multiselect(
                "Select metrics to view", options, default=[]
            )

            # average over repeats
            average_repeats = st.checkbox("Average over repeats", value=False)

            # Calculate metrics (per repeat)
            compound_metrics = calculate_compound_metrics(long_df, selected_compound, average_repeats=average_repeats)

            # Calculate AUC
            auc_df = calculate_compound_auc(long_df, selected_compound, average_repeats=average_repeats)

            # display df
            display_cols = ['strain_id']
            if not average_repeats:
                display_cols.append('repeat')

            display_df = compound_metrics[display_cols].copy()

            if 'min' in selected_options:
                display_df[f'{selected_compound}_min'] = compound_metrics['min_conc']
            if 'mean' in selected_options:
                display_df[f'{selected_compound}_mean'] = compound_metrics['mean_conc']
            if 'max' in selected_options:
                display_df[f'{selected_compound}_max'] = compound_metrics['max_conc']
            if 'final' in selected_options:
                display_df[f'{selected_compound}_final'] = compound_metrics['final_conc']
            if 'AUC' in selected_options:
                display_df = display_df.merge(auc_df, on=display_cols, how='left')

            display_df_app = display_df.drop(columns=['gene_expression_factor', 'translation_efficiency'],
                                             errors='ignore')

            st.write(f"Selected metrics for {selected_compound}", display_df_app)

            st.subheader("Bar Plots")

            # Converting repeat to string for proper grouping
            if 'repeat' in display_df_app.columns:
                display_df_app['repeat'] = display_df_app['repeat'].astype(str)

            # Metrics to plot
            plot_options = [col for col in display_df_app.columns if col not in ['strain_id', 'repeat']]
            selected_plot_metrics = st.multiselect("Select metrics to plot", plot_options, default=[])

            # one plot per selected metric
            for metric in selected_plot_metrics:
                fig, ax = plt.subplots(figsize=(10, 5))
                sns.barplot(
                    data=display_df_app,
                    x='strain_id',
                    y=metric,
                    hue='repeat' if 'repeat' in display_df_app.columns else None,
                    dodge=True
                )
                ax.set_title(f"{metric} for {selected_compound}")
                ax.set_xlabel("Strain ID")
                ax.set_ylabel(metric)
                plt.xticks(rotation=45)
                plt.tight_layout()
                st.pyplot(fig)

        with tabs[2]:
            selected_compound = st.session_state['selected_compound']

            # user can select more than one strain via multiselect
            strain_ids = long_df['strain_id'].unique()
            selected_strains = st.multiselect("Select Strain(s)", strain_ids, default=[strain_ids[0]])

            # checkbox to average over repeats
            average_repeats = st.checkbox("Average across repeats?", value=False)

            # Loop over each selected strain
            for strain_id in selected_strains:
                st.write(f"### Strain: {strain_id}")
                plot_metric_timeseries(
                    long_df,
                    strain_id=strain_id,
                    metric=selected_compound,
                    average_repeats=average_repeats
                )


