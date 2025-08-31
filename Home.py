import os
import json
import yaml
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

from wildtype_utils import (
    melt_baseline_timeseries_columns,
    plot_max_generation_distribution,
    plot_generation_distribution,
    average_concentration_time_series_binned,
    plot_baseline_heatmap,
    plot_avg_concentration_heatmap,
    mean_max_generation,
    summary_compound_concentrations,
)
from variant_utils import (
    calculate_compound_metrics,
    calculate_compound_auc,
    plot_metric_timeseries,
)

st.set_page_config(page_title="E. coli Whole Cell Data Exploration", layout="centered")

st.title("E. coli Whole-Cell Model Data Exploration~~~")
st.markdown("""
Welcome to the **E. coli Whole Cell Modelling Data Analysis App**, created to analyse high-throughput time-series WCM simulations of *engineered E. coli* strains expressing a heterologous product

This app allows you to explore the metabolic burden during heterologous expression in E. coli comparing **wildtype** and **variant** simulations.

---

### Pages:

- **Baseline Analysis**  
  Explore wildtype simulation runs, including viability analysis, time-series visualisations, and summary statistics inc. AUC, and mean, max, and min concentrations 

- **Variant Analysis**  
 Explore variant simulation runs, including viability analysis, summary metrics, and time-series visualisations
  

---

### How to Use:
1. Upload your `.pkl` data file and `config.yaml` file on each page.
2. Explore simulations using plots and tables.
3. Use filters to narrow by simulation, metabolite, or repeats

Enjoy your Analysis! ☺
---

""")


