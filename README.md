# Streamlit application for analysing E. coli WCM data

This repository contains an applications created to analyse high-throughput time-series WCM simulations of engineered E. coli strains expressing a heterologous product. 

## Application Structure
This application is structured in 3 pages:
1. **Home page**: containing instructures and overview
2. **Wildtype Strain Analysis**: Lifespan, heatmaps, summary statistics, time-series plots
3. **Variant Strain Analysis**: Lifespan, summary statistics, time-series plots

##  Instructions
1. clone respository "gitclone https://github.com/shreyashrestha24/streamlit_application.git" "cd streamlit_application"
2. Install required Python packages (Python 3.12) "pip install streamlit, pandas, numpy, scipy, matplotlib, seaborn, pyyaml"
3. Run Streamlit app in terminal with command "streamlit run Home.py"
4. Upload pre-processed simulation data (pickle/.pkl format). Simulation pickle files should be filtered for relevant columns before uploading due to Streamlit's 200MB limit. Other metrics like tRNA charging ratio should be calculated too.

Upload the config.yaml file containing columns of interest for wildtype and variant analysis and specify which columns are metabolites (IND, ATP) and non-metabolites (growth, tRNA charging)

# 6. Enjoy your analysis :)

This application was developed in PyCharm v2025.1.2

