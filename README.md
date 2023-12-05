# The impact of Covid-19 vaccination in Aotearoa New Zealand: a modelling study



# 1. Overview

This repository contains the code for the article 'The impact of Covid-19 vaccination in Aotearoa New Zealand: a modelling study', which shows the effect of different vaccination strategies (counterfactuals to what occurred in reality) on the evolution of the pandemic.

The code includes a simple approximate Bayesian computation (ABC) method to fit the model to epidemiological data, and account for the effect of uncertain parameters. The code also inclues a scenario simulation feature allowing to change several model parameters and simulate their effect on daily infections, reported cases, hospital admissions, hospital occupancy, and fatalities.

An earlier version of the code on this repo (see [here](https://gitlab.com/tpmcovid/ode-model) was peer reviewed by Ning Hua at Precision Driven Health. Documentation of the code following this review can be found in the `doc` folder on this repo.  

Only code and other small files are kept in this repo. The data and outputs produced can be kept locally, and are not automatically pushed to the repo. The raw data cannot be shared publicly due to confidentiality issues. 


# 2. Prerequisites and installation
This project was coded using version 2021b of Matlab. Some of the features might not work well if using previous versions of the software.

To install the COVID-19 ODE model from a terminal or command line:
1) Open the terminal and navigate to the desired destination folder
2) Type `git clone https://github.com/SamikDatta/covid19_vaccination.git` and press enter
3) All project folders and files should now be in your chosen repository

To install the COVID-19 ODE model from the Matlab software
1) Open Matlab
2) On the **Home** tab, click **New > Project > From Git**.
3) Enter the HTTPS repository path `https://github.com/SamikDatta/covid19_vaccination.git` into the **Repository path** field
4) In the **Sandbox** field, select the working folder where you want to put the retrieved files for your new project
5) Click **Retrieve**
6) All folders and files should now be in your chosen repository

Here is a list and description of the folders that will get cloned in the chosen repository together with the script in the top-level directory and this `README.md`:
* **functions**: this folder contains all functions used to run the main script
* **data**: this is where the datafiles needed to run the script are kept, see the table below for a full list
* **results**: this folder will only contain a README file at installation, but will be filled with simulation results and graph after running the main script
* **postProcessing**: this folder contains the functions needed to process the outputs from running the ODE model, as described further in the ReadMe in that subfolder
* **doc**: code peer review documentation on an earlier version of this repo.
* **latex**: latex source code generated for the Tables in the article.


To run the COVID-19 ODE model, you will also have to make sure the following data files are in the indicated folders:

|Filename|Description|
|--------------------|-----|
|*epidata_by_age_and_vax_[DD-MMM-YYYY].mat*| Number of reported Covid-19 cases, Covid-19 hospital admissions and Covid-19 deaths on each date by 10-year age group, vaccination status  (not available publicly). |
|*vaccine_data_national_[YYYY-MM-DD].mat*| Cumulative total number of 1st, 2nd, 3rd and 4th or subsequent vaccine doses given by data and 5-year age band (not available publicly). |
|*vaccine_data_Maori_[YYYY-MM-DD].mat*| Cumulative number of 1st, 2nd, 3rd and 4th or subsequent vaccine doses given by data and 5-year age band for people whose prioritised ethnicity was Māori (not available publicly). |
|*vaccine_data_EuropeanOther_[YYYY-MM-DD].mat*| Cumulative number of 1st, 2nd, 3rd and 4th or subsequent vaccine doses given by data and 5-year age band for people whose prioritised ethnicity was European or other (not available publicly). |
|*reshaped_b2_projections_final_[YYYY-MM-DD].csv*| This was constructed by combining historical vaccine uptake data with future uptake projections provided by the Ministry of Health. As these projections have now been superceded by more recent uptake data, this file is now obsolete and not needed. |
|*therapeutics_by_age_[DD-MMM-YYYY].mat*| MOH data for the proportion of cases with an antiviral prescription in each group over time (not available publicly). |
|*covid-cases-in-hospital-counts-location-16-Aug-2023.xlsx*| Data from file *covid-cases-in-hospital-counts-location.xlsx* file in the *nz-covid-data/cases/* folder of the [MOH Github repo](https://github.com/minhealthnz/nz-covid-data/tree/main/cases).|
|*border_incidence.xlsx*|*data*| Data collected by M. Plank and stored in [this GitHub page](https://github.com/michaelplanknz/modelling-ba5-in-nz). |
|*nzcontmatrix.xlsx*| This is taken from the paper 'Projecting social contact matrices in 152 countries using contact surveys and demographic data' by Prem et al. (specifically the overall age-and-location-specific contact matrix for New Zealand), and is available [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697). |
|*popsize_national.xlsx*| National population size estimates in 5-year age bands (HSU population 2022). |
|*popproj_national2018-21.xlsx*| National population size estimates in 5-year age bands (StatsNZ projected population for 2021 - see [here](https://www.tepou.co.nz/resources/dhb-population-profiles-2021-2031-pdf) ). |
|*HSU_by_age_eth.csv*| Population size estmiates by prioritised ethnicity in 5-year age bands (HSU population 2022) |
|*popproj2018-21.csv*| Population size estmiates by prioritised ethnicity in 5-year age bands (StatsNZ projected population for 2021 - see [here](https://www.tepou.co.nz/resources/dhb-population-profiles-2021-2031-pdf) ). |
|*actual_deaths_by_age_and_sex.csv*| Number of recorded Covid-19 deaths by sex and 1-year age band for 1 Jan 2022 - 30 Jun 2023 (not available publicly). |
|*maori_outcomes_by_age.csv*| Number of recorded Māori Covid-19 hospitalisations and deaths by 5-year age band for 1 Jan 2022 - 30 Jun 2023  (not available publicly). |
|*nz-complete_cohort-life-tables-1876-2021.csv*| Cohort life tables used for calculating YLL [published by StatsNZ](https://www.stats.govt.nz/information-releases/new-zealand-cohort-life-tables-march-2023-update/). |


# 3. How to run the ODE model

For all code, set the Matlab current directory to the Gitlab repo root folder.

First, the model is run in two stages. The script `fitModel.m` runs the ABC fitting algorithm for parameter estimation and saves the accepted parameter combinations and associated model output in the output file:
* *results_Uni_Filtered100_DDMMMFit_DDMMMRun.mat*. 

Second, the script `runScenarios.m` reads in the results of `fitModel.m` and runs a number of alternative scenarios for each of the accepted parameter combinations. This produces the following output files, which are described in more detail in the following sections:
* *results_Uni_Filtered100_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_FilteredBest_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.png*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit_ageSplit.png*

## Model initialisation: 
Both `fitModel.m` and `runScenarios.m` begin with a common initialisation section that specifies the names of the parameters to be fitted bythe ABC algorithm (`getParsToFit.m`) and the path and file names for input data (`getDataFileNames.m`), and inputs all data and defines fixed parameter values. Before running these scripts, **ensure the correct filenames have been defined** in `getDataFileNames.m`, and that the corresponding files are all in the *data* folder. For a complete list of the required files and where to get them, please see the "Prerequisites" section of this README.


## Model fitting
`fitModel.m` contains the simple ABC method used to narrow down the prior distribution of fitted parameters to the posterior distribution of values which resulted in the 1% best fitting model results (see Supplementary section S9 and Tables S1-S2 of the article for details). 

The fitting procedure can take several hours depending on the value chosen for `nSamples`. If using all 16 default fitted parameters, we recommend having at least `nSamples` = 15,000 to ensure enough combinations of parameter values get tested. Note that using 15,000 will take approximately 10 hours to run, but this estimate might vary depending on the available processing power.

The following outputs will be produced and saved in the `results` folder:
* *results_Uni_Filtered100_DDMMMFit_DDMMMRun.mat* : a Matlab structure containing  a feeld for each of the model output. Each field contains a matrix with 16 columns corresponding to 16 age groups and a line for each day simulated. For a description of the outputs, please see the preamble of the `functions/extractEpiVarsCompact.m` function. The output file is labeled with the datestamp of the line data used for the fit, and the date in which the fit was run.
* *violinPlots_DDMMMFit_DDMMMRun.png* : a figure representing the marginal posterior distributions of the fitted parameters on a 0-1 scale. This figure can be used to gauge which parameter values gave the best 1% fit to the data.
* *CHR_CFR.png* : a figure comparing the age-specific case hospitalisation ratio and case fatality ratio according to model versus data.

## Running alternative scenarios
`runScenarios.m` reads in the fitted posterior distribution for parameters saved by a previous run of `fitModel.m` (see above). Note the variable `fitLbl` in `runScenarios.m` specifies the datestamp on the file containing the output of the model fitting routine (this is the date that `fitModel.m` was run).  **Caution: if any of the underlying model parameters are changed, it is necessary to re-run the fitting procedure.**

Scenarios are run from the pre-loaded posterior parameter sets and using the different scenarios for vaccination rates, which are defined in the last section of `getBasePar.m`.
Once the scenario parameter values have been defined in `getBasePar.m`, the script will run each scenario using all the accepted samples (1% best) from the posterior distribution of parameter values that have been produced earlier.
The following output file will then be saved in the `results` folder:
* *results_Uni_Filtered100_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for each of the posterior parameter sets, labeled by scenario name and datestamp of the data used for the model fitting. Note that this is an intermediary result file, which will be processed in the next section.

The scenario results are then processed to produce the 95% bands (the 95% of accepted model simulations with the smallest value of the distance function) data and best-fit (the single model simulation the smallest value of the distance function) data. 
The following output files will be produced and saved in the `results` folder:
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for the 95% best fitting parameter sets out of the posterior parameter sets.
* *results_Uni_FilteredBest_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for the single best fitting parameter set.

Finally, plots are created for the aggregated and age-split model results in each scenario.
This will generate two figures in the `results` folder in .png format:
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.png* : a tiled plot of aggregated model outputs and data for daily infections, daily reported cases, daily hospital admissions, hospital occupancy, and daily deaths
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit_ageSplit.png* : a tiled plot of age-split model outputs and data for daily reported cases, daily hospital admissions, daily deaths, ratio of hospital admissions / cases, case ascertainment rate, and ratio of cumulative infections per total population


## Sensitivity analysis on population size data

If the `popSizeFname` field in `getDataFileNames.m` is changed from "popsize_national.xlsx" to "popproj_national2018-21.xlsx" (an alternative population projection dataset from Stats NZ -- see article for details), the sensitivity run will be performed. n this caes, relevant intermediary and output filenames are appended with `_sensitivity`. 

# 4. Post processing

The top-level script `writeModelOutput.m` reads in the results of the scenario models (output by `runScenarios.m`) and generates a single .csv file (`results/model_output.csv`) containing the the number of infections, number of first infections, number of admissions and number of deaths in each of the 16 age groups, for each of the 150 accepted parameter combinations, and for each of the 9 model scenarios. All values are aggregated over the relevant time period. 

The script `scenario_analysis.m` reads in this .csv file and calculates the median and 95% CrI values reported in Tables 1 and S6 of the article (latex source code for the Tables is generated and saved in `latex/results_table.tex`). This includes the calculation of YLL from cohort life tables, and the calculation of differences (Delta) between scenarios for the same parameter combination (see article Methods for details).   


# 5. Other 

Other graphs shown in the article may be produced by running the following scripts that are in the top-level directory of the repo:
* plotAntivirals.m -- plot of the  time series for the proportion of cases with an antiviral prescriptipon in each age group.
* plotREIt_prior.m -- plot of samples from the prior distribution of the time-varying reproduction number excluding immunity.


# 6. Acknowledgements 
The authors acknowledge the role of the New Zealand Ministry of Health, StatsNZ, and the Institute of Environmental Science and Research in supplying data in support of this work, and the Covid-19 Modelling Government Steering for helping to design modelling questions and to interpret model outputs. The authors are grateful to Ning Hua and Rachel Owens at Precision Driven Health for code peer review, and to Nigel French, Emily Harvey, Markus Luczak-Roesch, Dion O'Neale, Matt Parry and Patricia Priest for feedback on the model. The authors acknowledge the contributions of Rachelle Binny, Shaun Hendy, Kannan Ridings, Nicholas Steyn and Leighton Watson to a previous model from which this model was developed.

