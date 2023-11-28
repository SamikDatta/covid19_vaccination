# The impact of Covid-19 vaccination in Aotearoa New Zealand: a modelling study



# Overview

This repository contains the code for the article 'The impact of Covid-19 vaccination in Aotearoa New Zealand: a modelling study', which shows the effect of different vaccination strategies (counterfactuals to what occurred in reality) on the evolution of the pandemic.

The code includes a simple approximate Bayesian approximation (ABC) method to fit the model to epidemiological data, and account for the effect of uncertain parameters. The code also inclues a scenario simulation feature allowing to change several model parameters and simulate their effect on daily infections, reported cases, hospital admissions, hospital occupancy, and fatalities.

An earlier version of the code on this repo (see [here](https://gitlab.com/tpmcovid/ode-model) was peer reviewed by Ning Hua at Precision Driven Health. Documentation of the code following this review can be found in the `doc` folder on this repo.  

Only code and other small files are kept in this repo. The data and outputs produced can be kept locally, and are not automatically pushed to the repo. The raw data cannot be shared publicly due to confidentiality issues. 


# 2. Prerequisites and installation
This project was coded under version 2021b of Matlab. Some of the features might not work well if using previous versions of the software.

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

Here is a list and description of the folders that will get cloned in the chosen repository together with the main file `naiveABCSEIR.m` and this `README.md`:
* **functions**: this folder contains all functions used to run the main script
* **data**: this is where the datafiles needed to run the script are kept, see the table below for a full list
* **results**: this folder will only contain a README file at installation, but will be filled with simulation results and graph after running the main script
* **postProcessing**: this folder contains the functions needed to process the outputs from running the ODE model, as described further in the ReadMe in that subfolder
* **processMOHdata**: this folder contains the functions needed to process the line data on COVID-related cases, hospitalisations and deaths, which is used for the model fitting procedure and plots, as described further in Section 3 of this README
* **processVaxData**: this folder contains the functions needed to process and check the vaccination data, as described in Section 4 of this README
* **resources/project**


To run the COVID-19 ODE model, you will also have to make sure the following data files are in the indicated folders:

|Filename|Folder|How to get|
|--------------------|---|-----|
|*TPM_comm_cases_info_[YYYY-MM-DD].csv*|*processMOHdata*| This is the unit report data, and is sourced directly from the Ministry of Health, and is not publicly available. |
|*TPM_vaccine_[YYYY-MM-DD].csv*|*processVaxData*| This is the vaccination data, and is sourced directly from the Ministry of Health, and is not publicly available. |
|*epidata_by_age_and_vax_[DD-MMM-YYYY].mat*|*data*| See Section 3 of the README for the main Covid-19 ODE model [here](https://gitlab.com/tpmcovid/ode-model) |
|*vaccine_data_national_[YYYY-MM-DD].mat*|*data*| See Section 4 of the README for the main Covid-19 ODE model [here](https://gitlab.com/tpmcovid/ode-model) |
|*reshaped_b2_projections_final_[YYYY-MM-DD].csv*|*data*| This was constructed by combining historical vaccine uptake data with future uptake projections provided by the Ministry of Health. As these projections have now been superceded by more recent uptake data, this file is now obsolete and not needed. |
|*therapeutics_by_age_[DD-MMM-YYYY].mat*|*data*| ? |
|*hospOccDataTotals_MOH_[YYYY-MM-DD].xlsx*|*data*|Data from file *covid-cases-in-hospital-counts-location.xlsx* file in the *nz-covid-data/cases/* folder of the [MOH Github repo](https://github.com/minhealthnz/nz-covid-data/tree/main/cases)|
|*popsize_national.xlsx*|*data*| ? |
|*nzcontmatrix.xlsx*|*data*| This is taken from the paper 'Projecting social contact matrices in 152 countries using contact surveys and demographic data' by Prem et al. (specifically the overall age-and-location-specific contact matrix for New Zealand), and is available [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697). |
|*border_incidence.xlsx*|*data*| Data collected by M.Plank and stored in [this GitHub page](https://github.com/michaelplanknz/modelling-ba5-in-nz) |



# 3. How to process the unit record data
The model takes in a processed version of the *TPM_comm_cases_info_YYYY-MM-DD.csv* file, which you should have previously downloaded to the *processMOHdata* folder. To process it, follow these steps:
1. In the *processMOHdata* folder, open the *epi_data_cleaning.m* Matlab file.
2. Update the variable `readDate` to the datestamp at the end of the *TPM_comm_cases_info_YYYY-MM-DD.csv* file, which should be in the same folder
3. Run *epi_data_cleaning.m*
4. You should now have a *epidata_by_age_and_vax_DD-MMM-YYYY.mat* file in the `data` folder
5. In the main file *naiveABCSEIR.m*, change the `dateLbl` variable at the top of the script to the same datestamp as the one in the *data/epidata_by_age_and_vax_DD-MMM-YYYY.mat* file you just produced.

# 4. How to process the vaccination data
1. Make sure the vaccine line data from the sftp (*TPM_vaccine_YYYY-MM-DD.csv*) is in the `processVaxData` folder
2. If updating the vax data file to a more recent one, check that the latest *TPM_vaccine_YYYY-MM-DD.csv* file isn't smaller than the previous one (if it is this is a red flag that some records have disappeared, as happened previously!)
3. In `processVaxDataNational.m`, set the `datelbl_administered` variable to be the date stamp on the data file, in the form YYYY-MM-DD
4. Run `processVaxDataNational.m`, which will save a *vaccine_data_YYYY-MM-DD.csv* file in the `data` folder. This will contain a time series for the cumulative number of kth doses (1, 2, 3, 4+) given to people in each 5-year age band (the last age band being 75+ years)
5. It is strongly recommended to set `check45doses = 1`, this will plot timeseries of the cumulative doses (up to 5+ doses), which will help spot any inconsistencies in the data.
6. In the main file *naiveABCSEIR.m*, change the `dateVax` variable at the top of the script to the same datestamp as the one in the *data/vaccine_data_YYYY-MM-DD.csv* file you just produced.

# 5. How to run the ODE model

For all code, set the Matlab current directory to the Gitlab repo root folder 'ode-model'.

The main file `naiveABCSEIR.m` is divided into six sections, which can either be run all together by running the entire script, or separately. Running the entire script will produce the following outputs, which are described in more detail in the following sections:
* *results_Uni_Filtered100_DDMMMFit_DDMMMRun.mat*
* *violinPlots_DDMMMFit_DDMMMRun.png*
* *results_Uni_Filtered100_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_FilteredBest_SCENARIO NAME_DD-MMM-YYYY_fit.mat*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.png*
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit_ageSplit.png*

## **Section 1 - parameter initialisation**: 
This is where all the data is input and all the fixed and fitted model parameters are defined. Before running the script, **ensure the correct filenames have been defined** in this section, and that the corresponding files are all in the *data* folder. For a complete list of the required files and where to get them, please see the "Prerequisites" section of this README.


## **Section 2 - model fitting**
This section contains the simple ABC method used to narrow down the prior distribution of fitted parameters to the posterior distribution of values which resulted in the 1% best fitting model results. 
For a list of assumptions on the prior distributions, please see Table 1 of the *C_ODE_model_assumptions.pdf* document in this repo.

The fitting procedure can take several hours depending on the value chosen for `nSamples`. If using all 16 default fitted parameters, we recommend having at least `nSamples` = 15,000 to ensure enough combinations of parameter values get tested. Note that using 15,000 will take approximately 10 hours to run, but this estimate might vary depending on the available processing power.

The following outputs will be produced and saved in the `results` folder:
* *results_Uni_Filtered100_DDMMMFit_DDMMMRun.mat* : a Matlab structure containing  a feeld for each of the model output. Each feeld contains a matrix with 16 columns corresponding to 16 age groups and a line for each day simulated. For a description of the outputs, please see the preamble of the `functions/extractEpiVarsCompact.m` function. The output file is labeled with the datestamp of the line data used for the fit, and the date in which the fit was run.
* *violinPlots_DDMMMFit_DDMMMRun.png* : a figure representing the posterior distributions of the fitted parameters on a 0-1 scale. This figure can be used to gauge which parameter values gave the best 1% fit to the data.

## **Section 3 - loading of previous model fitting results (optional)**
If results from a previous model fitting are already available, it's possible to skip the first two sections and load the fitted posterior parameter set in this section. **Caution: if any of the model parameters have been changed, it is strongly recommended to re-run the fitting procedure.**
To load a previous model fit, change the `load100` parameter to `1` and update the `fIn` parameter with the path to the corresponding *results_Uni_Filtered100_DDMMMFit_DDMMMRun.mat* file.

In this section, it is also possible to display a table of posterior values corresponding to "best fitting" posterior, together with the interquartile range associated to each fitted parameter. To display this table, change the `showPosteriorStats` parameter to `1`.

## **Section 4 - checking vaccination rates (optional)**

This section plots the vaccination rates for the different scenarios presented in the paper.

## **Section 5 - scenario simulations**
Using the scenario parameter values defined in `getBasePar.m`, this is where scenarios are run from the pre-loaded posterior parameter sets. To define the different scenario, open the `getBasePar.m` function and scroll down to the last section of the script. There parameters that can currently be customised are the following:

* `policyDate`: date where a new policy affecting transmission takes place
* `cRampDeltaPolicy`: transmission multiplier associated with a policy change. If set to 1, no change will happen. It is both possible to run a single transmission multiplier for each policy, or a range of three transmission multipliers (low-middle-high) for each policy, as described in the script. However, running the latter kind of scenarios will only produce .csv result files, and no plots.
* `VOC2date`: date when a new variant of concern is assumed to become predominant in the country
* `VOC2active`: boolean variable, set to 1 to include a variant of concern on the chosen date
* `seasonMultAmp`: transmission multiplier associated with seasonality. It is set, by default, to 0.1, meaning that we assume a peak of 10% increase in transmission during winter months compared to transmission on 1st April, and a trough of 10% decrease in transmission during summer months. Note that the seasonality wave is set to be active from 1 April 2023, and any model results previous to this date will not include a seasonality effect.
* `antiviralsEffectIHRmult`: multiplier associated with the effect of antivirals. If set to 0, the antivirals are assumed to give no protection from hospitalisation, if set to 1, the antivirals are assumed to give perfect protection. Note that the effect of antivirals on fatality rate is a fitted parameter and is defined in `getBasePar.m`.

It is possible to test multiple scenarios at a time, by adding values to test to each parameter vector. 

Once the scenario parameter values have been defined in `getBasePar.m`, the script will run each scenario using all 1% best (posterior) sets of fitted parameter that have been produced or loaded earlier.
The following output file will then be saved in the `results` folder:
* *results_Uni_Filtered100_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for each of the posterior parameter sets, labeled by scenario name and datestamp of the data used for the model fitting. Note that this is an intermediary result file, which will be processed in the next section.

## **Section 6 - scenario result post-processing**
This is where the scenario results are processed to produce the 95% bands data and best-fit data. 
The following output files will be produced and saved in the `results` folder:
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for the 95% best fitting parameter sets out of the posterior parameter sets.
* *results_Uni_FilteredBest_SCENARIO NAME_DD-MMM-YYYY_fit.mat* : Matlab structure containing the results of the scenario simulations for the single best fitting parameter set.

If the `cRampDeltaPolicy` parameter was set to have a low-middle-high value for each policy tested, this section will merge the results for each of those three levels into a single result file.

## **Section 7 - plots**
This is where the aggregated and age-split plots for the scenario outputs are created. To run this, make sure that the variable `folderIn` points to the correct folder (if running the entire main script, or if the results to plot are kept in the `results` folder, simply set `folderIn = "results/"`). 

This script will output two figure in the `results` folder in .png format:
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit.png* : a tiled plot of aggregated model outputs and data for daily infections, daily reported cases, daily hospital admissions, hospital occupancy, and daily deaths
* *results_Uni_Filtered95_SCENARIO NAME_DD-MMM-YYYY_fit_ageSplit.png* : a tiled plot of age-split model outputs and data for daily reported cases, daily hospital admissions, daily deaths, ratio of hospital admissions / cases, case ascertainment rate, and ratio of cumulative infections per total population

Note that the plotting functions are not currently compatible with results produced using the low-middle-high values for the `cRampDeltaPolicy` scenarios parameter, this section will produce an error if trying to run it with those results, but this will not affect the rest of the model run.

# 6. Performing the sensitivity run

If you replace the population size file "popsize_national.xlsx" with "popproj_national2018-21.xlsx" (an alternative dataset from Stats NZ), the sensitivity run will be performed. This can be done by commenting line 26 and uncommenting line 27 of 'naiveABCSEIR.m'.



# Acknowledgements 
The authors acknowledge the role of the New Zealand Ministry of Health, StatsNZ, and the Institute of Environmental Science and Research in supplying data in support of this work, and the Covid-19 Modelling Government Steering for helping to design modelling questions and to interpret model outputs. The authors are grateful to Ning Hua and Rachel Owens at Precision Driven Health for code peer review, and to Nigel French, Emily Harvey, Markus Luczak-Roesch, Dion O'Neale, Matt Parry and Patricia Priest for feedback on the model. The authors acknowledge the contributions of Rachelle Binny, Shaun Hendy, Kannan Ridings, Nicholas Steyn and Leighton Watson to a previous model from which this model was developed.

