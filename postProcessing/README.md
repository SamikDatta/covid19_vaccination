# Post processing of results files

This folder contains some scripts used to process the .mat files produced by `naiveABCSEIR.m` into more readable file formats.

It contains the following functions:
- `writeMOHspreadsheetsTrajAndBands.m`
- `writeRMarkdownSpreadsheets.m`
- `writeRMarkdownSpreadsheetsCompare.m`
- `writeCSVsForMike.m`

## How to process results files into Excel spreadsheets (for TAS dashboard)
1. Open `writeMOHspreadsheetsTrajAndBands.m`
2. Change the variable `folderIn` to the path to the .mat results files produced by `naiveABCSEIR.m`. This script only reads in the results file *results_Uni_Filtered95_[SCENARIONAME]_[DD-MM-YYYY]_fit.mat*, labeled with the scenario name and the date of the data the model was fitted to 
3. Change the variable `folderOut` to the path where you want the Excel spreadsheets to be saved
4. Change the variable `scenario_labels_fIn` to the scenario name labels used in the .mat results files. By default, these will be *scenario1*, *scenario2*, etc.
5. Change the variable `scenario_labels_fOut` to the scenario name labels you'd like to have in the Excel spreadsheets' filenames
6. Run `writeMOHspreadsheetsTrajAndBands.m`

The script will produce five Excel spreadsheets in the folder specified:
- *infections_with_bands_[SCENARIONAME]_[DD-MM-YYYY]_fit.xlsx*
- *cases_with_bands_[SCENARIONAME]_[DD-MM-YYYY]_fit.xlsx*
- *newAdmissions_with_bands_[SCENARIONAME]_[DD-MM-YYYY]_fit.xlsx*
- *hospOcc_with_bands_[SCENARIONAME]_[DD-MM-YYYY]_fit.xlsx*
- *deaths_with_bands_[SCENARIONAME]_[DD-MM-YYYY]_fit.xlsx*

## How to process results files into Excel spreadsheets (for producing R Markdown dashboards, e.g. `reports/2023-05-01-May-2023-isolation-scenarios.Rmd`)
1. Open `writeRMarkdownSpreadsheets.m`
2. Change the variable `filename_base` to the base folder you want the Excel files to be saved in. Excel files for each scenario will be saved into a separate subfolder within this.
3. Change the variable `scenarios` to the names of the scenarios you have run (i.e. the name which distinguishes the outputs from the ODE model).
4. Change the variables `filenameBands` and `filenameBestFit` to match the format that the ODE model outputs are in, so that combined with the scenario names the correct files are accessed within each iteration of theloop.
5. Run `writeRMarkdownSpreadsheets.m`

The script will produce 13 Excel spreadsheets in each subfolder for one scenario:
- *all_infections_with_bands.xlsx*
- *cases_with_bands.xlsx*
- *cum_all_infections_with_bands.xlsx*
- *cum_first_infections_with_bands.xlsx*
- *cum_reinfections_with_bands.xlsx*
- *deaths_with_bands.xlsx*
- *first_infections_with_bands.xlsx*
- *hospOcc_with_bands.xlsx*
- *never_inf_with_bands.xlsx* (note that this has not been calculated correctly for never-infected, and `functions/epiVarsCompact` needs updating to calculate this.)
- *newAdmissions_with_bands.xlsx*
- *popn_sizes.xlsx*
- *reinfections_with_bands.xlsx*
- *weighted_sus_with_bands.xlsx*

## How to process results files into Excel spreadsheets (for producing R Markdown comparison documents, e.g. `reports/2023-04-28 Compare scenarios update.qmd`)
1. Open `writeRMarkdownSpreadsheetsCompare.m`
2. Change the variable `scenarios` to the names of the scenarios you have run (i.e. the name which distinguishes the outputs from the ODE model).
3. Change the variables `saveStr` and `saveStr_age` to the names of the Excel files you want to produce for the total and age-split numbers respectively.
4. Change the variable `filenameBands` to match the format that the ODE model outputs are in, so that combined with the scenario names the correct file is accessed within each iteration of theloop.
5. Run `writeRMarkdownSpreadsheetsCompare.m`

The script will produce two Excel spreadsheets in the `results` folder:
- filename `XXX.xlsx` (whatever you have made the variable `saveStr`).
- filename `YYY.xlsx` (whatever you have made the variable `saveStr_age`).

## How to process results files into CSV for comparing vaccination scenarios
1. Open `writeCSVsForMike.m`
2. Change the variable `filenameAll` to match the format that the ODE model outputs are in, so that combined with the scenario names the correct file is accessed within each iteration of theloop.
3. Change the variable `saveStr` to the name of the CSV file you want to produce for the differences between scenarios.
5. Run `writeCSVsForMike.m`
