# Post processing of results files

This folder contains some scripts used to process the .mat files produced by `naiveABCSEIR.m` into more readable file formats.

It contains the following functions:
- `writeRMarkdownSpreadsheetsCompare.m`

## How to process results files into Excel spreadsheets (for producing R Markdown for Figure 4)
1. Open `writeRMarkdownSpreadsheetsCompare.m`
2. Change the variable `scenarios` to the names of the scenarios you have run (i.e. the name which distinguishes the outputs from the ODE model).
3. Change the variables `saveStr` and `saveStr_age` to the names of the Excel files you want to produce for the total and age-split numbers respectively.
4. Change the variable `filenameBands` to match the format that the ODE model outputs are in, so that combined with the scenario names the correct file is accessed within each iteration of the loop.
5. Run `writeRMarkdownSpreadsheetsCompare.m`.

The script will produce two Excel spreadsheets in the `results` folder:
- filename `XXX.xlsx` (whatever you have made the variable `saveStr`).
- filename `YYY.xlsx` (whatever you have made the variable `saveStr_age`).

## How to produce Figure 4

1. Open the file `Figure 4 markdown.qmd`.
2. Change the variables `filename_data` and `filename_age` to match the filenames you chose in the script `writeRMarkdownSpreadsheetsCompare.m`. This should then run to produce a HTML file containing Figure 4.
