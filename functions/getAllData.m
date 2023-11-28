function dataComb = getAllData(myDataPath, dataFileNames)
% Function that reads in data on daily cases, daily hospital admissions,
% daily deaths, hospital occupancy, and border incidence. Data is split by
% vaccination status and age group
% INPUTS:
% - myDataPath: path to folder containing data
% - dataFileNames: structure with the following fields
%       - epiDataFname: filename of cases, hosp. admissions and deaths datafile
%       - hospOccFname: filename of hosp occupancy datafile
%       - borderIncFname: filename of border incidence datafile
% OUTPUT:
% - dataComb: table containing all data combined


% Importing line data on COVID cases, hosp. admissions and deaths
fName = myDataPath + dataFileNames.epiDataFname;
readIn = load(fName);
epiData = readIn.outTab;

% Create new fields for total daily cases and deaths
epiData.nCasesData = sum(epiData.nCases_v0 + epiData.nCases_v1 + epiData.nCases_v2 + epiData.nCases_v3, 2); 
epiData.nHospData = sum(epiData.nHosp_strict_v0_byDateOfAdmission + epiData.nHosp_strict_v1_byDateOfAdmission + epiData.nHosp_strict_v2_byDateOfAdmission + epiData.nHosp_strict_v3_byDateOfAdmission, 2);
epiData.nDeathsData = sum(epiData.nDeaths_strict_v0_byDateOfDeath + epiData.nDeaths_strict_v1_byDateOfDeath + epiData.nDeaths_strict_v2_byDateOfDeath + epiData.nDeaths_strict_v3_byDateOfDeath, 2); 

% Importing line data on hospital occupancy
fName = myDataPath + dataFileNames.hospOccFname;
hospData = readtable(fName, "Sheet", "NZ total");
hospData.Properties.VariableNames = ["date", "hospOccTotalMOH"];
hospData = hospData(datenum(hospData.date) >= datenum('25JAN2022'), :);
hospData.hospOccTotalMOH = str2double(hospData.hospOccTotalMOH);

% Importing data on incidence in border workers
fName = myDataPath + dataFileNames.borderIncFname;
borderData = readtable(fName);
borderData.WeekEnding = datetime(string(borderData.WeekEnding));
borderData = borderData(~isnat(borderData.WeekEnding), :);
borderData.date = borderData.WeekEnding-7;          % Tests week ending X roughly indicative of new infections around X-7 

% Merge into one table 
tmp = outerjoin(epiData, hospData, 'Keys', 'date', 'MergeKeys', true);
dataComb = outerjoin(tmp, borderData, 'Keys', 'date', 'MergeKeys', true);

end
