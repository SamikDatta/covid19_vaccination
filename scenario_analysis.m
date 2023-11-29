clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with Matlab functions
addpath('functions');

% Quantiles for summary stats [median + 95% CI]
Alpha = 0.05;
qt = [Alpha/2, 0.5, 1-Alpha/2];    


% Specify file names:
modelCase = "base";              % set this to "base" to analyse the base run with HSU pop data, or "sensitivity" to analyse the sensitivity analysis eith StatsNZ ERP data

if modelCase == "base"
    fNameModel = 'results/model_output.csv ';   % Model results file name to input
    fNameOut = 'latex/results_table.tex';       % File name for saving table latex
    fNamePop = 'data/HSU_by_age_eth.csv';       % HSU population file name
elseif modelCase == "sensitivity"
    fNameModel = 'results/model_output_sensitivity_ERP.csv ';  % Model results file name to input
    fNameOut = 'latex/results_table_sensitivity_ERP.tex';   % File name for saving table latex
    fNamePop = 'data/popproj2018-21.csv';                   % StatsNZ projected population file name
end

fNameDeaths = 'data/actual_deaths_by_age_and_sex.csv';   % Actual deaths by age (1 year bands) and sex file name
fNameMaori = 'data/maori_outcomes_by_age';              % Actual Maori hospitalisations and deaths by age (5 year bands) file name
fNameLifeTables = 'data/nz-complete-cohort-life-tables-1876-2021.csv';   % Life tables file name  https://www.stats.govt.nz/information-releases/new-zealand-cohort-life-tables-march-2023-update/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in model results
inTab = readtable(fNameModel);
nRows = height(inTab);

% Import life tables and actual male and female deaths in one year age band, exM and exF are life expectancy at age for male and female respectively
[age, exM, exF] = importLifeTables(fNameLifeTables);

deathsTab = readtable(fNameDeaths);
assert(isequal(age, deathsTab.age));

% Read in Maori deaths in model age bands
ageModel = (0:5:75)';
maoriOutcomes = readtable(fNameMaori);

% Import population data
[popData, ethNames] = getPopData(fNamePop);
% Create 2nd version of pop size data with over 75s pooled
pop16 = poolOver75s(popData, ethNames);


% Define compartor scenario for each scenario to comapre to:
scenarioLabels = categorical(["Baseline", "No vaccine", "No AVs", "No vaccine or AVs", "No vaccine in U60s", "10% drop in rates", "20-25-year-old rates", "Maori rates", "Euro/other rates"]);
nScenarios = length(scenarioLabels);
scenarioBase = repmat("Baseline", 1 , nScenarios);
scenarioBase(scenarioLabels == "No vaccine or AVs") = "No AVs";
scenarioBase(scenarioLabels == "Maori rates") = "Euro/other rates";


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert model output table into a more usable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdl.scenario = scenarioLabels(inTab.scenario)';
mdl.simulation = inTab.Simulation;
mdl.nInf = table2array(inTab(:, 3:18));
mdl.nInf1 = table2array(inTab(:, 19:34));
mdl.nHosp = table2array(inTab(:, 35:50));
mdl.nDeaths = table2array(inTab(:, 51:66));
mdl.peakOcc = inTab.peakOcc;
mdl = struct2table(mdl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate some rates for each row in table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% age-specific attack rates 
mdl.AR = mdl.nInf1 ./ pop16.Total';

% age-specific hops and deaths per 1st infection
mdl.I1HR = mdl.nHosp./mdl.nInf1;
mdl.I1FR = mdl.nDeaths./mdl.nInf1;

% age-specific per capita hosp and mortality rates 
mdl.hospRate = mdl.nHosp./pop16.Total';
mdl.deathRate = mdl.nDeaths./pop16.Total';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate aggregates and YLL for each row in table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdl.nInfTot = sum(mdl.nInf, 2);
mdl.nInf1Tot = sum(mdl.nInf1, 2);
mdl.nHospTot = sum(mdl.nHosp, 2);
mdl.nDeathsTot = sum(mdl.nDeaths, 2);

mdl.hospRateTot = mdl.nHospTot/sum(pop16.Total);
mdl.deathRateTot = mdl.nDeathsTot/sum(pop16.Total);

mdl.ARTot = mdl.nInf1Tot ./ sum(pop16.Total);


% calculate YLL for each row in table
exM = exM(:, 2);             % get lower and upper ends of CI
exF = exF(:, 2);   



% Calculate YLL based on actual deaths data and model deaths (for each row in table) 
mdl.YLL = calcYLL(age, exM, exF,  deathsTab.nDeathsActualM,  deathsTab.nDeathsActualF, ageModel, mdl.nDeaths')';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Deltas relative to comparator scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdl.dInf = nan(nRows, 1);
mdl.dHosp = nan(nRows, 1);
mdl.dDeaths = nan(nRows, 1);
mdl.dAR = nan(nRows, 1);
mdl.dHospRateTot = nan(nRows, 1);
mdl.dDeathRateTot = nan(nRows, 1);
mdl.dYLL = nan(nRows, 1);

mdl.dHospMaori1 = nan(nRows, 1);
mdl.dDeathsMaori1 = nan(nRows, 1);
mdl.dHospMaori2 = nan(nRows, 1);
mdl.dDeathsMaori2 = nan(nRows, 1);

for iScenario = 1:nScenarios
    currFlag = mdl.scenario == scenarioLabels(iScenario);
    baseFlag = mdl.scenario == scenarioBase(iScenario);
    mdl.dInf(currFlag) = mdl.nInfTot(currFlag) - mdl.nInfTot(baseFlag);
    mdl.dHosp(currFlag) = mdl.nHospTot(currFlag) - mdl.nHospTot(baseFlag);
    mdl.dDeaths(currFlag) = mdl.nDeathsTot(currFlag) - mdl.nDeathsTot(baseFlag);
    mdl.dAR(currFlag) = mdl.ARTot(currFlag) - mdl.ARTot(baseFlag);
    mdl.dHospRateTot(currFlag) = mdl.hospRateTot(currFlag) - mdl.hospRateTot(baseFlag);
    mdl.dDeathRateTot(currFlag) = mdl.deathRateTot(currFlag) - mdl.deathRateTot(baseFlag);
    mdl.dYLL(currFlag) = mdl.YLL(currFlag) - mdl.YLL(baseFlag);
    if scenarioLabels(iScenario) == "Maori rates"
        % Method 1 - calculate Delta in age-specific hosp rates and apply
        % to Maori pop size
        dHR = mdl.hospRate(currFlag, :)-mdl.hospRate(baseFlag, :);
        dFR = mdl.deathRate(currFlag, :)-mdl.deathRate(baseFlag, :);
        mdl.dHospMaori1(currFlag, :) = sum(dHR.*(pop16.Maori'), 2);
        mdl.dDeathsMaori1(currFlag, :) = sum(dFR.*(pop16.Maori'), 2);
        % Method 2 - calculate age-specific ratio of Euro to Maori rates and apply to actual Maori deaths
        rHR = mdl.hospRate(baseFlag, :)./mdl.hospRate(currFlag, :);
        rFR = mdl.deathRate(baseFlag, :)./mdl.deathRate(currFlag, :);
        mdl.dHospMaori2(currFlag, :) = sum((1-rHR).*maoriOutcomes.nMaoriHosp', 2);
        mdl.dDeathsMaori2(currFlag, :) = sum((1-rFR).*maoriOutcomes.nMaoriDeaths', 2);
        
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate table of summary stats for each scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.scenario = scenarioLabels';
for iScenario = 1:nScenarios
    currFlag = mdl.scenario == scenarioLabels(iScenario);
    results.nInfTot(iScenario, :) = quantile(mdl.nInfTot(currFlag), qt);
    results.nHospTot(iScenario, :) = quantile(mdl.nHospTot(currFlag), qt);
    results.nDeathsTot(iScenario, :) = quantile(mdl.nDeathsTot(currFlag), qt);
    results.ARTot(iScenario, :) = quantile(mdl.ARTot(currFlag), qt);
    results.YLL(iScenario, :) = quantile(mdl.YLL(currFlag), qt);
    results.peakOcc(iScenario, :) = quantile(mdl.peakOcc(currFlag), qt);
    results.dInf(iScenario, :) = quantile(mdl.dInf(currFlag), qt);
    results.dHosp(iScenario, :) = quantile(mdl.dHosp(currFlag), qt);
    results.dDeaths(iScenario, :) = quantile(mdl.dDeaths(currFlag), qt);
    results.dAR(iScenario, :) = quantile(mdl.dAR(currFlag), qt);
    results.dHospRateTot(iScenario, :) = quantile(mdl.dHospRateTot(currFlag), qt);
    results.dDeathRateTot(iScenario, :) = quantile(mdl.dDeathRateTot(currFlag), qt);
    results.dYLL(iScenario, :) = quantile(mdl.dYLL(currFlag), qt);
    results.dHospMaori1(iScenario, :) = quantile(mdl.dHospMaori1(currFlag), qt);
    results.dDeathsMaori1(iScenario, :) = quantile(mdl.dDeathsMaori1(currFlag), qt);
    results.dHospMaori2(iScenario, :) = quantile(mdl.dHospMaori2(currFlag), qt);
    results.dDeathsMaori2(iScenario, :) = quantile(mdl.dDeathsMaori2(currFlag), qt);
end
results = struct2table(results);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate latex table output and plot figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeTableLatex(results, fNameOut);


% Plotting
[~, scenToPlot] = sort(results.nDeathsTot(:, 2), 'ascend');     % Want plot to be sorted in order of median deaths
scenToPlot = scenToPlot(ismember(scenToPlot, 1:7));     % Only include scenarios 1-7
nBars = length(scenToPlot);

figure(1)
bar(1:nBars, results.nDeathsTot(scenToPlot, 2))
hold on
er = errorbar(1:nBars, results.nDeathsTot(scenToPlot, 2), results.nDeathsTot(scenToPlot, 2)-results.nDeathsTot(scenToPlot, 1), results.nDeathsTot(scenToPlot, 3)-results.nDeathsTot(scenToPlot, 2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
yline(3196, 'k--')
ylabel('deaths')
h = gca;
h.XTickLabel = results.scenario(scenToPlot);


