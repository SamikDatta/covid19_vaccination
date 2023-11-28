%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            ODE model for Omicron outbreak in New Zealand
%       with parameter uncertainty quantification via simple ABC
% SENSITIVITY RUN FOR DIFFERENT popSizeFname - different population data used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

myDataPath = 'data/';
addpath('functions');

%% 1. Global variables initialisation
% List of parameters to fit to data
parsToFit = {'dateSeed', 'Cstart', 'Cramp', 'rampDays', 'rampStart', ...
    'pTestMult', 'IFR', 'IHR', 'relaxAlpha', 'MRampDays', 'vocWane', ...
    'waneRate', 'Cramp2', 'ramp2Days', 'ramp2Start', 'aViralEffect'};
nParsToFit = length(parsToFit);

% Filenames of input datafiles
dateLbl = "13-Aug-2023";                    % Datestamp of line data used for fitting and plotting
runDate = string(datetime("today", 'Format', 'ddMMM'));        % Set this to today's date if running whole script from start to end, or to a previously saved date if only rerunning from sec. 3 without redoing model fit 

epiDataFname = "epidata_by_age_and_vax_" + dateLbl + ".mat";    % Line data
vaxDataFname = "vaccine_data_national_2023-06-06";              % Vax data
vaxProjFname = "reshaped_b2_projections_final_2022-07-13.csv";  % Vax projections (currently not used)
AVdataFname = "therapeutics_by_age_14-Aug-2023.mat";            % Antiviral data
hospOccFname = "covid-cases-in-hospital-counts-location-16-Aug-2023.xlsx";           % Only used for plotting
popSizeFname = "popsize_national.xlsx";                         % NZ population structure
% popSizeFname = "popproj_national2018-21.xlsx";                % NZ population structure (alternative from Stats NZ)
CMdataFname = "nzcontmatrix.xlsx";                              % Prem contact matrix
borderIncFname = "border_incidence.xlsx";

% Import data used for fitting and plotting
dataComb = getAllData(myDataPath, epiDataFname, hospOccFname, borderIncFname);
tMaxData = datenum(dataComb.date(find(~isnan(dataComb.nCasesData), 1, 'last')));

% Create structure of base parameters that do not change from one
% realisation to the next - run for period with data (tMaxData)
parBase = getBasePar(tMaxData, myDataPath, popSizeFname, CMdataFname, ...
    vaxDataFname, vaxProjFname, AVdataFname);


%% 2. Model fitting

% Define which data to fit (currently daily cases, daily deaths, daily new
% hosp. admissions, daily incidence per capita, age breakdown of daily
% cases, age breakdown of daily hosp. admissions):
dataVars = {'nCasesData','nDeathsData','nHospData', 'NationalBorder', ...
    'CasesByAge', 'HospByAge'};

% Define how many random parameter samples to use for the model fitting
nSamples = 15000; %15000 is our usual nSample, takes about 10h to run
pUni = rand(nSamples, length(parsToFit));
dUni = zeros(nSamples, 1);

% Time vector
t = parBase.tBase;

% For each fitted parameter sample, solve the ODE and calculate distance
% function:
parfor iSample = 1:nSamples
    fprintf('Sample %i of %i\n', iSample, nSamples)

    % Solve ODE with current fitted parameter sample
    ThetaTemp = array2table(pUni(iSample ,:), 'VariableNames', parsToFit);
    parInd(iSample) = getParUnified(ThetaTemp, parBase);

    % Create temporary structure with merged fields from parBase and parInd
    % because parfor won't allow index variables in an anonymous
    % function:
    parTemp = catstruct(parBase, parInd(iSample));

    % Initial conditions and ODE options
    IC = getIC(parTemp);
    odeOptions = odeset('NonNegative', ones(size(IC))');

    % Solve ODE
    [~, Y] = ode45(@(t, y)myODEs2(t, y, parTemp), t, IC, odeOptions);

    % Retain a compact structure of variables rather than the full Y
    epiVarsCompact(iSample) = extractEpiVarsCompact(t, Y, parTemp);


    % Calculate distance function from ODE results
    [~, ~, dist] = calcErrorWithAgeCurrent(t, epiVarsCompact(iSample), dataComb, parTemp, dataVars);
    errWeights = ones(1, length(dataVars));
    dUni(iSample, 1) = sum(errWeights .* dist, 2);

end

% Filter Uni to get 'posterior'
qtol = 0.01; % proportion to retain 1% best
condition = dUni <= quantile(dUni, qtol);

% Overwrite epiVarsCompact with retained posterior sets
epiVarsCompact = epiVarsCompact(condition);

% Keep copy of retained. 100 refers to (approx) '100% region of posterior'
epiVarsCompact100 = epiVarsCompact;
pUniFiltered100 = pUni(condition,:);
dUniFiltered100 = dUni(condition);

% Plot violin plots of posterior distribution of parameters
f = figure(1);
f.Position = [50 200 1000 400];
vs = violinplot(pUniFiltered100, cellstr(parsToFit));
if parBase.sensitivity_flag == 0
    saveas(f, sprintf('results/violinPlots_%sFit_%sRun.png', ...
        string(datetime(dateLbl, 'Format', 'ddMMM')), ...
        string(datetime("today", 'Format', 'ddMMM'))));
else % add '_sensitivity' to filename if needed
    saveas(f, sprintf('results/violinPlots_%sFit_%sRun_sensitivity.png', ...
        string(datetime(dateLbl, 'Format', 'ddMMM')), ...
        string(datetime("today", 'Format', 'ddMMM'))));
end

% Save filtered posterior
Theta = array2table(pUniFiltered100, 'VariableNames', parsToFit);

% Save results of all samples run into a mat file
if parBase.sensitivity_flag == 0
fOut = sprintf('results/results_Uni_Filtered100_%sFit_%sRun.mat', ...
    string(datetime(dateLbl, 'Format', 'ddMMM')), ...
    string(datetime("today", 'Format', 'ddMMM')));
else % add '_sensitivity' to filename if needed
fOut = sprintf('results/results_Uni_Filtered100_%sFit_%sRun_sensitivity.mat', ...
    string(datetime(dateLbl, 'Format', 'ddMMM')), ...
    string(datetime("today", 'Format', 'ddMMM')));
end
save(fOut, 't', 'epiVarsCompact', 'Theta', 'parBase', ...
    'pUniFiltered100', 'dUniFiltered100', '-v7.3');

%% 3. Posterior distribution loading
% (can start from here if fitting already run)
load100 = 1;
if load100 == 1
    fIn = sprintf('results/results_Uni_Filtered100_%sFit_%sRun.mat', string(datetime(dateLbl, 'Format', 'ddMMM')), runDate);
    load(fIn);
    fprintf('\nLoading filtered 100...\nOriginal number samples assuming 1%% acceptance rate = %i\n', ...
        size(epiVarsCompact,2)/0.01)
end

% Change to true for a table of posterior values for all fitted parameters
showPosteriorStats = 1;
if showPosteriorStats == 1
    [~, mi] = min(dUniFiltered100); % Index of best fit run
    getStatIndex = mi; % Change to whichever run to run posterior stats on
    posteriorStats = getPosteriorStats(pUniFiltered100, dUniFiltered100, parsToFit, getStatIndex);
    disp(posteriorStats)
end


%% 4. Plotting vaccination rates

clear parInd

% End of simulation date
tModelRunTo = datenum('06-Aug-2023'); % needs to be 6th August 2023 or later, otherwise crashes (related to antivirals)

% Structure of base parameters based on running into future
parBase = getBasePar(tModelRunTo, myDataPath, popSizeFname, CMdataFname, ...
    vaxDataFname, vaxProjFname, AVdataFname);
t = parBase.tBase;


%% Let us check the vaccination rates for all scenarios

all_dates = datetime(parBase.tBase,'ConvertFrom','datenum');

vac_table = table('Size', [parBase.nScenarios 6], ...
    'VariableTypes', ["string", repmat("cell", 1, 5)], ...
    'VariableNames', ["scenario", "first dose", "second dose", ...
    "third dose", "fourth dose", "total doses"]);
vac_table.scenario = ["baseline", "No vaccination", "No antivirals", "No vaccination or antivirals", ...
    "No vaccination in U60s", "90% Vaccination", "reduced coverage in older ages", ...
    "Maori vaccination", "European / Other vaccination"]';

ThetaTemp = array2table(pUniFiltered100(1, :), 'VariableNames', parsToFit);

for i = 1:parBase.nScenarios

    parScen = parBase.scenarios(i, :);
    parInd = getParUnified(ThetaTemp, parBase, parScen);
    vac_table.("first dose"){i} = parInd.nDoses1Smoothed;
    vac_table.("second dose"){i} = parInd.nDoses2Smoothed;
    vac_table.("third dose"){i} = parInd.nDoses3Smoothed;
    vac_table.("fourth dose"){i} = parInd.nDoses4Smoothed;
    vac_table.("total doses"){i} = vac_table.("first dose"){i} + vac_table.("second dose"){i} + ...
        vac_table.("third dose"){i} + vac_table.("fourth dose"){i};

end

% Plot cumulative doses for each scenario

legLbls = ["Baseline", "10% drop", "20-25 yo rates", "Euro rates", "Māori rates"];
f = figure(2);  
f.Position = [100 100 1200 800];
titles = {'0-5', '5-10', '10-15', '15-20', '20-25', '25-30', ...
    '30-35', '35-40', '40-45', '45-50', '50-55', '55-60', '60-65', ...
    '65-70', '70-75', '75+'};
tiledlayout(4, 4);
for ag = 1:16
    nexttile
    hold on
    title(titles(ag))
    plot(all_dates, cumsum(vac_table.("total doses"){1}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){6}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){7}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){9}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){8}(:, ag))/parBase.popCount(ag) );
    hold off
    xlim([datetime(2021, 1, 1), datetime(2023, 7, 1)] )
    ylim([0 4])
    grid on
    if ag == 1
        legend(legLbls, 'Location', 'NorthWest');
    end
    if mod(ag, 4) == 1
        ylabel('doses per capita')
    end
end

% Pull some summary ethnicity stats for paper
ti = find(t == datenum('01FEB2022'));
AM = sum(vac_table.("second dose"){8}(1:ti,:))
AE = sum(vac_table.("second dose"){9}(1:ti,:))
pM1 = sum(AM(4:13))/sum(parBase.popCount(4:13))
pM2 = sum(AM(14:end))/sum(parBase.popCount(14:end))
pE1 = sum(AE(4:13))/sum(parBase.popCount(4:13))
pE2 = sum(AE(14:end))/sum(parBase.popCount(14:end))

BM = sum(vac_table.("third dose"){8}(1:ti,:))
BE = sum(vac_table.("third dose"){9}(1:ti,:))
qM1 = sum(BM(4:13))/sum(parBase.popCount(4:13))
qM2 = sum(BM(14:end))/sum(parBase.popCount(14:end))
qE1 = sum(BE(4:13))/sum(parBase.popCount(4:13))
qE2 = sum(BE(14:end))/sum(parBase.popCount(14:end))


%% 5. Scenario simulations

%% Now onto the actual runs

nScenarios = parBase.nScenarios; % number of scenarios

for i = 1:nScenarios

    % Table containing scenario-specific parameters
    parScen = parBase.scenarios(i, :);

    % Run current scenario for all best 1% fitted runs
    parfor iSample = 1:size(pUniFiltered100, 1)

        fprintf('Scenario %i of %i, sample %i of %i\n', i, nScenarios, iSample, size(pUniFiltered100,1))

        % solve ODE with optimal (not super efficient but OK)
        ThetaTemp = array2table(pUniFiltered100(iSample, :), 'VariableNames', parsToFit);

        % parameter structure for ith realisation
        parInd = getParUnified(ThetaTemp, parBase, parScen);

        parTemp = catstruct(parBase, parInd);

        if i == 3 || i == 4 % for scenarios 3 and 4, no antivirals, so reduce effect to zero
            parTemp.antiviralsEffectIFRmult = 0;
        end

        % Initial condition
        IC = getIC(parTemp);
        odeOptions = odeset('NonNegative', ones(size(IC))' );
        [~, Y] = ode45(@(t, y) myODEs2(t, y, parTemp), t, IC, odeOptions);

        % retain a compact structure of variables rather than the full Y
        epiVarsCompact(iSample) = extractEpiVarsCompact(t, Y, parTemp);
        epiVarsCompact_vaxSplitInf(iSample) = extractEpiVarsCompact_vaxSplitInf(t, Y, parTemp);

    end

    epiVarsCompact100Future = epiVarsCompact;
    epiVarsCompact_vaxSplitInf100Future = epiVarsCompact_vaxSplitInf;


    Theta = array2table(pUniFiltered100, 'VariableNames', parsToFit);

    if parBase.sensitivity_flag == 0
        fOut = sprintf('results/results_Uni_Filtered100_%s_%s_fit.mat', parBase.scenario_names{i}, dateLbl);
    else % add '_sensitivity' to filename if needed
        fOut = sprintf('results/results_Uni_Filtered100_%s_%s_fit_sensitivity.mat', parBase.scenario_names{i}, dateLbl);
    end

    save(fOut, 't', 'epiVarsCompact', 'Theta', 'parBase', 'epiVarsCompact100Future', 'epiVarsCompact_vaxSplitInf100Future', 'pUniFiltered100', 'dUniFiltered100','-v7.3');
end


%% 6. Get and save scenario bands and bestFit results

if size(parBase.cRampDeltaPolicy, 1) == 1

    labels = ["95", "Best"];
    conditions = [0.95, 0.0001];

    for i = 1:parBase.nScenarios

        % Import data for scenario i
        scenName = parBase.scenario_names{i};
        if parBase.sensitivity_flag == 0
            fIn = sprintf('results/results_Uni_Filtered100_%s_%s_fit.mat', ...
                parBase.scenario_names{i},dateLbl);
        else
            fIn = sprintf('results/results_Uni_Filtered100_%s_%s_fit_sensitivity.mat', ...
                parBase.scenario_names{i},dateLbl);
        end

        load(fIn);

        for ii = 1:length(labels)

            % --Filter uni
            if labels(ii) ~= "Best"
                condition = dUniFiltered100 <= quantile(dUniFiltered100, conditions(ii));
            else % add '_sensitivity' to filename if needed
                [~, condition] = min(dUniFiltered100);
            end

            epiVarsCompact = epiVarsCompact100Future(condition);
            epiVarsCompact_vaxSplitInf = epiVarsCompact_vaxSplitInf100Future(condition);
            pUniFiltered = pUniFiltered100(condition,:);
            dUniFiltered = dUniFiltered100(condition);

            % save filtered uni
            Theta = array2table(pUniFiltered, 'VariableNames', parsToFit);
            if parBase.sensitivity_flag == 0
                fOut = sprintf('results/results_Uni_Filtered%s_%s_%s_fit.mat', labels(ii), scenName, dateLbl);
            else % add '_sensitivity' to filename if needed
                fOut = sprintf('results/results_Uni_Filtered%s_%s_%s_fit_sensitivity.mat', labels(ii), scenName, dateLbl);
            end
            save(fOut, 't', 'epiVarsCompact', 'epiVarsCompact_vaxSplitInf', 'Theta', 'parBase', 'pUniFiltered', 'dUniFiltered', '-v7.3');
        end

    end

    % Merge results on different levels of transmission increase if several
    % were run, i.e. collate results for the lower and upper bounds of a
    % transmission increase into a single result file. NOTE: plots won't work
    % if this is the case. Only use to produce results spreadsheets
elseif size(parBase.cRampDeltaPolicy, 1) > 1
    for seasi = 1:length(seasonMultAmp)
        for policyi = 1:size(cRampDeltaPolicy, 2)
            mergedt = zeros();
            mergedEpiVarsBest = struct([]);
            mergedEpiVarsBands = struct([]);
            for leveli = 1:size(cRampDeltaPolicy, 1)
                fIn = sprintf("results/results_Uni_Filtered100_%.1fseasonMult_%+.1f%%transInc_%s_fit", seasonMultAmp(seasi), 100*(cRampDeltaPolicy(leveli, policyi)-1), dateLbl);
                load(fIn, 't', 'epiVarsCompact');
                mergedt = [mergedt, t];
                mergedEpiVarsBands = [mergedEpiVarsBest, epiVarsCompact(dUniFiltered100 <= quantile(dUniFiltered100, 0.95))];
                mergedEpiVarsBest = [mergedEpiVarsBands, epiVarsCompact(dUniFiltered100 <= quantile(dUniFiltered100, 0.0000001))];
            end
            if parBase.sensitivity_flag == 0
                fOutBest = sprintf('results/results_Uni_FilteredBest_%s_fit.mat', policyNames(policyi));
                fOutBands = sprintf('results/results_Uni_Filtered95_%s_fit.mat', policyNames(policyi));
            else % add '_sensitivity' to filename if needed
                fOutBest = sprintf('results/results_Uni_FilteredBest_%s_fit_sensitivity.mat', policyNames(policyi));
                fOutBands = sprintf('results/results_Uni_Filtered95_%s_fit_sensitivity.mat', policyNames(policyi));
            end

            save(fOutBest, 'mergedt', 'mergedEpiVarsBest', '-v7.3');
            save(fOutBands, 'mergedt', 'mergedEpiVarsBands', '-v7.3');
        end
    end
end


%% 7. Plots

overwriteFig = true;   % If true, any previous figure with the same name will be overwritten
folderIn = 'results/'; % Folder where results are stored
dataComb = getAllData(myDataPath, epiDataFname, hospOccFname, borderIncFname); % data for plotting

% Plot aggregated and age-split plots for each scenario
for i = 1:parBase.nScenarios
    if parBase.sensitivity_flag == 0
        filenameBands = sprintf('%sresults_%s_%s_fit.mat', folderIn, ...
            ['Uni_Filtered95_', parBase.scenario_names{i}],dateLbl);
        filenameBestFit = sprintf('%sresults_%s_%s_fit.mat', folderIn, ...
            ['Uni_FilteredBest_', parBase.scenario_names{i}], dateLbl);
    else % add '_sensitivity' to filename if needed
        filenameBands = sprintf('%sresults_%s_%s_fit_sensitivity.mat', folderIn, ...
            ['Uni_Filtered95_', parBase.scenario_names{i}],dateLbl);
        filenameBestFit = sprintf('%sresults_%s_%s_fit_sensitivity.mat', folderIn, ...
            ['Uni_FilteredBest_', parBase.scenario_names{i}], dateLbl);
    end

    % Aggregated plots
    plotTrajModelTiled(filenameBands, filenameBestFit, dateLbl, ...
        dateLbl, overwriteFig, dataComb)

    % Age-split plots for cases, hosp.adm. and deaths
    plotTrajModelAgeSplit(filenameBands, filenameBestFit, dataComb, ...
        dateLbl, parsToFit, i, overwriteFig)

end
