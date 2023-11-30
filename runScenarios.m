%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            ODE model for Omicron outbreak in New Zealand
%       with parameter uncertainty quantification via simple ABC
%                      PART II. Scenario simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('functions');

% For reproducibility
rng(34604);

dateLbl = "13-Aug-2023";                    % Datestamp of line data used for fitting and plotting
fitLbl = "30Nov";                           % Datestamp on the file containing the results of the model fitting orutine (date the model fitting routine was run)

%% 1. Initialisation

% Get a list of names of the parameters to be fitted to data
parsToFit = getParsToFit();

% Get data file names
[myDataPath, dataFileNames] = getDataFileNames(dateLbl);

% Import data used for fitting and plotting
dataComb = getAllData(myDataPath, dataFileNames);
tMaxData = datenum(dataComb.date(find(~isnan(dataComb.nCasesData), 1, 'last')));

% Create structure of base parameters that do not change from one
% realisation to the next - run for period with data (tMaxData)
parBase = getBasePar(tMaxData, myDataPath, dataFileNames);






%% 2. Posterior distribution loading
fIn = sprintf('results/results_Uni_Filtered100_%sFit_%sRun.mat', string(datetime(dateLbl, 'Format', 'ddMMM')), fitLbl);
load(fIn);
fprintf('\nLoading filtered 100...\nOriginal number samples assuming 1%% acceptance rate = %i\n', ...
    size(epiVarsCompact,2)/0.01)

% Change to true for a table of posterior values for all fitted parameters
showPosteriorStats = 1;
if showPosteriorStats == 1
    [~, mi] = min(dUniFiltered100); % Index of best fit run
    getStatIndex = mi; % Change to whichever run to run posterior stats on
    posteriorStats = getPosteriorStats(pUniFiltered100, dUniFiltered100, parsToFit, getStatIndex);
    disp(posteriorStats)
end




%% 3. Scenario simulations

% End of simulation date
tModelRunTo = datenum('06-Aug-2023'); % needs to be 6th August 2023 or later, otherwise crashes (related to antivirals)

% Structure of base parameters based on running into future
parBase = getBasePar(tModelRunTo, myDataPath, dataFileNames);


% Plot vaccine doses vs time for each scenario
makeVaccinationPlots(parBase, parsToFit)

nScenarios = parBase.nScenarios; % number of scenarios

t = parBase.tBase;
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


%% 4. Get and save scenario bands and bestFit results

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


%% 5. Plots

overwriteFig = true;   % If true, any previous figure with the same name will be overwritten
folderIn = 'results/'; % Folder where results are stored
dataComb = getAllData(myDataPath, dataFileNames); % data for plotting

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
