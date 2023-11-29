%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            ODE model for Omicron outbreak in New Zealand
%       with parameter uncertainty quantification via simple ABC
%                        PART I. Model fitting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('functions');

% For reproducibility
rng(34604);

dateLbl = "13-Aug-2023";                    % Datestamp of line data used for fitting and plotting


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
% Save model output, parameter values and distance function values for accepted simulations
save(fOut, 't', 'epiVarsCompact', 'Theta', 'parBase', ...
    'pUniFiltered100', 'dUniFiltered100', '-v7.3');

