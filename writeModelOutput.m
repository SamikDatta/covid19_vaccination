%% Analysis to compare scenarios

% Base folder to start in - make this the Gitlab repo root folder 'ode-model'

clearvars % remove all variables
close all % close all figures

% Latency parameter needed for analysis
parBase.tE = 1; % Latent period (set manually so model fit not needed)

% Scenario parameters

scenarios = [1:9]';
last_date = '30-June-2023'; % date to save to

% Get info from first file

filenameAll = sprintf('results/results_Uni_Filtered100_scenario%d_13-Aug-2023_fit.mat', scenarios(1));

% Load first file
load(filenameAll, 't'); % load time vector
ind = find(t == datenum(last_date)); % index of the desired date
num_simulations = 150; % manually setting number of runs per scenario


%% Code for generating necessary numbers for table

% Totals

saveStr = 'results/model_output.csv';

% Delete if present
if exist(saveStr, 'file') == 2
    delete(saveStr);
end


% set up table

varNames = ["scenario", "Simulation", "totalInfections", "totalFirstInfections", "totalAdmissions", "totalDeaths"];
varTypes = ["string", "double", repmat("cell", 1, numel(varNames)-2)];
outTab = table('Size', [num_simulations*size(scenarios, 1) numel(varNames)], 'VariableTypes', varTypes, ...
    'VariableNames', varNames);

outTab.scenario = repelem(scenarios, num_simulations); % assigning scenarios
outTab.Simulation = repmat([1:num_simulations]', numel(scenarios), 1); % assigning simulation numbers

% Now loop through scenarios

row_number = 0; % starting value for row

for j = 1:numel(scenarios)

    % filenames
    filenameAll = sprintf('results/results_Uni_Filtered100_scenario%d_13-Aug-2023_fit.mat', scenarios(j));

    % Get variables -- bands
    load(filenameAll, 'epiVarsCompact');

    [newDailyCases0_95, newDailyCases1_95, newDailyCases2_95, newDailyCases3_95, newDailyCasesr_95, ...
        newDailyHosp0_95, newDailyHosp1_95, newDailyHosp2_95, newDailyHosp3_95, newDailyHospr_95, ...
        newDailyDeaths0_95, newDailyDeaths1_95, newDailyDeaths2_95, newDailyDeaths3_95, newDailyDeathsr_95, ...
        Hocc_95, ~, ~, E1_95, E2_95] = getVarsToPlot(epiVarsCompact);


    % -- cumulative infections

    val_range = 1/parBase.tE*(cumsum(E1_95) + cumsum(E2_95)); % by age
    cum_infections = reshape(val_range(ind, :, :), 16, num_simulations)'; % reshape to row = simulation, column = age group

    % -- cumulative first infections

    val_range = 1/parBase.tE*cumsum(E1_95); % by age
    cum_first_infections = reshape(val_range(ind, :, :), 16, num_simulations)'; % reshape to row = simulation, column = age group

    % -- cumulative cases

    val_range = cumsum(newDailyHosp0_95 + newDailyHosp1_95 + newDailyHosp2_95 + newDailyHosp3_95 + newDailyHospr_95);
    cum_admissions = reshape(val_range(ind, :, :), 16, num_simulations)'; % reshape to row = simulation, column = age group

    % -- cumulative deaths

    val_range = cumsum(newDailyDeaths0_95 + newDailyDeaths1_95 + newDailyDeaths2_95 + newDailyDeaths3_95 + newDailyDeathsr_95);
    cum_deaths = reshape(val_range(ind, :, :), 16, num_simulations)'; % reshape to row = simulation, column = age group
    

    % Now assign to table in loop (is there a better way?)

    for k = 1:num_simulations
        
        outTab.totalInfections{row_number+k} = cum_infections(k, :);
        outTab.totalFirstInfections{row_number+k} = cum_first_infections(k, :);
        outTab.totalAdmissions{row_number+k} = cum_admissions(k, :);
        outTab.totalDeaths{row_number+k} = cum_deaths(k, :);

    end

    row_number = row_number + num_simulations; % update row number

end


% Output totals file

writetable(outTab, saveStr);

