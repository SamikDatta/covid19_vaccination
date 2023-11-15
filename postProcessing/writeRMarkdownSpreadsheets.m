%% Pulling out info for spreadsheets to read into R Markdown documents

% Base folder to start in - make this the Gitlab repo root folder 'ode-model'

clearvars % remove all variables
close all % close all figures

% Base folder name
filename_base = sprintf('results_for_vaccination_scenarios_paper/'); % edit this as required

% Latency parameter needed for analysis
parBase.tE = 1; % Latent period (set manually so model fit not needed)

% Set future parameters

scenarios = num2str([1:9]'); % scenarios
first_date = '05-Mar-2021'; % Date to save from
last_date = '06-August-2023'; % date to save to


% Which example runs to output
% run_pick = randsample(size(E1_95, 3), 10); % 10 random runs
run_pick = [11:18, 32, 75]; % specify manually


%% Code for producing Excel spreadsheets (with bands)

for j = 1:numel(scenarios)

    % Produce folder if needed
    scenario_str = sprintf('results/%s/%s', filename_base, scenarios(j));
    if ~exist(scenario_str , 'dir')
        mkdir(scenario_str);
    end

    % filenames
    filenameBands = sprintf('results_Uni_Filtered95_scenario%s_13-Aug-2023_fit.mat', scenarios(j));
    filenameBestFit = sprintf('results_Uni_FilteredBest_scenario%s_13-Aug-2023_fit.mat', scenarios(j));

    if j == 1 % first time round
        load(filenameBestFit, 't'); % load time vector
        % ind = find(t == datenum(first_date)):find(t == datenum(last_date)); % time range
        ind = 1:numel(t); % use every timestep
        basic_table = table; % construct basic table
        basic_table.date = datetime(t(ind)', 'ConvertFrom', 'datenum'); % add dates (same for all scenarios)
    end

    % Get variables -- best fit
    load(filenameBestFit, 'epiVarsCompact');

    [newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, newDailyCasesr, ...
        newDailyHosp0, newDailyHosp1, newDailyHosp2, newDailyHosp3, newDailyHospr, ...
        newDailyDeaths0, newDailyDeaths1, newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, ...
        Hocc, S, Sw, E1, E2] = getVarsToPlot(epiVarsCompact);

    % Get variables -- bands
    load(filenameBands, 'epiVarsCompact');

    % Get variables for plotting
    [newDailyCases0_95, newDailyCases1_95, newDailyCases2_95, newDailyCases3_95, newDailyCasesr_95, ...
        newDailyHosp0_95, newDailyHosp1_95, newDailyHosp2_95, newDailyHosp3_95, newDailyHospr_95, ...
        newDailyDeaths0_95, newDailyDeaths1_95, newDailyDeaths2_95, newDailyDeaths3_95, newDailyDeathsr_95, ...
        Hocc_95, S_95, Sw_95, E1_95, E2_95] = getVarsToPlot(epiVarsCompact);

    % Save time series of main variables for MOH

    % -- population sizes over time

    outTab = basic_table;
    popnSize = epiVarsCompact(1).N;
    outTab.popnSize = popnSize(ind, :);
    tot_popn_size = sum(popnSize(ind, :), 2); % needed for other tables

    fOut = sprintf('%s/popn_sizes.xlsx', scenario_str);
    writetable(outTab, fOut);

    % -- all infections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestInfections = 1/parBase.tE*(E1(ind, :) + E2(ind, :)); % total daily infections
    infections_by_age = 1/parBase.tE*(E1_95(ind, :, :) + E2_95(ind, :, :)); % by age
    infections_total = squeeze(sum(infections_by_age, 2)); % sum along ages
    outTab.upper95TotalInfections = max(infections_total, [], 2);
    outTab.lower95TotalInfections = min(infections_total, [], 2);
    outTab.upper95Infections = max(infections_by_age, [], 3);
    outTab.lower95Infections = min(infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(E1_95(ind, :, run_pick(k)) + E2_95(ind, :, run_pick(k)));
    end

    fOut = sprintf('%s/all_infections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- cumulative all infections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestCumInfections = 1/parBase.tE*(cumsum(E1(ind, :)) + cumsum(E2(ind, :))); % get cumulative infections
    cum_infections_by_age = 1/parBase.tE*(cumsum(E1_95(ind, :, :)) + cumsum(E2_95(ind, :, :))); % by age
    cum_infections_total = squeeze(sum(cum_infections_by_age, 2)); % sum along ages
    outTab.upper95CumTotalInfections = max(cum_infections_total, [], 2);
    outTab.lower95CumTotalInfections = min(cum_infections_total, [], 2);
    outTab.upper95CumInfections = max(cum_infections_by_age, [], 3);
    outTab.lower95CumInfections = min(cum_infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(cumsum(E1_95(ind, :, run_pick(k))) + cumsum(E2_95(ind, :, run_pick(k))));
    end

    fOut = sprintf('%s/cum_all_infections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- first infections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestInfections = 1/parBase.tE*(E1(ind, :));
    infections_by_age = 1/parBase.tE*E1_95(ind, :, :); % by age
    infections_total = squeeze(sum(infections_by_age, 2)); % sum along ages
    outTab.upper95TotalInfections = max(infections_total, [], 2);
    outTab.lower95TotalInfections = min(infections_total, [], 2);
    outTab.upper95Infections = max(infections_by_age, [], 3);
    outTab.lower95Infections = min(infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(E1_95(ind, :, run_pick(k)));
    end

    fOut = sprintf('%s/first_infections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- cumulative first infections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestCumInfections = 1/parBase.tE*(cumsum(E1(ind, :))); % get cumulative infections
    cum_infections_by_age = 1/parBase.tE*cumsum(E1_95(ind, :, :)); % by age
    cum_infections_total = squeeze(sum(cum_infections_by_age, 2)); % sum along ages
    outTab.upper95CumTotalInfections = max(cum_infections_total, [], 2);
    outTab.lower95CumTotalInfections = min(cum_infections_total, [], 2);
    outTab.upper95CumInfections = max(cum_infections_by_age, [], 3);
    outTab.lower95CumInfections = min(cum_infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(cumsum(E1_95(ind, :, run_pick(k))));
    end

    fOut = sprintf('%s/cum_first_infections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- reinfections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestInfections = 1/parBase.tE*(E2(ind, :));
    infections_by_age = 1/parBase.tE*E2_95(ind, :, :); % by age
    infections_total = squeeze(sum(infections_by_age, 2)); % sum along ages
    outTab.upper95TotalInfections = max(infections_total, [], 2);
    outTab.lower95TotalInfections = min(infections_total, [], 2);
    outTab.upper95Infections = max(infections_by_age, [], 3);
    outTab.lower95Infections = min(infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(E2_95(ind, :, run_pick(k)));
    end

    fOut = sprintf('%s/reinfections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- cumulative reinfections

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestCumInfections = 1/parBase.tE*(cumsum(E2(ind, :))); % get cumulative infections
    cum_infections_by_age = 1/parBase.tE*cumsum(E2_95(ind, :, :)); % by age
    cum_infections_total = squeeze(sum(cum_infections_by_age, 2)); % sum along ages
    outTab.upper95CumTotalInfections = max(cum_infections_total, [], 2);
    outTab.lower95CumTotalInfections = min(cum_infections_total, [], 2);
    outTab.upper95CumInfections = max(cum_infections_by_age, [], 3);
    outTab.lower95CumInfections = min(cum_infections_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = 1/parBase.tE*(cumsum(E2_95(ind, :, run_pick(k))));
    end

    fOut = sprintf('%s/cum_reinfections_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);



    % -- weighted susceptibility to infection

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestWeightedSus = Sw;
    weighted_sus_by_age = Sw_95; % by age
    weighted_sus_total = squeeze(sum(weighted_sus_by_age, 2)); % sum along ages
    outTab.upper95TotalWeightedSus = max(weighted_sus_total, [], 2);
    outTab.lower95TotalWeightedSus = min(weighted_sus_total, [], 2);
    outTab.upper95NewDailyWeightedSus = max(weighted_sus_by_age, [], 3);
    outTab.lower95NewDailyWeightedSus = min(weighted_sus_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = Sw_95(ind, :, run_pick(k));
    end

    fOut = sprintf('%s/weighted_sus_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- cases

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestNewDailyCases = newDailyCases0(ind, :) + newDailyCases1(ind, :) + ...
        newDailyCases2(ind, :) + newDailyCases3(ind, :) + newDailyCasesr(ind, :);
    cases_by_age = newDailyCases0_95(ind, :, :) + newDailyCases1_95(ind, :, :) + ...
        newDailyCases2_95(ind, :, :) + newDailyCases3_95(ind, :, :) + newDailyCasesr_95(ind, :, :); % by age
    cases_total = squeeze(sum(cases_by_age, 2)); % sum along ages
    outTab.upper95TotalCases = max(cases_total, [], 2);
    outTab.lower95TotalCases = min(cases_total, [], 2);
    outTab.upper95NewDailyCases = max(cases_by_age, [], 3);
    outTab.lower95NewDailyCases = min(cases_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = newDailyCases0_95(ind, :, run_pick(k)) + newDailyCases1_95(ind, :, run_pick(k)) + ...
            newDailyCases2_95(ind, :, run_pick(k)) + newDailyCases3_95(ind, :, run_pick(k)) + ...
            newDailyCasesr_95(ind, :, run_pick(k));
    end

    fOut = sprintf('%s/cases_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- new admissions

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestNewDailyHosp = newDailyHosp0(ind, :) + newDailyHosp1(ind, :) + ...
        newDailyHosp2(ind, :) + newDailyHosp3(ind, :) + newDailyHospr(ind, :);
    admissions_by_age = newDailyHosp0_95(ind, :, :) + ...
        newDailyHosp1_95(ind, :, :) + newDailyHosp2_95(ind, :, :) + ...
        newDailyHosp3_95(ind, :, :) + newDailyHospr_95(ind, :, :); % by age
    admissions_total = squeeze(sum(admissions_by_age, 2)); % sum along ages
    outTab.upper95TotalHosp = max(admissions_total, [], 2);
    outTab.lower95TotalHosp = min(admissions_total, [], 2);
    outTab.upper95NewDailyHosp = max(admissions_by_age, [], 3);
    outTab.lower95NewDailyHosp = min(admissions_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = newDailyHosp0_95(ind, :, run_pick(k)) + ...
            newDailyHosp1_95(ind, :, run_pick(k)) + newDailyHosp2_95(ind, :, run_pick(k)) + ...
            newDailyHosp3_95(ind, :, run_pick(k)) + newDailyHospr_95(ind, :, run_pick(k));
    end

    fOut = sprintf('%s/newAdmissions_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- hospital occupancy

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestHospOccupancy = Hocc(ind, :);
    occupancy_by_age = Hocc_95(ind, :, :); % by age
    occupancy_total = squeeze(sum(occupancy_by_age, 2)); % sum along ages
    outTab.upper95TotalHospOccupancy = max(occupancy_total, [], 2);
    outTab.lower95TotalHospOccupancy = min(occupancy_total, [], 2);
    outTab.upper95HospOccupancy = max(occupancy_by_age, [], 3);
    outTab.lower95HospOccupancy = min(occupancy_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = Hocc_95(ind, :, run_pick(k));
    end

    fOut = sprintf('%s/hospOcc_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);


    % -- daily deaths

    outTab = basic_table;
    % best then upper then lower for each
    outTab.bestNewDailyDeaths = newDailyDeaths0(ind, :) + newDailyDeaths1(ind, :) + ...
        newDailyDeaths2(ind, :) + newDailyDeaths3(ind, :) + newDailyDeathsr(ind, :);
    deaths_by_age = newDailyDeaths0_95(ind, :, :) + ...
        newDailyDeaths1_95(ind, :, :) + newDailyDeaths2_95(ind, :, :) + ...
        newDailyDeaths3_95(ind, :, :) + newDailyDeathsr_95(ind, :, :); % by age
    deaths_total = squeeze(sum(deaths_by_age, 2)); % sum along ages
    outTab.upper95TotalDeaths = max(deaths_total, [], 2);
    outTab.lower95TotalDeaths = min(deaths_total, [], 2);
    outTab.upper95NewDailyDeaths = max(deaths_by_age, [], 3);
    outTab.lower95NewDailyDeaths = min(deaths_by_age, [], 3);

    for k = 1:length(run_pick)
        F0 = sprintf('run_%02d', k); % name run
        outTab.(F0) = newDailyDeaths0_95(ind, :, run_pick(k)) + ...
            newDailyDeaths1_95(ind, :, run_pick(k)) + newDailyDeaths2_95(ind, :, run_pick(k)) + ...
            newDailyDeaths3_95(ind, :, run_pick(k)) + newDailyDeathsr_95(ind, :, run_pick(k));
    end

    fOut = sprintf('%s/deaths_with_bands.xlsx', scenario_str);
    writetable(outTab, fOut);

end

