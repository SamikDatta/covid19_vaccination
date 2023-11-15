%% Analysis to compare scenarios

% Base folder to start in - make this the Gitlab repo root folder 'ode-model'

clearvars % remove all variables
close all % close all figures


% Latency parameter needed for analysis
parBase.tE = 1; % Latent period (set manually so model fit not needed)

% Set future parameters

scenarios = [1:9]';
date_of_change = '05-Mar-2021'; % Date of change from which we want comparisons


%% Code for generating necessary numbers for table

% Totals

saveStr = 'results/results_for_vaccination_scenarios_paper/compare_scenarios_20231016.xlsx';
saveStr_age = 'results/results_for_vaccination_scenarios_paper/compare_scenarios_20231016_age_split.xlsx';

% Delete if present
if exist(saveStr, 'file') == 2
    delete(saveStr);
end

% By age group

% Adding every two age groups: https://au.mathworks.com/matlabcentral/answers/521241-sum-every-n-columns-element-by-element

% Delete if present
if exist(saveStr_age, 'file') == 2
    delete(saveStr_age);
end

% set up table

varNames = ["scenario", "cumInfectionsMedian", "cumInfectionsLower", "cumInfectionsUpper", ...
    "cumInfectionsDiff", "cumInfectionsDiffLower", "cumInfectionsDiffUpper", ...
    "cumInfectionsDiffPerc", "cumInfectionsDiffPercLower", "cumInfectionsDiffPercUpper", ...
    "cumCasesMedian", "cumCasesLower", "cumCasesUpper", ...
    "cumCasesDiff", "cumCasesDiffLower", "cumCasesDiffUpper", ...
    "cumCasesDiffPerc", "cumCasesDiffPercLower", "cumCasesDiffPercUpper", ...
    "cumAdmissionsMedian", "cumAdmissionsLower", "cumAdmissionsUpper", ...
    "cumAdmissionsDiff", "cumAdmissionsDiffLower", "cumAdmissionsDiffUpper", ...
    "cumAdmissionsDiffPerc", "cumAdmissionsDiffPercLower", "cumAdmissionsDiffPercUpper", ...
    "cumDeathsMedian", "cumDeathsLower", "cumDeathsUpper", ...
    "cumDeathsDiff", "cumDeathsDiffLower", "cumDeathsDiffUpper", ...
    "cumDeathsDiffPerc", "cumDeathsDiffPercLower", "cumDeathsDiffPercUpper", ...
    "peakOccupancyMedian", "peakOccupancyLower", "peakOccupancyUpper", ...
    "peakOccupancyDiff", "peakOccupancyDiffLower", "peakOccupancyDiffUpper", ...
    "peakOccupancyDiffPerc", "peakOccupancyDiffPercLower", "peakOccupancyDiffPercUpper"];
varTypes = ["string", repmat("cell", 1, numel(varNames)-1)];
outTab = table('Size', [size(scenarios, 1) numel(varNames)], 'VariableTypes', varTypes, ...
    'VariableNames', varNames);

outTab.scenario = scenarios;

% Now loop through scenarios

row_number = 0; % start row counter
sheet_number = 0; % start sheet counter


for j = 1:numel(scenarios)

    row_number = row_number + 1; % increase counter by 1
    sheet_number = sheet_number + 1; % increase counter by 1

    % filenames
    filenameBands = sprintf('results_Uni_Filtered100_scenario%d_13-Aug-2023_fit.mat', scenarios(j));

    if j == 1 % first time round

        % Time range
        load(filenameBands, 't'); % load time vector
        ind = find(t == datenum(date_of_change)):numel(t); % for the whole time period
        weekly = 1:7:numel(ind); % getting weekly from this
        time_range = datetime(t(ind(weekly)), 'ConvertFrom','datenum');

    end

    % Get variables -- bands
    load(filenameBands, 'epiVarsCompact');

    [newDailyCases0_95, newDailyCases1_95, newDailyCases2_95, newDailyCases3_95, newDailyCasesr_95, ...
        newDailyHosp0_95, newDailyHosp1_95, newDailyHosp2_95, newDailyHosp3_95, newDailyHospr_95, ...
        newDailyDeaths0_95, newDailyDeaths1_95, newDailyDeaths2_95, newDailyDeaths3_95, newDailyDeathsr_95, ...
        Hocc_95, ~, ~, E1_95, E2_95] = getVarsToPlot(epiVarsCompact);

    % reset age table

    ageTab = []; % reset structure
    ageTab.dates = time_range';
    ageTab.scenario = repmat(scenarios(j), numel(time_range), 1);


    % -- cumulative infections

    val_range = 1/parBase.tE*(E1_95(ind, :, :) + E2_95(ind, :, :));
    val_range = reshape(cumsum(sum(val_range, 2)), size(val_range, 1), size(val_range, 3));
    slice = quantile(val_range, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
    outTab.cumInfectionsMedian{row_number} = slice(weekly, 2);
    outTab.cumInfectionsLower{row_number} = slice(weekly, 1);
    outTab.cumInfectionsUpper{row_number} = slice(weekly, 3);

    % ranges of differences
    if j == 1
        inf_zero = val_range;
    else
        inf_scen = val_range;
        diff_vec = inf_scen - inf_zero; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        %             outTab.cumInfectionsDiff{row_number} = outTab.cumInfectionsMedian{row_number} - inf_best_fit; % using best fit
        outTab.cumInfectionsDiff{row_number} = slice(weekly, 2); % using median
        outTab.cumInfectionsDiffLower{row_number} = slice(weekly, 1);
        outTab.cumInfectionsDiffUpper{row_number} = slice(weekly, 3);
        percent_vec = diff_vec*100 ./ inf_zero;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        outTab.cumInfectionsDiffPerc{row_number} = slice(weekly, 2); % using median
        outTab.cumInfectionsDiffPercLower{row_number} = slice(weekly, 1);
        outTab.cumInfectionsDiffPercUpper{row_number} = slice(weekly, 3);
    end

    % -- cumulative infections by age

    val_range = 1/parBase.tE*(E1_95(ind, :, :) + E2_95(ind, :, :));
    val_range = cumsum(squeeze(sum(reshape(val_range, size(val_range, 1), 2, 8, []), 2)));
    slice = quantile(val_range, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
    ageTab.cumInfectionsMedian = slice(weekly, :, 2);
    ageTab.cumInfectionsLower = slice(weekly, :, 1);
    ageTab.cumInfectionsUpper = slice(weekly, :, 3);

    % ranges of differences
    if j == 1
        inf_zero_age = val_range;
        ageTab.cumInfectionsDiff = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumInfectionsDiffLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumInfectionsDiffUpper = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumInfectionsDiffPerc = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumInfectionsDiffPercLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumInfectionsDiffPercUpper = NaN(size(ageTab.cumInfectionsMedian));
    else
        inf_scen = val_range;
        diff_vec = inf_scen - inf_zero_age; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        %             ageTab.cumInfectionsDiff = ageTab.cumInfections - inf_best_fit; % using best fit
        ageTab.cumInfectionsDiff = slice(weekly, :, 2); % using median
        ageTab.cumInfectionsDiffLower = slice(weekly, :, 1);
        ageTab.cumInfectionsDiffUpper = slice(weekly, :, 3);
        percent_vec = diff_vec*100 ./ inf_zero_age;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        ageTab.cumInfectionsDiffPerc = slice(weekly, :, 2); % using median
        ageTab.cumInfectionsDiffPercLower = slice(weekly, :, 1);
        ageTab.cumInfectionsDiffPercUpper = slice(weekly, :, 3);
    end

    % -- cumulative cases

    val_range = newDailyCases0_95(ind, :, :) + newDailyCases1_95(ind, :, :) + ...
        newDailyCases2_95(ind, :, :) + newDailyCases3_95(ind, :, :) + ...
        newDailyCasesr_95(ind, :, :);
    val_range = reshape(cumsum(sum(val_range, 2)), size(val_range, 1), size(val_range, 3));
    slice = quantile(val_range, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
    outTab.cumCasesMedian{row_number} = slice(weekly, 2);
    outTab.cumCasesLower{row_number} = slice(weekly, 1);
    outTab.cumCasesUpper{row_number} = slice(weekly, 3);

    % ranges of differences
    if j == 1
        cases_zero = val_range;
    else
        cases_scen = val_range;
        diff_vec = cases_scen - cases_zero; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        outTab.cumCasesDiff{row_number} = slice(weekly, 2); % using median
        outTab.cumCasesDiffLower{row_number} = slice(weekly, 1);
        outTab.cumCasesDiffUpper{row_number} = slice(weekly, 3);
        percent_vec = diff_vec*100 ./ cases_zero;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs in %
        outTab.cumCasesDiffPerc{row_number} = slice(weekly, 2); % using median
        outTab.cumCasesDiffPercLower{row_number} = slice(weekly, 1);
        outTab.cumCasesDiffPercUpper{row_number} = slice(weekly, 3);
    end

    % -- cumulative cases by age

    val_range = newDailyCases0_95(ind, :, :) + newDailyCases1_95(ind, :, :) + ...
        newDailyCases2_95(ind, :, :) + newDailyCases3_95(ind, :, :) + ...
        newDailyCasesr_95(ind, :, :);
    val_range = cumsum(squeeze(sum(reshape(val_range, size(val_range, 1), 2, 8, []), 2)));
    slice = quantile(val_range, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
    ageTab.cumCasesMedian = slice(weekly, :, 2);
    ageTab.cumCasesLower = slice(weekly, :, 1);
    ageTab.cumCasesUpper = slice(weekly, :, 3);

    % ranges of differences
    if j == 1
        cases_zero_age = val_range;
        ageTab.cumCasesDiff = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumCasesDiffLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumCasesDiffUpper = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumCasesDiffPerc = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumCasesDiffPercLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumCasesDiffPercUpper = NaN(size(ageTab.cumInfectionsMedian));
    else
        cases_scen = val_range;
        diff_vec = cases_scen - cases_zero_age; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        ageTab.cumCasesDiff = slice(weekly, :, 2); % using median
        ageTab.cumCasesDiffLower = slice(weekly, :, 1);
        ageTab.cumCasesDiffUpper = slice(weekly, :, 3);
        percent_vec = diff_vec*100 ./ cases_zero_age;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs in %
        ageTab.cumCasesDiffPerc = slice(weekly, :, 2); % using median
        ageTab.cumCasesDiffPercLower = slice(weekly, :, 1);
        ageTab.cumCasesDiffPercUpper = slice(weekly, :, 3);
    end

    % -- cumulative hospital admissions

    val_range = newDailyHosp0_95(ind, :, :) + newDailyHosp1_95(ind, :, :) + ...
        newDailyHosp2_95(ind, :, :) + newDailyHosp3_95(ind, :, :) + ...
        newDailyHospr_95(ind, :, :);
    val_range = reshape(cumsum(sum(val_range, 2)), size(val_range, 1), size(val_range, 3));
    slice = quantile(val_range, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
    outTab.cumAdmissionsMedian{row_number} = slice(weekly, 2);
    outTab.cumAdmissionsLower{row_number} = slice(weekly, 1);
    outTab.cumAdmissionsUpper{row_number} = slice(weekly, 3);

    % ranges of differences
    if j == 1
        admissions_zero = val_range;
    else
        admissions_scen = val_range;
        diff_vec = admissions_scen - admissions_zero; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        outTab.cumAdmissionsDiff{row_number} = slice(weekly, 2); % using median
        outTab.cumAdmissionsDiffLower{row_number} = slice(weekly, 1);
        outTab.cumAdmissionsDiffUpper{row_number} = slice(weekly, 3);
        percent_vec = diff_vec*100 ./ admissions_zero;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs in %
        outTab.cumAdmissionsDiffPerc{row_number} = slice(weekly, 2); % using median
        outTab.cumAdmissionsDiffPercLower{row_number} = slice(weekly, 1);
        outTab.cumAdmissionsDiffPercUpper{row_number} = slice(weekly, 3);
    end

    % -- cumulative hospital admissions by age

    val_range = newDailyHosp0_95(ind, :, :) + newDailyHosp1_95(ind, :, :) + ...
        newDailyHosp2_95(ind, :, :) + newDailyHosp3_95(ind, :, :) + ...
        newDailyHospr_95(ind, :, :);
    val_range = cumsum(squeeze(sum(reshape(val_range, size(val_range, 1), 2, 8, []), 2)));
    slice = quantile(val_range, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
    ageTab.cumAdmissionsMedian = slice(weekly, :, 2);
    ageTab.cumAdmissionsLower = slice(weekly, :, 1);
    ageTab.cumAdmissionsUpper = slice(weekly, :, 3);

    % ranges of differences
    if j == 1
        admissions_zero_age = val_range;
        ageTab.cumAdmissionsDiff = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumAdmissionsDiffLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumAdmissionsDiffUpper = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumAdmissionsDiffPerc = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumAdmissionsDiffPercLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumAdmissionsDiffPercUpper = NaN(size(ageTab.cumInfectionsMedian));
    else
        admissions_scen = val_range;
        diff_vec = admissions_scen - admissions_zero_age; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        ageTab.cumAdmissionsDiff = slice(weekly, :, 2); % using median
        ageTab.cumAdmissionsDiffLower = slice(weekly, :, 1);
        ageTab.cumAdmissionsDiffUpper = slice(weekly, :, 3);
        percent_vec = diff_vec*100 ./ admissions_zero_age;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs in %
        ageTab.cumAdmissionsDiffPerc = slice(weekly, :, 2); % using median
        ageTab.cumAdmissionsDiffPercLower = slice(weekly, :, 1);
        ageTab.cumAdmissionsDiffPercUpper = slice(weekly, :, 3);
    end

    % -- cumulative deaths

    val_range = newDailyDeaths0_95(ind, :, :) + newDailyDeaths1_95(ind, :, :) + ...
        newDailyDeaths2_95(ind, :, :) + newDailyDeaths3_95(ind, :, :) + ...
        newDailyDeathsr_95(ind, :, :);
    val_range = reshape(cumsum(sum(val_range, 2)), size(val_range, 1), size(val_range, 3));
    slice = quantile(val_range, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
    outTab.cumDeathsMedian{row_number} = slice(weekly, 2);
    outTab.cumDeathsLower{row_number} = slice(weekly, 1);
    outTab.cumDeathsUpper{row_number} = slice(weekly, 3);

    % ranges of differences
    if j == 1
        deaths_zero = val_range;
    else
        deaths_scen = val_range;
        diff_vec = deaths_scen - deaths_zero; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        outTab.cumDeathsDiff{row_number} = slice(weekly, 2); % using median
        outTab.cumDeathsDiffLower{row_number} = slice(weekly, 1);
        outTab.cumDeathsDiffUpper{row_number} = slice(weekly, 3);
        percent_vec = diff_vec*100 ./ deaths_zero;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs in %
        outTab.cumDeathsDiffPerc{row_number} = slice(weekly, 2); % using median
        outTab.cumDeathsDiffPercLower{row_number} = slice(weekly, 1);
        outTab.cumDeathsDiffPercUpper{row_number} = slice(weekly, 3);
    end

    % -- cumulative deaths by age

    val_range = newDailyDeaths0_95(ind, :, :) + newDailyDeaths1_95(ind, :, :) + ...
        newDailyDeaths2_95(ind, :, :) + newDailyDeaths3_95(ind, :, :) + ...
        newDailyDeathsr_95(ind, :, :);
    val_range = cumsum(squeeze(sum(reshape(val_range, size(val_range, 1), 2, 8, []), 2)));
    slice = quantile(val_range, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
    ageTab.cumDeathsMedian = slice(weekly, :, 2);
    ageTab.cumDeathsLower = slice(weekly, :, 1);
    ageTab.cumDeathsUpper = slice(weekly, :, 3);

    % ranges of differences
    if j == 1
        deaths_zero_age = val_range;
        ageTab.cumDeathsDiff = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumDeathsDiffLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumDeathsDiffUpper = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumDeathsDiffPerc = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumDeathsDiffPercLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.cumDeathsDiffPercUpper = NaN(size(ageTab.cumInfectionsMedian));
    else
        deaths_scen = val_range;
        diff_vec = deaths_scen - deaths_zero_age; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        ageTab.cumDeathsDiff = slice(weekly, :, 2); % using median
        ageTab.cumDeathsDiffLower = slice(weekly, :, 1);
        ageTab.cumDeathsDiffUpper = slice(weekly, :, 3);
        percent_vec = diff_vec*100 ./ deaths_zero_age;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs in %
        ageTab.cumDeathsDiffPerc = slice(weekly, :, 2); % using median
        ageTab.cumDeathsDiffPercLower = slice(weekly, :, 1);
        ageTab.cumDeathsDiffPercUpper = slice(weekly, :, 3);
    end

    % -- Peak hospital occupancy (use cumulative maximum)

    val_range = sum(Hocc_95(ind, :, :), 2); % sum across age groups
    val_range = reshape(val_range, size(val_range, 1), size(val_range, 3)); % reshape
    val_range = cummax(val_range); % get cumulative maximum for each run
    slice = quantile(val_range, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
    outTab.peakOccupancyMedian{row_number} = slice(weekly, 2);
    outTab.peakOccupancyLower{row_number} = slice(weekly, 1);
    outTab.peakOccupancyUpper{row_number} = slice(weekly, 3);

    % ranges of differences
    if j == 1
        occupancy_zero = val_range;
    else
        occupancy_scen = val_range;
        diff_vec = occupancy_scen - occupancy_zero; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs
        outTab.peakOccupancyDiff{row_number} = slice(weekly, 2); % using median
        outTab.peakOccupancyDiffLower{row_number} = slice(weekly, 1);
        outTab.peakOccupancyDiffUpper{row_number} = slice(weekly, 3);
        percent_vec = diff_vec*100 ./ occupancy_zero;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 2); % getting median and 95% CIs in %
        outTab.peakOccupancyDiffPerc{row_number} = slice(weekly, 2); % using median
        outTab.peakOccupancyDiffPercLower{row_number} = slice(weekly, 1);
        outTab.peakOccupancyDiffPercUpper{row_number} = slice(weekly, 3);
    end

    % -- Peak hospital occupancy by age

    % ageTab.peakOccupancy = cummax(Hocc(ind, :));
    val_range = Hocc_95(ind, :, :); % cumulative max across time by age group
    val_range = cummax(squeeze(sum(reshape(val_range, size(val_range, 1), 2, 8, []), 2)));
    slice = quantile(val_range, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
    ageTab.peakOccupancyMedian = slice(weekly, :, 2);
    ageTab.peakOccupancyLower = slice(weekly, :, 1);
    ageTab.peakOccupancyUpper = slice(weekly, :, 3);

    % ranges of differences
    if j == 1
        occupancy_zero_age = val_range;
        ageTab.peakOccupancyDiff = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.peakOccupancyDiffLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.peakOccupancyDiffUpper = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.peakOccupancyDiffPerc = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.peakOccupancyDiffPercLower = NaN(size(ageTab.cumInfectionsMedian));
        ageTab.peakOccupancyDiffPercUpper = NaN(size(ageTab.cumInfectionsMedian));
    else
        occupancy_scen = val_range;
        diff_vec = occupancy_scen - occupancy_zero_age; % difference with no change in behaviour
        slice = quantile(diff_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs
        ageTab.peakOccupancyDiff = slice(weekly, :, 2); % using median
        ageTab.peakOccupancyDiffLower = slice(weekly, :, 1);
        ageTab.peakOccupancyDiffUpper = slice(weekly, :, 3);
        percent_vec = diff_vec*100 ./ occupancy_zero_age;
        slice = quantile(percent_vec, [0.0275 0.5 0.975], 3); % getting median and 95% CIs in %
        ageTab.peakOccupancyDiffPerc = slice(weekly, :, 2); % using median
        ageTab.peakOccupancyDiffPercLower = slice(weekly, :, 1);
        ageTab.peakOccupancyDiffPercUpper = slice(weekly, :, 3);
    end

    % Output to separate sheet for each scenario

    writetable(struct2table(ageTab), saveStr_age, 'sheet', sheet_number);

end


% Output totals file

writetable(outTab, saveStr);

