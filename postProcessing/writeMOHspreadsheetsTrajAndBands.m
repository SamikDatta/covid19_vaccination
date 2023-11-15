% Code for spitting out TAS spreadsheets

folderIn = '../results/results_18JunFit_06JunVax_06JunAV_26JunRun';
folderOut = '../results/test';

% Produce folder if needed
if ~exist(folderOut, 'dir')
    mkdir(folderOut);
end

dateLblResults = "18-Jun-2023";            % label for data files - data up to this date was used for fitting
cropDate = "30-Sep-2023";

scenario_labels_fIn = {'scenario1', 'scenario2', 'scenario3'}; % scenarios
scenario_labels_fOut = {'noSeasonality', '10%seasonalityEffect', '20%seasonalityEffect'}; % scenarios to spit out

% Latency parameter needed for analysis
parBase.tE = 1; % Latent period (set manually so model fit not needed)

% bands
bandsFilenameForm = 'Uni_Filtered95_';

for j = 1:length(scenario_labels_fIn)
    
    
    % Get variables -- bands
    fIn = sprintf('%s/results_Uni_Filtered95_%s_%s_fit.mat', folderIn, scenario_labels_fIn{j}, dateLblResults);
    load(fIn, 't', 'epiVarsCompact');
    
    % Get variables for plotting
    [newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
        newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
        newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
        newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ~, ~, E1, ...
        E2, ~, ~] = getVarsToPlot(epiVarsCompact);
    
    
    date0 = datenum('20JAN2022');
    ind = t >= date0 & t <= datenum(cropDate);
    
    % Save time series of main variables for MOH
    
    % -- all infections
    outTab = [];
    outTab.date = datetime(t(ind)', 'ConvertFrom', 'datenum' );
    % best then upper then lower for each
    popnSize = epiVarsCompact(1).N;
    outTab.bestInfections = 1/parBase.tE*(E1(ind, :) + E2(ind, :))./sum(popnSize(ind, :), 2);
    outTab.upper95Infections = 1/parBase.tE*(max(E1_95(ind, :, :),[],3) + max(E2_95(ind, :, :),[],3))./sum(popnSize(ind, :), 2);
    outTab.lower95Infections = 1/parBase.tE*(min(E1_95(ind, :, :),[],3) + min(E2_95(ind, :, :),[],3))./sum(popnSize(ind, :), 2);
    
    outTab = struct2table(outTab);
    fOut = sprintf('%s/infections_with_bands_%s_%s_fit.xlsx', folderOut, scenario_labels_fOut{j}, dateLblResults);
    writetable(outTab, fOut);
    
    % -- cases
    outTab = [];
    outTab.date = datetime(t(ind)', 'ConvertFrom', 'datenum' );
    % best then upper then lower for each
    outTab.bestNewDailyCases0 = newDailyCases0(ind, :);
    outTab.upper95NewDailyCases0 = max(newDailyCases0_95(ind, :, :),[],3);
    outTab.lower95NewDailyCases0 = min(newDailyCases0_95(ind, :, :),[],3);
    
    outTab.bestNewDailyCases1 = newDailyCases1(ind, :);
    outTab.upper95NewDailyCases1 = max(newDailyCases1_95(ind, :, :),[],3);
    outTab.lower95NewDailyCases1 = min(newDailyCases1_95(ind, :, :),[],3);
    
    outTab.bestNewDailyCases2 = newDailyCases2(ind, :);
    outTab.upper95NewDailyCases2 = max(newDailyCases2_95(ind, :, :),[],3);
    outTab.lower95NewDailyCases2 = min(newDailyCases2_95(ind, :, :),[],3);
    
    outTab.bestNewDailyCases3 = newDailyCases3(ind, :);
    outTab.upper95NewDailyCases3 = max(newDailyCases3_95(ind, :, :),[],3);
    outTab.lower95NewDailyCases3 = min(newDailyCases3_95(ind, :, :),[],3);
    
    outTab.bestNewDailyCasesr = newDailyCasesr(ind, :);
    outTab.upper95NewDailyCasesr = max(newDailyCasesr_95(ind, :, :),[],3);
    outTab.lower95NewDailyCasesr = min(newDailyCasesr_95(ind, :, :),[],3);
    
    outTab = struct2table(outTab);
    fOut = sprintf('%s/cases_with_bands_%s_%s_fit.xlsx', folderOut, scenario_labels_fOut{j}, dateLblResults);
    writetable(outTab, fOut);
    
    % -- new hosp
    outTab = [];
    outTab.date = datetime(t(ind)', 'ConvertFrom', 'datenum' );
    % best then upper then lower for each
    outTab.bestNewDailyHosp0 = newDailyHosp0(ind, :);
    outTab.upper95NewDailyHosp0 = max(newDailyHosp0_95(ind, :, :),[],3);
    outTab.lower95NewDailyHosp0 = min(newDailyHosp0_95(ind, :, :),[],3);
    
    outTab.bestNewDailyHosp1 = newDailyHosp1(ind, :);
    outTab.upper95NewDailyHosp1 = max(newDailyHosp1_95(ind, :, :),[],3);
    outTab.lower95NewDailyHosp1 = min(newDailyHosp1_95(ind, :, :),[],3);
    
    outTab.bestNewDailyHosp2 = newDailyHosp2(ind, :);
    outTab.upper95NewDailyHosp2 = max(newDailyHosp2_95(ind, :, :),[],3);
    outTab.lower95NewDailyHosp2 = min(newDailyHosp2_95(ind, :, :),[],3);
    
    outTab.bestNewDailyHosp3 = newDailyHosp3(ind, :);
    outTab.upper95NewDailyHosp3 = max(newDailyHosp3_95(ind, :, :),[],3);
    outTab.lower95NewDailyHosp3 = min(newDailyHosp3_95(ind, :, :),[],3);
    
    outTab.bestNewDailyHospr = newDailyHospr(ind, :);
    outTab.upper95NewDailyHospr = max(newDailyHospr_95(ind, :, :),[],3);
    outTab.lower95NewDailyHospr = min(newDailyHospr_95(ind, :, :),[],3);
    
    outTab = struct2table(outTab);
    fOut = sprintf('%s/newAdmissions_with_bands_%s_%s_fit.xlsx', folderOut, scenario_labels_fOut{j}, dateLblResults);
    writetable(outTab, fOut);
    
    % -- hosp occ
    outTab = [];
    outTab.date = datetime(t(ind)', 'ConvertFrom', 'datenum' );
    % best then upper then lower for each
    outTab.bestHospOccupancy = Hocc(ind, :);
    outTab.upper95HospOccupancy = max(Hocc_95(ind, :, :),[],3);
    outTab.lower95HospOccupancy = min(Hocc_95(ind, :, :),[],3);
    
    outTab = struct2table(outTab);
    fOut = sprintf('%s/hospOcc_with_bands_%s_%s_fit.xlsx', folderOut, scenario_labels_fOut{j}, dateLblResults);
    writetable(outTab, fOut);
    
    % -- daily death
    outTab = [];
    outTab.date = datetime(t(ind)', 'ConvertFrom', 'datenum' );
    % best then upper then lower for each
    outTab.bestNewDailyDeaths0 = newDailyDeaths0(ind, :);
    outTab.upper95NewDailyDeaths0 = max(newDailyDeaths0_95(ind, :, :),[],3);
    outTab.lower95NewDailyDeaths0 = min(newDailyDeaths0_95(ind, :, :),[],3);
    
    outTab.bestNewDailyDeaths1 = newDailyDeaths1(ind, :);
    outTab.upper95NewDailyDeaths1 = max(newDailyDeaths1_95(ind, :, :),[],3);
    outTab.lower95NewDailyDeaths1 = min(newDailyDeaths1_95(ind, :, :),[],3);
    
    outTab.bestNewDailyDeaths2 = newDailyDeaths2(ind, :);
    outTab.upper95NewDailyDeaths2 = max(newDailyDeaths2_95(ind, :, :),[],3);
    outTab.lower95NewDailyDeaths2 = min(newDailyDeaths2_95(ind, :, :),[],3);
    
    outTab.bestNewDailyDeaths3 = newDailyDeaths3(ind, :);
    outTab.upper95NewDailyDeaths3 = max(newDailyDeaths3_95(ind, :, :),[],3);
    outTab.lower95NewDailyDeaths3 = min(newDailyDeaths3_95(ind, :, :),[],3);
    
    outTab.bestNewDailyDeathsr = newDailyDeathsr(ind, :);
    outTab.upper95NewDailyDeathsr = max(newDailyDeathsr_95(ind, :, :),[],3);
    outTab.lower95NewDailyDeathsr = min(newDailyDeathsr_95(ind, :, :),[],3);
    
    outTab = struct2table(outTab);
    fOut = sprintf('%s/deaths_with_bands_%s_%s_fit.xlsx', folderOut, scenario_labels_fOut{j}, dateLblResults);
    writetable(outTab, fOut);
    
end
