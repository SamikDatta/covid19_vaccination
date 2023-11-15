function [nObs, resid, dist] = calcErrorWithAgeCurrent(t, epiVarsCompact, dataComb, par, dataVars)

% calculate error function for model output (epiVarsCompact) and data
% (dataComb)
%minCaseFitDate = datenum('01MAR2022');  % only fit to cases after 1 march
minCaseFitDate = datenum('14FEB2022');  % only fit to cases after 14 feb OJM
minOtherFitDate = datenum('01FEB2022'); % only fit to other data after 1 feb
deathLag = 10;                      % exclude this number of data points at the end of the time series for deaths to account for date of death to report lag
hospLag = 40;                       % exclude this number of data points at the end of the time series for new admissions to account for report lag
latestData = datenum(dataComb.date( find(~isnan(dataComb.nHospData), 1, 'last') ));

% Get variables for comparing to data
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
    newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
    newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
    newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ~, ~, E1, ...
    E2] = getVarsToPlot(epiVarsCompact);

% Sum over immunity status
newDailyCases = newDailyCases0+newDailyCases1+newDailyCases2+newDailyCases3+newDailyCasesr;
newDailyHosp = newDailyHosp0+newDailyHosp1+newDailyHosp2+newDailyHosp3+newDailyHospr;
newDailyDeaths = newDailyDeaths0+newDailyDeaths1+newDailyDeaths2+newDailyDeaths3+newDailyDeathsr;

% Totals across all age groups:
incidenceRel_all = (1/par.tE * sum(E1+E2, 2)./sum(epiVarsCompact.N, 2));
newDailyCases_all = sum(newDailyCases, 2);
newDailyHosp_all = sum(newDailyHosp, 2);
hospOcc_all = sum(Hocc, 2);
newDailyDeaths_all = sum(newDailyDeaths, 2);

%vars are 
% daily cases
% daily deaths
% daily new admissions
% hospital occupancy
% daily incidence per capita
% age breakdown of cases

% Mike small eps:
% sml_incid = 5e-6;
% sml_cases = 10;
% sml_hospAdm = 0.1;
% sml_hospOcc = 0.5;
% sml_deaths = 0.01;
% sml_age_frac = 5e-5;
%
% OJM: just use 1 for all

nDist = length(dataVars);
% eps = ones(1,nDist); % assuming all are 1
eps = [1, 1, 1, 0.0001, 0.01, 0.01]; % sensible values
nObs = zeros(1, nDist);
dist  = zeros(1, nDist);

for i = 1:length(dataVars)
    dataVars(i);
    if strcmp(dataVars(i),'nCasesData')
        %model is in data and model is beyond min data
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minCaseFitDate ; 
        y1 = newDailyCases_all(ind1);
        %data is in model and data is beyond min data
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minCaseFitDate;
        dataRaw = table2array(dataComb(:,dataVars(i)));
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'nDeathsData')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate & t <= latestData-deathLag;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate & datenum(dataComb.date) <= latestData-deathLag;
        y1 = newDailyDeaths_all(ind1);
        %
        dataRaw = table2array(dataComb(:, dataVars(i)));
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'nHospData')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate & t <= latestData-hospLag;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate & datenum(dataComb.date) <= latestData-hospLag;
        y1 = newDailyHosp_all(ind1);
        dataRaw = table2array(dataComb(:, dataVars(i)));
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'Total_hosp_MOHweb')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate;
        y1 = hospOcc_all(ind1);
        dataRaw = table2array(dataComb(:, dataVars(i)));
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'hospOccTotalMOH') % same as previous but with updated MOH hosp.occ. data
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate;
        y1 = hospOcc_all(ind1);
        dataRaw = table2array(dataComb(:, dataVars(i)));
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'NationalBorder')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minCaseFitDate;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minCaseFitDate;
        y1 = incidenceRel_all(ind1);
        dataRaw = table2array(dataComb(:, dataVars(i)));
        dataRaw = dataRaw/1000/7;
        y2 = dataRaw(ind2);
    elseif strcmp(dataVars(i),'CasesByAge')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minCaseFitDate;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minCaseFitDate;
        newDailyCases10 = combineBands(newDailyCases);
        y1 = newDailyCases10(ind1, :);
        % Fit to fraction of cases in each age group rather than absolute numbers:
        y1 = y1./sum(y1, 2);
        dataRaw = dataComb.nCases_v0+dataComb.nCases_v1+dataComb.nCases_v2+dataComb.nCases_v3;
        dataRaw = dataRaw(ind2, :);
        y2 = dataRaw./sum(dataRaw, 2);
%         % Just fit to fraction 60+
%         y1 = sum(y1(:, 7:8) ,2)./sum(y1, 2);
%         dataRaw = dataComb.nCases_v0+dataComb.nCases_v1+dataComb.nCases_v2+dataComb.nCases_v3;
%         dataRaw = dataRaw(ind2,:);
%         y2 = sum(dataRaw(:, 7:8) ,2)./sum(dataRaw, 2);
    elseif strcmp(dataVars(i),'HospByAge')
        ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate & t <= latestData-hospLag;
        ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate & datenum(dataComb.date) <= latestData-hospLag;
        
        newDailyHosp10 = combineBands(newDailyHosp);
        y1 = newDailyHosp10(ind1, :);
        % Fit to fraction of admissions in each age group rather than absolute numbers:
        y1 = y1./sum(y1, 2);
        dataRaw = dataComb.nHosp_strict_v0_byDateOfAdmission + ...
            dataComb.nHosp_strict_v1_byDateOfAdmission + ...
            dataComb.nHosp_strict_v2_byDateOfAdmission + ...
            dataComb.nHosp_strict_v3_byDateOfAdmission;
        dataRaw = dataRaw(ind2, :);
        y2 = dataRaw./sum(dataRaw, 2);
    else
        disp('no model output identified')
        return
    end
    
    % Check that filtered timestamps for raw data and model data are equal
    modeltFiltered = t(ind1)';
    datatFiltered = datenum(dataComb.date(ind2));
    assert(isequal(modeltFiltered, datatFiltered))

    % drop nans
    keep = ~isnan(y2) == 1;
    y2 = y2(keep);
    y1 = y1(keep);

    smoothWindow = 7;
    dataMean = log(smoothdata(y2, 'movmean', smoothWindow)+eps(i));
    modelMean = log(y1 + eps(i));

    nObs(i) = length(y2);

    resid{i} = dataMean-modelMean; %useful for diagnostics. cell array coz different number data points per metric
    dist(i) = (1/nObs(i))*sum((dataMean-modelMean).^2); % divide by nObs to normalise
end

