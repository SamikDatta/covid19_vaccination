function [] = plotTrajModelAgeSplit(filenameBands, filenameBestFit, dataComb, ...
    dateLblResults, parsToFit, iscenario, overwriteFig)
% Function that plots a grid of plots (cases, hospital
% admissions, and deaths), which include a best fit
% line, confidence bands, real data, and vertical lines representing an
% intervention date for each 10year age band

% INPUTS:
% - filenameBands: "folder/filename.mat" for bands data
% - filenameBestFit: "folder/filename.mat" for best fit data
% - dataComb: data for plotting, as output by getAllData.m
% - dateLblResult: datetime variable with date to which data was fitted
% - parsToFit: list of parameter names, as defined at the top of the main file
% - iscenario: index of scenario to plot from the parBase.scenarios
% - overwriteFig: set to true to overwrite png plots in results

% OUTPUTS:
% - .png figures in results, if they don't already exist or if overwriteFig
%   set to true

% Set to 1 to plot best fit lines:
plotBestFit = 0;

% ----------------- IMPORT BANDS DATA -------------------------------------
load(filenameBands, 'epiVarsCompact', 'parBase', 'pUniFiltered');

% Popn size time series (for 1st realisation only)
popnSize = epiVarsCompact(1).N;

% Get variables for plotting
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
    newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
    newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
    newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ~, ~, E1, ...
    E2, ~, ~] = getVarsToPlot(epiVarsCompact);

% Sum over immunity status
newDailyCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3 + newDailyCasesr;
newDailyHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3 + newDailyHospr;
newDailyDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3 + newDailyDeathsr;

% Totals across all age groups & immunity status
incidenceRel_all = squeeze(1/parBase.tE*sum(E1+E2, 2)./sum(popnSize, 2))';       % incidence relative to popn size
newDailyCases_all = squeeze(sum(newDailyCases, 2))';
newDailyHosp_all = squeeze(sum(newDailyHosp, 2))';
newDailyDeaths_all = squeeze(sum(newDailyDeaths, 2))';

% Get pTest values from 95% parameter sets
ThetaTemp = array2table(pUniFiltered, 'VariableNames', parsToFit);
for thetai = 1:size(pUniFiltered, 1)
    parInd = getParUnified(ThetaTemp(thetai, :), parBase, parBase.scenarios(iscenario, :));
    pTestTS95(:, :, thetai) = parInd.pTestTS';
end

bandsData = max(0, cat(3, newDailyCases_all, ...
    newDailyHosp_all, newDailyDeaths_all, ...
    newDailyHosp_all ./ newDailyCases_all, ...
    squeeze(sum(pTestTS95.*((E1+E2)./sum(E1+E2, 2)), 2))', ...
    zeros(size(newDailyCases_all))));
bandsDataAge = max(0, cat(4, newDailyCases, newDailyHosp, newDailyDeaths, ...
    newDailyHosp ./ newDailyCases, pTestTS95, zeros(size(newDailyCases))));

% ----------------- IMPORT BEST FIT DATA ----------------------------------
load(filenameBestFit, 'epiVarsCompact', 'parBase', 't', 'pUniFiltered');
t = datetime(t,'ConvertFrom','datenum');
% tPlotRange = [t(days(datetime('01-Jan-2022') - t(1))), t(days(datetime('30-Sep-2023') - t(1)))];
tPlotRange = [t(days(datetime('01-Jan-2022') - t(1))), t(end)]; % just do until last day of simulation

% Popn size time series (for 1st realisation only)
popnSize = epiVarsCompact(1).N;

% Get variables for plotting
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
    newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
    newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
    newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ~, ~, E1, ...
    E2, ~, ~] = getVarsToPlot(epiVarsCompact);

% Sum over immunity status
newDailyCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3 + newDailyCasesr;
newDailyHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3 + newDailyHospr;
newDailyDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3 + newDailyDeathsr;

% Totals across all age groups & immunity status:
incidenceRel_all = squeeze(1/parBase.tE*sum(E1+E2, 2)./sum(popnSize, 2))';       % incidence relative to popn size
newDailyCases_all = squeeze(sum(newDailyCases, 2))';
newDailyHosp_all = squeeze(sum(newDailyHosp, 2))';
newDailyDeaths_all = squeeze(sum(newDailyDeaths, 2))';


% Get pTest values from 95% parameter sets
ThetaTemp = array2table(pUniFiltered, 'VariableNames', parsToFit);
parInd = getParUnified(ThetaTemp, parBase, parBase.scenarios(iscenario, :));
pTestTSbest = parInd.pTestTS';


bestFitData = cat(3, newDailyCases_all, ...
    newDailyHosp_all, ...
    newDailyDeaths_all, ...
    newDailyHosp_all ./ newDailyCases_all, ...
    sum(pTestTSbest.*((E1+E2)./sum(E1+E2, 2)), 2)', ...
    cumsum(sum(E1, 2), 1)', ...
    cumsum(sum(E1+E2, 2), 1)');

bestFitDataAge = cat(3, newDailyCases, newDailyHosp, newDailyDeaths, ...
    max(newDailyHosp ./ newDailyCases, 0), E1 + E2, cumsum(E1, 1), ...
    cumsum(E1+E2, 1), pTestTSbest);

% ------------------- IMPORT REAL DATA ------------------------------------

tData = dataComb.date;

realData = cat(2, dataComb.nCasesData, dataComb.nHospData, ...
    dataComb.nDeathsData, dataComb.nHospData ./ dataComb.nCasesData, ...
    zeros(size(dataComb.nCasesData)));

nCasesDataAge = dataComb.nCases_v0 + dataComb.nCases_v1 + ...
    dataComb.nCases_v2 + dataComb.nCases_v3;
nHospDataAge = dataComb.nHosp_strict_v0_byDateOfAdmission + dataComb.nHosp_strict_v1_byDateOfAdmission + ...
    dataComb.nHosp_strict_v2_byDateOfAdmission + dataComb.nHosp_strict_v3_byDateOfAdmission;
nDeathsDataAge = dataComb.nDeaths_strict_v0_byDateOfDeath + dataComb.nDeaths_strict_v1_byDateOfDeath + ...
    dataComb.nDeaths_strict_v2_byDateOfDeath + dataComb.nDeaths_strict_v3_byDateOfDeath;

hospCasesRatio = nHospDataAge ./ nCasesDataAge;
hospCasesRatio(hospCasesRatio == Inf) = NaN;

realDataAge = cat(3, nCasesDataAge, nHospDataAge, nDeathsDataAge, ...
    hospCasesRatio, zeros(size(nCasesDataAge)));



% ------------------------------ PLOTS ------------------------------------

plotTitle = split(filenameBands, '/');
plotTitle = plotTitle(2);

age_groups = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'};
subplotYaxis = {'new daily cases', 'new daily admissions', 'daily deaths', ...
    'hosp. admissions/cases', 'CAR', 'C. infections / pop.'};
smoothDays = [7, 14, 21, 14, 7];

f = figure;
f.Position = [50 0 1800 700];
tl = tiledlayout(6, 9);
% title(tl, plotTitle, 'Interpreter', 'none')

for i = 1:size(subplotYaxis, 2)
    for ag = 0:size(age_groups, 2)
        nexttile
        hold on
        if ag == 0
            title('all ages')
            fill([t, fliplr(t)], [min(bandsData(:, :, i),[], 1), ...
                fliplr(max(bandsData(:, :, i),[], 1))], ...
                'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            if i < 5
                if plotBestFit == 1 ; plot(t, median(bestFitData(:, :, i), 1), 'k-'); end
            elseif i == 5
                plot(t, median(bandsData(:, :, i), 1), 'k-')
            else
                plot(t, bestFitData(:, :, 6) ./ sum(popnSize, 2)', 'k-')
                plot(t, bestFitData(:, :, 7) ./ sum(popnSize, 2)', 'k--')
            end
            if i < 5; plot(tData, smoothdata(realData(:, i), 'movmean', 7), 'b-'); end
            
            % Remove data before first community seeds to calculate plot
            % y limits from bands:
            cropAt = datenum(parInd.dateSeed) - parBase.date0;
            if i < 6; ylim([min(bandsData(:, cropAt:end, i), [], 'all'), max(bandsData(:, cropAt:end, i), [], 'all')]); end

        else
            title(age_groups(ag))
            if i < 4
                fill([t, fliplr(t)], [min(sum(bandsDataAge(:, ag*2-1:ag*2, :, i), 2), [], 3)', ...
                    fliplr(max(sum(bandsDataAge(:, ag*2-1:ag*2, :, i), 2), [], 3)')], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                if plotBestFit == 1; plot(t, sum(bestFitDataAge(:, ag*2-1:ag*2, i), 2)', 'k-'); end
                plot(tData, smoothdata(realDataAge(:, ag, i), 'movmean', smoothDays(i)), 'b-')
            elseif i == 4    % daily hosp/cases graphs
                % Ratio of sum for current 10-year age band
                hospCasesi = sum(bandsDataAge(:, ag*2-1:ag*2, :, 2), 2) ./ sum(bandsDataAge(:, ag*2-1:ag*2, :, 1), 2);
                hospCasesi(isinf(hospCasesi) | isnan(hospCasesi)) = 0;
                fill([t, fliplr(t)], [min(hospCasesi, [], 3)', fliplr(max(hospCasesi, [], 3)')], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                if plotBestFit == 1
                    plot(t, sum(bestFitDataAge(:, ag*2-1:ag*2, 2), 2)' ./ ...
                        sum(bestFitDataAge(:, ag*2-1:ag*2, 1), 2)', 'k-')
                end
                plot(tData, smoothdata(realDataAge(:, ag, i), 'movmean', smoothDays(i)), 'b-')
                ylim([min(min(smoothdata(realDataAge(:, ag, i), 'movmean', smoothDays(i))), ...
                    min(hospCasesi(cropAt:end, :, :), [], 'all')), ...
                    max(max(smoothdata(realDataAge(:, ag, i), 'movmean', smoothDays(i))),...
                    max(hospCasesi(cropAt:end, :, :), [], 'all'))])
            elseif i == 5 % Plots of ratio of reported cases/infections
                pTestLB = max(0, (min(bandsDataAge(:, ag*2-1, :, i), [], 3) .* (E1(:, ag*2-1) + E2(:, ag*2-1)) + ...
                    min(bandsDataAge(:, ag*2, :, i), [], 3) .* (E1(:, ag*2) + E2(:, ag*2))) ./ ...
                    (E1(:, ag*2-1) + E2(:, ag*2-1) + E1(:, ag*2) + E2(:, ag*2)));
                pTestUB = max(0, (max(bandsDataAge(:, ag*2-1, :, i), [], 3) .* (E1(:, ag*2-1) + E2(:, ag*2-1)) + ...
                    max(bandsDataAge(:, ag*2, :, i), [], 3) .* (E1(:, ag*2) + E2(:, ag*2))) ./ ...
                    (E1(:, ag*2-1) + E2(:, ag*2-1) + E1(:, ag*2) + E2(:, ag*2)));
                pTestMed = max(0, (median(bandsDataAge(:, ag*2-1, :, i), 3) .* (E1(:, ag*2-1) + E2(:, ag*2-1)) + ...
                    median(bandsDataAge(:, ag*2, :, i), 3) .* (E1(:, ag*2) + E2(:, ag*2))) ./ ...
                    (E1(:, ag*2-1) + E2(:, ag*2-1) + E1(:, ag*2) + E2(:, ag*2)));
                fill([t, fliplr(t)], [pTestLB', fliplr(pTestUB')], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                plot(t, pTestMed', 'k-') % Median
            else
                plot(t, sum(bestFitDataAge(:, ag*2-1:ag*2, 6), 2) ./ ...
                    sum(popnSize(:, ag*2-1:ag*2), 2), 'k-')
                plot(t, sum(bestFitDataAge(:, ag*2-1:ag*2, 7), 2) ./ ...
                    sum(popnSize(:, ag*2-1:ag*2), 2), 'k--')
            end
        end
        hold off
        xlim(tPlotRange)
        if i == 5; ylim([0, 1]); end
        if i == 6; ylim([0, 1.5]); end
        ylabel(subplotYaxis(i))
        grid on
        grid minor
    end
end
% leg = legend('95% confidence band', 'smoothed data', ...
%         ['fitted to date ','(',char(dateLblResults),')']);
% leg.Layout.Tile = 'South';


% Save figure if it doesn't already exist or if overwrite flag is on
figLabel = split(plotTitle, '.mat');
figLabel = append('results/', figLabel{1}, '_ageSplit.png');
if exist(figLabel, 'file') == 0 || overwriteFig
    saveas(f, figLabel)
end

end