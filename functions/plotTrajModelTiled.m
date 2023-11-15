function [] = plotTrajModelTiled(filenameBands, filenameBestFit, ...
    dateLblResults, dateLblData, overwriteFig, dataComb)
% Function that plots a grid of plots (infections, cases, hospital
% admissions, hospital occupancy, and deaths), which include a best fit
% line, confidence bands, real data, and vertical lines representing an
% intervention date

% INPUTS:
% - filenameBands: "folder/filename.mat" for bands data
% - filenameBestFit: "folder/filename.mat" for best fit data
% - dateLblResult: datetime variable with date to which data was fitted
% - dateLblData: datetime variable with date of most recent data
% - overwriteFig: set to true to overwrite png plots in results
% - dataComb: data for plotting, as output by getAllData.m

% OUTPUTS:
% - .png figures in results, if they don't already exist or if overwriteFig
%   set to true


% Set to 1 to plot best fit lines:
plotBestFit = 0;

% ----------------- IMPORT BANDS DATA -------------------------------------
load(filenameBands, 'epiVarsCompact', 'parBase');

% Popn size time series (for 1st realisation only)
popnSize = epiVarsCompact(1).N;

% Get variables for plotting
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
    newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
    newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
    newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ~, ~, E1, ...
    E2, ~, ~] = getVarsToPlot(epiVarsCompact);

E1bands = squeeze(sum(E1, 2));
E2bands = squeeze(sum(E2, 2));

% Sum over immunity status
newDailyCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3 + newDailyCasesr;
newDailyHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3 + newDailyHospr;
newDailyDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3 + newDailyDeathsr;

% Totals across all age groups & immunity status:
incidenceRel_all = squeeze(1/parBase.tE*sum(E1+E2, 2)./sum(popnSize, 2))';       % incidence relative to popn size
newDailyCases_all = squeeze(sum(newDailyCases, 2))';
newDailyHosp_all = squeeze(sum(newDailyHosp, 2))';
hospOcc_all = squeeze(sum( Hocc, 2))';
% hospCum_all = squeeze(sum( cumsum(newDailyHosp), 2))';
newDailyDeaths_all = squeeze(sum(newDailyDeaths, 2))';

nTraj = length(epiVarsCompact);

bandsData = cat(3, incidenceRel_all(1:nTraj, :)*1e5, ...
    newDailyCases_all(1:nTraj, :), newDailyHosp_all(1:nTraj, :), ...
    hospOcc_all(1:nTraj, :), newDailyDeaths_all(1:nTraj, :));

% ----------------- IMPORT BEST FIT DATA ----------------------------------
load(filenameBestFit, 'epiVarsCompact', 'parBase', 't');
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

% plot of 1st/2nd infections
% figure 
% plot(t, sum(E1, 2), t, sum(E2, 2))
% xlim(tPlotRange)
% ylim([0, inf])
% xlabel("date")
% ylabel("infections")
% legend({'first infections', 'second infections'})
% figure
% fill([t, fliplr(t)], max(0, [min(E2bands./(E1bands+E2bands),[],2)', ...
%         fliplr(max(E2bands./(E1bands+E2bands),[],2)')]), ...
%         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% hold on
% plot(t, sum(E2, 2)./sum(E1+E2, 2))
% xline(datetime('now'), 'r--')
% hold off
% grid on
% xlim(tPlotRange)
% ylim([0, inf])
% xlabel("date")
% ylabel("proportion of new infections that are reinfections")
% legend({'', '', 'today'})

% Sum over immunity status
newDailyCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3 + newDailyCasesr;
newDailyHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3 + newDailyHospr;
newDailyDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3 + newDailyDeathsr;

% Totals across all age groups & immunity status:
incidenceRel_all = squeeze(1/parBase.tE*sum(E1+E2, 2)./sum(popnSize, 2))';       % incidence relative to popn size
newDailyCases_all = squeeze(sum(newDailyCases, 2))';
newDailyHosp_all = squeeze(sum(newDailyHosp, 2))';
hospOcc_all = squeeze(sum( Hocc, 2))';
% hospCum_all = squeeze(sum( cumsum(newDailyHosp), 2))';
newDailyDeaths_all = squeeze(sum(newDailyDeaths, 2))';

nTraj = length(epiVarsCompact);

bestFitData = cat(3, incidenceRel_all(1:nTraj, :)*1e5, ...
    newDailyCases_all(1:nTraj, :), newDailyHosp_all(1:nTraj, :), ...
    hospOcc_all(1:nTraj, :), newDailyDeaths_all(1:nTraj, :));


% ------------------- IMPORT REAL DATA ------------------------------------

tData = dataComb.date;
tDataMax = days(datetime(dateLblData) - tData(1)); %days(datetime('02-Dec-2022') - tData(1));

realData = cat(2, dataComb.NationalBorder*100/7, dataComb.nCasesData, ...
    dataComb.nHospData, dataComb.hospOccTotalMOH, dataComb.nDeathsData);

% hospOccData =  cumsum(dataComb.nHospData-dataComb.nDisc_strict);
% realData2 = cat(2, dataComb.NationalBorder*100/7, dataComb.nCasesData, ...
%     dataComb.nHospData, max(0, hospOccData), dataComb.nDeathsData);


% ------------------------------ PLOTS ------------------------------------

plotTitle = split(filenameBands, '/');
plotTitle = plotTitle(2);

subplotTitles = {'Infections', 'Cases', 'Hosp. Admissions', ...
    'Hosp. Occupancy', 'Deaths'};
subplotYaxis = {'new daily infections per 100,000', 'new daily cases', ...
    'new daily admissions for COVID', 'hospital occupancy for COVID', ...
    'daily deaths'};

smoothDays = [0, 7, 14, 7, 21];

f = figure;
f.Position = [100 100 1000 800];
tl = tiledlayout(3, 2);
title(tl, plotTitle, 'Interpreter', 'none')
subtitle(tl, 'Aggregated over all ages')

lightBLUE = [0.356862745098039,0.811764705882353,0.956862745098039];
darkBLUE = [0.0196078431372549,0.0745098039215686,0.670588235294118];
red = [1, 0, 0];
green = [0, 0.5, 0];
colorGRADIENTflexible = @(i) green + (red-green)*i;

for i = 1:size(bandsData, 3)
    nexttile
    title(subplotTitles(i))
    hold on
    fill([t, fliplr(t)], [min(bandsData(:, :, i),[],1), ...
        fliplr(max(bandsData(:, :, i),[],1))], ...
        'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    plot(t, bandsData(1:20, :, i), '-k', 'Color', '#8b8b8b')
    if plotBestFit == 1 
        plot(t, median(bestFitData(:, :, i),1), 'k-', 'LineWidth', 2)
    end
    plot(tData(1:tDataMax), realData(1:tDataMax, i), 'b.')
    if i > 1
        plot(tData(1:tDataMax), smoothdata(realData(1:tDataMax, i), 'movmean', smoothDays(i)), 'b-', 'LineWidth', 2)
    end
    
    xline(datetime(dateLblResults), 'r--');
    % Change of criteria for antiviral - start of divergence in deaths
%     xline(datetime("14-Sep-2022"), '--', 'Color', '#028A25')


    hold off
    xlim(tPlotRange)
    ylim([0 inf])
    ylabel(subplotYaxis(i))
    grid on
    grid minor
end

% leg = legend('95% confidence band', 'best fit', 'data', 'data (av.)', ...
%                 ['first intervention date ','(',char(datesScenarioChange(1)),')'], ...
%                 ['fitted to date ','(',char(dateLblResults),')']);
% leg = legend('95% confidence band', 'best fit', 'data', ...
%                 ['first intervention date ','(',char(datesScenarioChange(1)),')'], ...
%                 ['second intervention date ','(',char(datesScenarioChange(2)),')'], ...
%                 ['fitted to date ','(',char(dateLblResults),')']);
if plotBestFit == 1
    leg = legend('95% confidence band', 'example trajectories', '', '', '', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
        'best fit', 'data', 'data (av.)', ...
        ['fitted to date ','(',char(dateLblResults),')']);
else
    leg = legend('95% confidence band', 'example trajectories', '', '', '', ...
        '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ...
        'data', 'data (av.)', ...
        ['fitted to date ','(',char(dateLblResults),')']);
end
leg.Layout.Tile = 'East';


% Save figure if it doesn't already exist or if overwrite flag is on
figLabel = split(plotTitle, '.mat');
figLabel = append('results/', figLabel{1}, '.png');
if exist(figLabel, 'file') == 0 || overwriteFig
    saveas(f, figLabel)
end

end