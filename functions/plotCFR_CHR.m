function plotCFR_CHR(epiVarsCompact, dataComb)


dateRange = [ datetime(2022, 1, 25), datetime(2023, 6, 30) ];



% Calclate model CHR and CFR
nReps = length(epiVarsCompact);

[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, S, Sw, E1, E2, I, A] = getVarsToPlot(epiVarsCompact);

nCases = newDailyCases0+newDailyCases1+newDailyCases2+newDailyCases3+newDailyCasesr;
nHosp =  newDailyHosp0+newDailyHosp1+newDailyHosp2+newDailyHosp3+newDailyHospr;
nDeaths = newDailyDeaths0+newDailyDeaths1+newDailyDeaths2+newDailyDeaths3+newDailyDeathsr;

inRangeFlag = t >= datenum(dateRange(1)) & t <= datenum(dateRange(2));
aggCases = squeeze( sum(nCases(inRangeFlag, :, :), 1 ));
aggHosp = squeeze( sum(nHosp(inRangeFlag, :, :), 1 ));
aggDeaths = squeeze( sum(nDeaths(inRangeFlag, :, :), 1));

coarseCases = aggCases(1:2:15, :)+aggCases(2:2:16, :);
coarseHosp = aggHosp(1:2:15, :)+aggHosp(2:2:16, :);
coarseDeaths = aggDeaths(1:2:15, :)+aggDeaths(2:2:16, :);

CHR = coarseHosp./coarseCases;
CFR = coarseDeaths./coarseCases;

qt = [0.025, 0.5, 0.975];
Q = quantile(CHR, qt, 2);
CHR_low = Q(:, 1)';
CHR_med = Q(:, 2)';
CHR_hi = Q(:, 3)';
Q = quantile(CFR, qt, 2);
CFR_low = Q(:, 1)';
CFR_med = Q(:, 2)';
CFR_hi = Q(:, 3)';



% Calculate data CHR and CFR

inRangeFlag = dataComb.date >= dateRange(1) & dataComb.date <= dateRange(2);
dataInRange = dataComb(inRangeFlag, :);
aggCasesData = sum(dataInRange.nCases_v0+dataInRange.nCases_v1+dataInRange.nCases_v2+dataInRange.nCases_v3);
aggHospData = sum(dataInRange.nHosp_strict_v0_byDateOfAdmission+dataInRange.nHosp_strict_v1_byDateOfAdmission+dataInRange.nHosp_strict_v2_byDateOfAdmission+dataInRange.nHosp_strict_v3_byDateOfAdmission);
aggDeathsData = sum(dataInRange.nDeaths_strict_v0_byDateOfDeath+dataInRange.nDeaths_strict_v1_byDateOfDeath+dataInRange.nDeaths_strict_v2_byDateOfDeath+dataInRange.nDeaths_strict_v3_byDateOfDeath);

CHRData = aggHospData./aggCasesData;
CFRData = aggDeathsData./aggCasesData;

age = 0:10:70;
tickLabels = {'0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70+'};

h = figure;
h.Position = [  560   342   848   606];
subplot(2, 2, 1)
errorbar(age, CHR_med, CHR_med-CHR_low, CHR_hi-CHR_med, 'bo-')
hold on
plot(age, CHRData, 'ro-')
xlim([-5 75])
h = gca;
h.XTick = age;
h.XTickLabel = tickLabels;
grid on
legend('model', 'data', 'Location', 'NorthWest')
ylabel('CHR')
title('(a)')
subplot(2, 2, 2)
errorbar(age, CFR_med, CFR_med-CFR_low, CFR_hi-CFR_med, 'bo-')
hold on
plot(age, CFRData, 'ro-')
xlim([-5 75])
h = gca;
h.XTick = age;
h.XTickLabel = tickLabels;
grid on
ylabel('CFR')
title('(b)')

subplot(2, 2, 3)
errorbar(age, CHR_med, CHR_med-CHR_low, CHR_hi-CHR_med, 'bo-')
hold on
plot(age, CHRData, 'ro-')
xlim([-5 75])
h = gca; h.YScale = 'log';
h.XTick = age;
h.XTickLabel = tickLabels;
grid on
xlabel('age (years)')
ylabel('CHR')
title('(c)')
subplot(2, 2, 4)
errorbar(age, CFR_med, CFR_med-CFR_low, CFR_hi-CFR_med, 'bo-')
hold on
plot(age, CFRData, 'ro-')
xlim([-5 75])
h = gca; h.YScale = 'log';
h.XTick = age;
h.XTickLabel = tickLabels;
grid on
ylabel('CFR')
title('(d)')
xlabel('age (years)')

saveas(h, 'results/CHR_CFR.png');

