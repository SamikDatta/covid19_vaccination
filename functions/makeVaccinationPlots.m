function makeVaccinationPlots(parBase, parsToFit)

t = parBase.tBase;

all_dates = datetime(parBase.tBase,'ConvertFrom','datenum');

vac_table = table('Size', [parBase.nScenarios 6], ...
    'VariableTypes', ["string", repmat("cell", 1, 5)], ...
    'VariableNames', ["scenario", "first dose", "second dose", ...
    "third dose", "fourth dose", "total doses"]);
vac_table.scenario = ["baseline", "No vaccination", "No antivirals", "No vaccination or antivirals", ...
    "No vaccination in U60s", "90% Vaccination", "reduced coverage in older ages", ...
    "Maori vaccination", "European / Other vaccination"]';

nParsToFit = length(parsToFit);
ThetaTemp = array2table(zeros(1, nParsToFit), 'VariableNames', parsToFit);      % Just use a dummy value of Theta=0 for all fitted parameters as these do not affect the vaccine plots

for i = 1:parBase.nScenarios

    parScen = parBase.scenarios(i, :);
    parInd = getParUnified(ThetaTemp, parBase, parScen);
    vac_table.("first dose"){i} = parInd.nDoses1Smoothed;
    vac_table.("second dose"){i} = parInd.nDoses2Smoothed;
    vac_table.("third dose"){i} = parInd.nDoses3Smoothed;
    vac_table.("fourth dose"){i} = parInd.nDoses4Smoothed;
    vac_table.("total doses"){i} = vac_table.("first dose"){i} + vac_table.("second dose"){i} + ...
        vac_table.("third dose"){i} + vac_table.("fourth dose"){i};

end

% Plot cumulative doses for each scenario

legLbls = ["Baseline", "10% drop", "20-25 yo rates", "Euro rates", "Māori rates"];
f = figure(1);  
f.Position = [100 100 1200 800];
titles = {'0-5', '5-10', '10-15', '15-20', '20-25', '25-30', ...
    '30-35', '35-40', '40-45', '45-50', '50-55', '55-60', '60-65', ...
    '65-70', '70-75', '75+'};
tiledlayout(4, 4);
for ag = 1:16
    nexttile
    hold on
    title(titles(ag))
    plot(all_dates, cumsum(vac_table.("total doses"){1}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){6}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){7}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){9}(:, ag))/parBase.popCount(ag), ...
        all_dates, cumsum(vac_table.("total doses"){8}(:, ag))/parBase.popCount(ag) );
    hold off
    xlim([datetime(2021, 1, 1), datetime(2023, 7, 1)] )
    ylim([0 4])
    grid on
    if ag == 1
        legend(legLbls, 'Location', 'NorthWest');
    end
    if mod(ag, 4) == 1
        ylabel('doses per capita')
    end
end

drawnow


% Print some summary ethnicity stats for paper
ti = find(t == datenum('01FEB2022'));
AM = sum(vac_table.("second dose"){8}(1:ti,:))
AE = sum(vac_table.("second dose"){9}(1:ti,:))
pM1 = sum(AM(4:13))/sum(parBase.popCount(4:13))
pM2 = sum(AM(14:end))/sum(parBase.popCount(14:end))
pE1 = sum(AE(4:13))/sum(parBase.popCount(4:13))
pE2 = sum(AE(14:end))/sum(parBase.popCount(14:end))

BM = sum(vac_table.("third dose"){8}(1:ti,:))
BE = sum(vac_table.("third dose"){9}(1:ti,:))
qM1 = sum(BM(4:13))/sum(parBase.popCount(4:13))
qM2 = sum(BM(14:end))/sum(parBase.popCount(14:end))
qE1 = sum(BE(4:13))/sum(parBase.popCount(4:13))
qE2 = sum(BE(14:end))/sum(parBase.popCount(14:end))
