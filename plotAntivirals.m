clear 
close all

smoothWindow = 28;

load('data/therapeutics_by_age_14-Aug-2023.mat')

t = outTab.date;
pTherap = outTab.nTreated./outTab.nCases;
pTherapSmoothed = movmean(outTab.nTreated./outTab.nCases, smoothWindow, 1);


colOrd = colororder;
colOrd = [colOrd; [0 0 0]; [1 0 1]; [0 1 1]; ];

figure(1)
set(gcf, 'DefaultAxesColorOrder', colOrd);
set(gcf, 'DefaultAxesLineStyleOrder', {'--','-'});
plot(t, pTherapSmoothed(:, 4:end))
xlim([datetime(2022, 3, 1), datetime(2023, 6, 30) ])
ylabel('proportion of cases treated')
legend(["15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90+"], 'Location', 'NorthWest')

