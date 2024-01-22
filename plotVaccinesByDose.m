clear
close all

addpath('functions');

dateLbl = "NA";     % no date label needed as epi data is not going to be used in this script

% Get data file names
[myDataPath, dataFileNames] = getDataFileNames(dateLbl);

% Read in national vaccination data 
tbl = readtable(myDataPath + dataFileNames.vaxDataFname);

% Read in HSU population size data by age (19 x 5 year age bands)
popSize = readmatrix(myDataPath + dataFileNames.popSizeFname);

% Calculate totla pop size and in various age categories
popSizeTot = sum(popSize(:, 2));
popSizeOver5 = sum(popSize(2:end, 2));
popSizeOver16 = 4/5*popSize(4, 2) + sum(popSize(5:end, 2));
popSizeOver50 = sum(popSize(11:end, 2));


% Calculate total 1st, 2nd, 3rd and 4th or subsequent doses across all ages:
tbl.doses1 = sum(table2array(tbl(:, 2:17)), 2);
tbl.doses2 = sum(table2array(tbl(:, 18:33)), 2);
tbl.doses3 = sum(table2array(tbl(:, 34:49)), 2);
tbl.doses4s = sum(table2array(tbl(:, 50:65)), 2);




% Plot figure for Supp Material
h = figure(1);
h.Position = [ 677   331   920   499];
plot(tbl.dates, tbl.doses1/popSizeTot, tbl.dates, tbl.doses2/popSizeTot, tbl.dates, tbl.doses3/popSizeTot, tbl.dates, tbl.doses4s/popSizeTot, 'LineWidth', 2)
ylabel('doses per capita')
legend('1st dose', '2nd dose', '3rd dose', '4th or subsequent dose', 'Location', 'northwest')
yline(popSizeOver5/popSizeTot, 'k--', 'HandleVisibility', 'off')
yline(popSizeOver16/popSizeTot, 'k--', 'HandleVisibility', 'off')
yline(popSizeOver50/popSizeTot, 'k--', 'HandleVisibility', 'off')
grid on


