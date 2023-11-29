function [age, exM, exF] = importLifeTables(fName)

% Import life tables and extract life expectancy at age  (expected years of life remaining) in 2022
% age contains a vector of ages
% exM contains a matrix of life expectancies for each age for males in 2022
% exF contains a matrix of life expectancies for each age for males in 2022
% In each of exM and exF, the three columns are the 5th percentile, the
% median and the 95th percentile



year0 = 2022;

opts = detectImportOptions(fName);
opts = setvartype(opts, {'percentile', 'sex'}, 'categorical');
tbl = readtable(fName, opts);

ind = tbl.yearofbirth+tbl.age == year0 | (tbl.yearofbirth == max(tbl.yearofbirth) & tbl.age == 0);
ageMMed = tbl.age(ind & tbl.sex == "male" & tbl.percentile == "median");
exMMed = tbl.ex(ind & tbl.sex == "male"  & tbl.percentile == "median");
ageMLow = tbl.age(ind & tbl.sex == "male"  & tbl.percentile == "5");
exMLow = tbl.ex(ind & tbl.sex == "male" & tbl.percentile == "5");
ageMHi = tbl.age(ind & tbl.sex == "male" & tbl.percentile == "95");
exMHi = tbl.ex(ind & tbl.sex == "male" & tbl.percentile == "95");
exM = [exMLow, exMMed, exMHi];

ageFMed = tbl.age(ind & tbl.sex == "female" & tbl.percentile == "median");
exFMed = tbl.ex(ind & tbl.sex == "female"  & tbl.percentile == "median");
ageFLow = tbl.age(ind & tbl.sex == "female"  & tbl.percentile == "5");
exFLow = tbl.ex(ind & tbl.sex == "female" & tbl.percentile == "5");
ageFHi = tbl.age(ind & tbl.sex == "female" & tbl.percentile == "95");
exFHi = tbl.ex(ind & tbl.sex == "female" & tbl.percentile == "95");
exF = [exFLow, exFMed, exFHi];

assert(isequal(ageMLow, ageMHi, ageMMed, ageFLow, ageFHi, ageFMed));
age = ageMMed;

[age, si] = sort(age, 'ascend');
exM = exM(si, :);
exF = exF(si, :);

