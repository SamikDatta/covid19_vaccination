function [nDeathsActualM, nDeathsActualF] = getActualDeaths(fName, age)


date0 = datetime(2022, 1, 1);       % start date for analysis  25JAN2022
date1 = datetime(2023, 6, 30);

% Load case linelist data 
load(fName, 'cases');

ageExt = [age; 130];
inFlag = cases.DIED_DT >= date0 & cases.DIED_DT <= date1 & cases.Died == "Yes" & ismember(cases.COD_SUMMARY, categorical(["COVID as contributory", "COVID as underlying", "Not available"]));
nDeathsActualM = histcounts(cases.Age(inFlag & cases.SEX == "Male"), ageExt)';
nDeathsActualF = histcounts(cases.Age(inFlag & cases.SEX == "Female"), ageExt)';
nDeathsActualU = histcounts(cases.Age(inFlag & cases.SEX == "Unknown"), ageExt)';
nDeathsActualX = sum(inFlag & ~(cases.Age >= ageExt(1) & cases.Age <= ageExt(end)));      % deaths with no valid age (0-130) recorded 
nDeathsActualU = nDeathsActualU + mnrnd(nDeathsActualX, (nDeathsActualM+nDeathsActualF+nDeathsActualU)/sum(nDeathsActualM+nDeathsActualF+nDeathsActualU))';      % distribute ageless deaths according to distribution of deaths of known age


% Distribute deaths of unknown sex in proportion to known male and female
% deaths (splitting 50-50 in age bands where there were 0 recorded male and female deaths):
pM = nDeathsActualM./(nDeathsActualM+nDeathsActualF);
pM(isnan(pM)) = 0.5;
extraDeathsM = binornd(nDeathsActualU, pM);
extraDeathsF = nDeathsActualU - extraDeathsM;

nDeathsActualM = nDeathsActualM + extraDeathsM;
nDeathsActualF = nDeathsActualF + extraDeathsF;



