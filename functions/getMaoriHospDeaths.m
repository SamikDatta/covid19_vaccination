function [nMaoriHosp, nMaoriDeaths] = getMaoriHospDeaths(fName, age)


date0 = datetime(2022, 1, 1);       % start date for analysis  25JAN2022
date1 = datetime(2023, 6, 30);

% Load case linelist data 
load(fName, 'cases');

ageExt = [age; 130];
inFlag = cases.ADMISSION_DT >= date0 & cases.ADMISSION_DT <= date1 & cases.COVID_RELATED_HOSPITALISATION == "1" & cases.ETHNICITYPRIORITISEDMPAMEU == "Maori";
nMaoriHosp = histcounts(cases.Age(inFlag), ageExt);
inFlag = cases.DIED_DT >= date0 & cases.DIED_DT <= date1 & cases.Died == "Yes" & ismember(cases.COD_SUMMARY, categorical(["COVID as contributory", "COVID as underlying", "Not available"])) & cases.ETHNICITYPRIORITISEDMPAMEU == "Maori";
nMaoriDeaths = histcounts(cases.Age(inFlag), ageExt);


