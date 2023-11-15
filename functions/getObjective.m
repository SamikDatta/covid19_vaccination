function obj = getObjective(t, epiVarsCompactTemp, parTemp)
% Function that computes the distance (objective) to the data for a given 
% parameter set
% INPUT:
% - t: time vector
% - epiVarsCompactTemp: structure containing ODE model solution
% - parTemp: parameter structure obtained to solve ODE model
% OUTPUT:
% - obj: sum of errors for given parameter set

% vars are 
% daily cases
% daily deaths
% daily new hospital admissions
% daily incidence per capita
% age breakdown of cases
% age breakdown of hospitalisations
dataVars = {'nCasesData','nDeathsData','nHospData', 'NationalBorder', 'CasesByAge', 'HospByAge'};

% calculate error function on data. extra outputs (nObs, resid, dist) for diagnostics
[~, ~, dist] = calcErrorWithAgeCurrent(t, epiVarsCompactTemp, dataComb, parTemp, dataVars);
errWeights = ones(1, length(dataVars)); 
obj = sum(errWeights .* dist, 2);


end