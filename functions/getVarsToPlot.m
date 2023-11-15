function [newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, ...
    newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, ...
    newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, ...
    newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, ...
    S, Sw, E1, E2, I, A] = getVarsToPlot(epiVarsCompact)

% Extract relevant time series by concatenating across (accepted) realisations
% Returns a 3D array (stack of matrices) where each matrix is a n x m
% matrix where n = number of time points and m = 16 is the number of age
% groups


% Daily cases
newDailyCases0 = cat(3, epiVarsCompact.C0);
newDailyCases0(2:end, :, :) = diff(newDailyCases0);        % diff of cumulative cases (cumsum of newDailyCases should be same as epiVarsCompact.C)
newDailyCases1 = cat(3, epiVarsCompact.C1);
newDailyCases1(2:end, :, :) = diff(newDailyCases1);        % diff of cumulative cases (cumsum of newDailyCases should be same as epiVarsCompact.C)
newDailyCases2 = cat(3, epiVarsCompact.C2);
newDailyCases2(2:end, :, :) = diff(newDailyCases2);        % diff of cumulative cases (cumsum of newDailyCases should be same as epiVarsCompact.C)
newDailyCases3 = cat(3, epiVarsCompact.C3);
newDailyCases3(2:end, :, :) = diff(newDailyCases3);        % diff of cumulative cases (cumsum of newDailyCases should be same as epiVarsCompact.C)
newDailyCasesr = cat(3, epiVarsCompact.Cr);
newDailyCasesr(2:end, :, :) = diff(newDailyCasesr);        % diff of cumulative cases (cumsum of newDailyCases should be same as epiVarsCompact.C)

% Daily new hospital admissions
newDailyHosp0 = cat(3, epiVarsCompact.H0occ) + cat(3, epiVarsCompact.H0dis);
newDailyHosp0(2:end, :, :) = diff(newDailyHosp0);
newDailyHosp1 = cat(3, epiVarsCompact.H1occ) + cat(3, epiVarsCompact.H1dis);
newDailyHosp1(2:end, :, :) = diff(newDailyHosp1);
newDailyHosp2 = cat(3, epiVarsCompact.H2occ) + cat(3, epiVarsCompact.H2dis);
newDailyHosp2(2:end, :, :) = diff(newDailyHosp2);
newDailyHosp3 = cat(3, epiVarsCompact.H3occ) + cat(3, epiVarsCompact.H3dis);
newDailyHosp3(2:end, :, :) = diff(newDailyHosp3);
newDailyHospr = cat(3, epiVarsCompact.Hrocc) + cat(3, epiVarsCompact.Hrdis);
newDailyHospr(2:end, :, :) = diff(newDailyHospr);

% Daily deaths
newDailyDeaths0 = cat(3, epiVarsCompact.F0);
newDailyDeaths0(2:end, :, :) = diff(newDailyDeaths0);       % diff of cumulative deaths
newDailyDeaths1 = cat(3, epiVarsCompact.F1);
newDailyDeaths1(2:end, :, :) = diff(newDailyDeaths1);       % diff of cumulative deaths
newDailyDeaths2 = cat(3, epiVarsCompact.F2);
newDailyDeaths2(2:end, :, :) = diff(newDailyDeaths2);       % diff of cumulative deaths
newDailyDeaths3 = cat(3, epiVarsCompact.F3);
newDailyDeaths3(2:end, :, :) = diff(newDailyDeaths3);       % diff of cumulative deaths
newDailyDeathsr = cat(3, epiVarsCompact.Fr);
newDailyDeathsr(2:end, :, :) = diff(newDailyDeathsr);       % diff of cumulative deaths

% Incidence of new infections
E1 = cat(3, epiVarsCompact.E1);
% incidence1 = 1/parBase.tE * E1;
% incidence1(2:end, :, :) = incidence1(2:end, :, :) + diff(E1);      % Incidence is dE/dt + 1/tE * E
E2 = cat(3, epiVarsCompact.E2);
% incidence2 = 1/parBase.tE * E2;
% incidence2(2:end, :, :) = incidence2(2:end, :, :) + diff(E2);

% Hosipital occupancy
Hocc = cat(3, epiVarsCompact.H0occ) + cat(3, epiVarsCompact.H1occ) + cat(3, epiVarsCompact.H2occ) + cat(3, epiVarsCompact.H3occ) + cat(3, epiVarsCompact.Hrocc);

% S, Sw, I and A variables    
S = cat(3, epiVarsCompact.S);
Sw = cat(3, epiVarsCompact.Sw);
I = cat(3, epiVarsCompact.I);
A = cat(3, epiVarsCompact.A);

end
