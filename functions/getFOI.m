function [FOI] = getFOI(t, N, I, A, C, VEt, par)
% Function to get force of infection at time t

% Get instantaneous control function
Ct = interp1(par.date0 + (0:par.tEnd), par.Ct, t);

% Seed cases
% par.seedDur/4 - two sd either side of the mean so that seeding happens
% over approx. 7 days around dateSeed.
communitySeeds = par.initialExp * normpdf(t-par.dateSeed, 0, par.seedDur/4);
borderSeeds = par.borderSeeds * normcdf(t-par.borderTime) * par.popCount/sum(par.popCount);
Iseed = communitySeeds + borderSeeds;

% Get next generation matrix
[~, NGMclin] = getNGMtimeDep(t, par);

% Column vector of force of infection on each age group:
FOI = par.R0 * Ct / par.tI * ...
    (NGMclin * (Iseed * par.tI + sum((1-VEt) .* (I + par.cSub * A), 2))) ./ N;

end
