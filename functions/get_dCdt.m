function dCdt = get_dCdt(newCasesToBe, C, par)

dCdt = zeros(par.nAgeGroups, par.nCaseComp);

% Progression through unobserved->observed case compartments
dCdt(:, 1) = newCasesToBe - 2/par.tLatentToTest * C(:, 1);
dCdt(:, 2) = 2/par.tLatentToTest * C(:, 1) -  2/par.tLatentToTest * C(:, 2);
dCdt(:, 3) = 2/par.tLatentToTest * C(:, 2);



