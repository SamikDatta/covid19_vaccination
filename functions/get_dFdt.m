function dFdt = get_dFdt(newDeathsToBe, F, par)

dFdt = zeros(par.nAgeGroups, par.nDeathComp);

% Progression through unobserved->observed death compartments
dFdt(:, 1) = newDeathsToBe - 2/par.tLatentToTest * F(:, 1);
dFdt(:, 2) = 2/par.tLatentToTest * F(:, 1) -  2/par.tLatentToTest * F(:, 2);
dFdt(:, 3) = 2/par.tLatentToTest * F(:, 2) - 1/par.tTestToHosp * F(:,3);
dFdt(:, 4) = 1/par.tTestToHosp * F(:,3) - 2/par.tDeath .* F(:, 4);
dFdt(:, 5) = 2/par.tDeath .* F(:, 4) - 2/par.tDeath * F(:, 5);
dFdt(:, 6) = 2/par.tDeath * F(:, 5);


