function dHdt = get_dHdt(newHospToBe, H, par)

dHdt = zeros(par.nAgeGroups, par.nHospComp);

% Progression through unobserved->observed hospital compartments
dHdt(:, 1) = newHospToBe - 2/par.tLatentToTest * H(:, 1);
dHdt(:, 2) = 2/par.tLatentToTest * H(:, 1) -  2/par.tLatentToTest * H(:, 2);
dHdt(:, 3) = 2/par.tLatentToTest * H(:, 2) - 1/par.tTestToHosp * H(:,3);
dHdt(:, 4) = 1/par.tTestToHosp * H(:,3) - 1./par.tLOS .* H(:, 4);
dHdt(:, 5) = 1./par.tLOS .* H(:, 4);

