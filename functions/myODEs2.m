function dydt = myODEs2(t, y, par)



% Matrices for each infection state whose rows represent age group and
% columns represent immunity compartment (i.e. vax status, prior infection,
% degree of waning)
[N, V, S, E, I, A, R, C0, C1, C2, C3, Cr, H0, H1, H2, H3, Hr, F0, F1, F2, F3, Fr] = extractEpiVarsFull(y, par);

% Interpolate between original and VOC VE parameters
Alpha = normcdf( t-par.vocWaneDate(1) , 0, par.vocWaneWindow);
VEi = (1-Alpha)*par.VEi + Alpha*par.VEi_VOC;
VEs = (1-Alpha)*par.VEs + Alpha*par.VEs_VOC;
VEt = (1-Alpha)*par.VEt + Alpha*par.VEt_VOC;
VEh = (1-Alpha)*par.VEh + Alpha*par.VEh_VOC;
VEf = (1-Alpha)*par.VEf + Alpha*par.VEf_VOC;



% Get force of infection
FOI = getFOI(t, N, I, A, C0+C1+C2+C3+Cr, VEt, par);

               
% Waning flows in and out of susceptible compartments  
waneMat_StoS_dynamic = par.waneNet_StoS;
waneMat_StoS_dynamic(1:10, :) = par.waneRate_StoS * waneMat_StoS_dynamic(1:10, :);
% OJM: implement variants here via sum of normpdf at different
% times...accelerates waning. vocWaneDate -> vector

if ~par.VOC2active
    % One variant
    waneMat_StoS_dynamic(11:14, :) = (par.waneRate_StoS + par.vocWaneAmount * sum(normpdf( t-par.vocWaneDate , 0, par.vocWaneWindow))) * waneMat_StoS_dynamic(11:14, :);
    waneRate_RtoS_dynamic = (par.waneRate_RtoS + par.vocWaneAmount*par.waneRate_RtoS/par.waneRate_StoS * sum(normpdf( t-par.vocWaneDate , 0, par.vocWaneWindow)));
else
    % Two variants:
    waneMat_StoS_dynamic(11:14, :) = (par.waneRate_StoS + par.vocWaneAmount * sum(normpdf( t-par.vocWaneDate , 0, par.vocWaneWindow)) + par.vocWaneAmount2 * sum(normpdf( t-par.vocWaneDate2 , 0, par.vocWaneWindow))) * waneMat_StoS_dynamic(11:14, :);
    waneRate_RtoS_dynamic = (par.waneRate_RtoS + par.vocWaneAmount*par.waneRate_RtoS/par.waneRate_StoS * sum(normpdf( t-par.vocWaneDate , 0, par.vocWaneWindow)) + par.vocWaneAmount2*par.waneRate_RtoS/par.waneRate_StoS * sum(normpdf( t-par.vocWaneDate2, 0, par.vocWaneWindow)));
end

dSdt_wane = S * waneMat_StoS_dynamic + waneRate_RtoS_dynamic * R * par.waneNet_RtoS;

[nDose1, nDose2, nDose3, nDose4] = getDosesPerUnitTime(t, par);

% Proportion vaccinated out of the uninfected compartment(s) is the number of
% nth doses administered divided by the total number of people who have not
% yet received their n dose (which also inclues people who have been
% previously infected). For 4th (or subsequent) doses it's the number of doses divided by
% number of ppl who've had at least 3 doses
pDose1_uninf = min(1, max(0, nDose1 ./ (N-V(:, 1))   ));
pDose2_uninf = min(1, max(0, nDose2 ./ (V(:, 1)-V(:, 2))   ));
pDose3_uninf = min(1, max(0, nDose3 ./ (V(:, 2)-V(:, 3))   ));
pDose4_uninf = min(1, max(0, nDose4 ./ V(:, 3)   ));

% Proportion vaccinated out of the previously infected compartment is the
% total number of doses minus those doses applied to the
% not-previously-infected compartmemnt(s) (see above), divided by the total number of people in all 4 previously infected compartments S(:, end-3:end) 
% (this assumes noone is vaccinated while in the E, I or R compartments)
pDose1_inf = min(1, max(0, nDose1 .* (N-V(:, 1)-S(:, 1))      ./((N-V(:, 1) ).*sum(S(:, end-3:end), 2))  ));
pDose2_inf = min(1, max(0, nDose2 .* (V(:, 1)-V(:, 2)-S(:, 2))           ./((V(:, 1)-V(:, 2) ).*sum(S(:, end-3:end), 2)  )  ));
pDose3_inf = min(1, max(0, nDose3 .* (V(:, 2)-V(:, 3)-sum(S(:, 3:6), 2)) ./((V(:, 2)-V(:, 3) ).*sum(S(:, end-3:end), 2)  )  ));
pDose4_inf = min(1, max(0, nDose4 .* (V(:, 3)-sum(S(:, 7:10), 2)) ./( V(:, 3) .*sum(S(:, end-3:end), 2)  )  ));

% Overall proportion leaving each susceptible compartment as a result of vaccination 
pDoseOut_byComp = [pDose1_uninf, pDose2_uninf, pDose3_uninf.*[1 1 1 1], pDose4_uninf.*[1 1 1 1], (pDose1_inf+pDose2_inf+pDose3_inf+pDose4_inf).*[1 1 1 1] ];


dSdt_vax = (S.*pDoseOut_byComp) * par.vaxNet;


% ODEs for transmission & population dynamics
dNdt = calcPopnDynamics(N, par.popnBirthRate, par.popnDeathRate, par.popnAgeingRate);
dSdt = -(FOI*(1-VEi)) .* S + dSdt_wane + dSdt_vax               + calcPopnDynamics(S, [par.popnBirthRate, zeros(1, par.nSusComp-1)], par.popnDeathRate, par.popnAgeingRate);
dEdt =  (FOI*(1-VEi)) .* S                      - 1/par.tE * E  + calcPopnDynamics(E, zeros(1, par.nSusComp), par.popnDeathRate, par.popnAgeingRate);
dIdt = 1/par.tE * (    par.pClin.*(1-VEs)./(1-VEi)) .* E - 1/par.tI * I  + calcPopnDynamics(I, zeros(1, par.nSusComp), par.popnDeathRate, par.popnAgeingRate);
dAdt = 1/par.tE * (1 - par.pClin.*(1-VEs)./(1-VEi)) .* E - 1/par.tI * A  + calcPopnDynamics(A, zeros(1, par.nSusComp), par.popnDeathRate, par.popnAgeingRate);
dRdt = 1/par.tI * (I+A) - waneRate_RtoS_dynamic * R                     + calcPopnDynamics(R, zeros(1, par.nSusComp), par.popnDeathRate, par.popnAgeingRate);
dVdt = [nDose1, nDose2, nDose3]                                     + calcPopnDynamics(V, zeros(1, 3), par.popnDeathRate, par.popnAgeingRate);


% Testing and disease model
% newXXXToBei is the number new people per day exiting the latent period
% who will eventually experience outcome XXX (case/hosp/death) and whose
% immunity status is i = 0, 1, 2 or 3+ doses or reinfection. Each variable
% is a column vector of 16 age groups

[pTestClinNow, pTestSubNow] = getPTestPerUnitTime(t, par.date0, ...
    par.tEnd, par.pTestClin, par.pTestSub);
newCasesToBe_all = 1/par.tE * (pTestClinNow .* par.pClin .* ...
    (1-VEs)./(1-VEi) + pTestSubNow .* (1 - par.pClin.*(1-VEs)./(1-VEi))) .* E;

newCasesToBe0 = newCasesToBe_all(:, 1);
newCasesToBe1 = newCasesToBe_all(:, 2);
newCasesToBe2 = sum(newCasesToBe_all(:, 3:6), 2);
newCasesToBe3 = sum(newCasesToBe_all(:, 7:10), 2);
newCasesToBer = sum(newCasesToBe_all(:, 11:14), 2);

[antiviralsTreatedPropNow] = getAntiviralsTreatedPerUnitTime(t, ...
    par.date0, par.tEnd, par.daily_treatcaseratio);

newHospToBe_all = 1/par.tE * (1-VEh)./(1-VEi) .* par.IHR .* E .* ...
    (1 - antiviralsTreatedPropNow .* par.antiviralsEffectIHRmult);
newHospToBe0 = newHospToBe_all(:, 1);
newHospToBe1 = newHospToBe_all(:, 2);
newHospToBe2 = sum(newHospToBe_all(:, 3:6), 2);
newHospToBe3 = sum(newHospToBe_all(:, 7:10), 2);
newHospToBer = sum(newHospToBe_all(:, 11:14), 2);

newDeathsToBe_all = 1/par.tE * (1-VEf)./(1-VEi) .* par.IFR .* E .* ...
    (1 - antiviralsTreatedPropNow .* par.antiviralsEffectIFRmult);
newDeathsToBe0 = newDeathsToBe_all(:, 1);
newDeathsToBe1 = newDeathsToBe_all(:, 2);
newDeathsToBe2 = sum(newDeathsToBe_all(:, 3:6), 2);
newDeathsToBe3 = sum(newDeathsToBe_all(:, 7:10), 2);
newDeathsToBer = sum(newDeathsToBe_all(:, 11:14), 2);

% Progression through unobserved->observed compartments
dC0dt = get_dCdt(newCasesToBe0, C0, par);
dC1dt = get_dCdt(newCasesToBe1, C1, par);
dC2dt = get_dCdt(newCasesToBe2, C2, par);
dC3dt = get_dCdt(newCasesToBe3, C3, par);
dCrdt = get_dCdt(newCasesToBer, Cr, par);

dH0dt = get_dHdt(newHospToBe0, H0, par);
dH1dt = get_dHdt(newHospToBe1, H1, par);
dH2dt = get_dHdt(newHospToBe2, H2, par);
dH3dt = get_dHdt(newHospToBe3, H3, par);
dHrdt = get_dHdt(newHospToBer, Hr, par);

dF0dt = get_dFdt(newDeathsToBe0, F0, par);
dF1dt = get_dFdt(newDeathsToBe1, F1, par);
dF2dt = get_dFdt(newDeathsToBe2, F2, par);
dF3dt = get_dFdt(newDeathsToBe3, F3, par);
dFrdt = get_dFdt(newDeathsToBer, Fr, par);


% Restructure as column vector
dydt = reshape([dNdt, dVdt, dSdt, dEdt, dIdt, dAdt, dRdt(:, 1:end-1), dC0dt, dC1dt, dC2dt, dC3dt, dCrdt, dH0dt, dH1dt, dH2dt, dH3dt, dHrdt, dF0dt, dF1dt, dF2dt, dF3dt, dFrdt], length(y), 1);

end

