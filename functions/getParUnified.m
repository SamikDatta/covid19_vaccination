function par = getParUnified(Theta, parBase, scenPar)
% Function to create structure of parameters that can change in the fitting
% procedures or in scenarios
% INPUTS:
% - Theta: table of parameter values to be fitted, or posterior
% - parBase: structure of fixed parameters, defined in getParBase.m
% - scenPar: if running scenarios, table of scenario-specific parameters
% OUTPUT:
% - par: structure of parameters


% Functions to add a zero-centred random variable to selected parameters
% according to random deviates [0,1] specified in input parameter Theta
plusMinus = @(z, r)(2*r*z - r);                     % uniform perturbation between +/- r, z is a random deviate [0,1]
plusMinusInt = @(z, r)(floor((2*r+1)*z) - r );      % uniform perturbation on integers between +/- r, z is a random deviate [0,1]


%--------------------- Seed date --------------------
par.dateSeed = datenum('19JAN2022') + plusMinusInt(Theta.dateSeed, 3);


% -------------------- Control function parameters -----------------------

% Starting value of Ct (fitted) - between 0.58-0.78
par.Ct = (0.68 + plusMinus(Theta.Cstart, 0.1)) * ones(1, parBase.tEnd+1);

CtRampStarts = [datenum('10-Mar-2022') + plusMinusInt(Theta.rampStart, 5), ...
    datenum('15-Sep-2022') + plusMinusInt(Theta.ramp2Start, 5), ...
    parBase.date0 + parBase.tEnd];
CtRampDays = [55 + plusMinusInt(Theta.rampDays, 20), ...
    10 + plusMinusInt(Theta.ramp2Days, 9), 0];

% Ct after 1st ramp up (0.89-1.31) 
% and after 2nd ramp up (0.89-1.31 * 1.1-1.3)
CtRamp = [1.1 + plusMinus(Theta.Cramp, 0.21), ...
    (1.1 + plusMinus(Theta.Cramp, 0.21)) * (1.2 + plusMinus(Theta.Cramp2, 0.1)), ...
    (1.1 + plusMinus(Theta.Cramp, 0.21)) * (1.2 + plusMinus(Theta.Cramp2, 0.1))];

if exist('scenPar', 'var')
    % If running a transmission change scenario, dummy third element of 
    % Ct-specific parameters is set to the appropriate value
    CtRampStarts(3) = datenum(scenPar.policyDate);
    CtRampDays(3) = 30;
    CtRamp(3) = CtRamp(2) * scenPar.cRampDeltaPolicy;
end

% Adding Ct ramp-ups at each date:
for pci = 1:length(CtRampStarts)
    ti = datenum(CtRampStarts(pci)) - parBase.date0;
    par.Ct(ti:ti+CtRampDays(pci)-1) = linspace(par.Ct(ti), CtRamp(pci), CtRampDays(pci));
    par.Ct(ti+CtRampDays(pci):end) = CtRamp(pci);
    par.Ct = par.Ct(1:parBase.tEnd+1); % Making sure size stays the same
end

%%%%% WINTER MODELLING (added on 30Mar)
% Default parameter (weak seasonality effect), change to 0.2 for strong
% effect, or 0 for no seasonality
seasonMultAmpl = 0; 
if exist('scenPar', 'var')
        seasonMultAmpl = scenPar.seasonMultAmp;
end
seasonStart = datenum("01-Apr-2023"); % seasonality model starts on 1Apr23
CtSeasonMult = 1 + seasonMultAmpl * (parBase.tBase >= seasonStart) .* ...
    sin(2 * pi * (parBase.tBase - seasonStart) / 365);

% Added season multiplier
par.Ct = par.Ct .* CtSeasonMult;

CtFigure = 0;
if CtFigure == 1
    figure
    plot(datetime(parBase.tBase(302:end), 'ConvertFrom', 'datenum'), par.Ct(302:end), 'LineWidth', 2)
    hold on
    grid on
    grid minor
    xlim(datetime([parBase.tBase(302), parBase.tBase(end)], 'ConvertFrom', 'datenum'))
    ylabel('Control function C(t)')
end


% --------------------- Testing and lag parameters -----------------------
% Overall scaling constant for all testing parameters
par.pTestClin = parBase.pTestClin0 .* (1 + plusMinus(Theta.pTestMult, 0.2));

% Subclinical cases have a lower testing probability
par.pTestSub = parBase.subClinPtestMult .* par.pTestClin;                           


% ------------- Disease Rate Multipliers --------------
% Overall scaling constants for IHR and IFR
par.IHRmult = 1 + plusMinus(Theta.IHR, 0.5);  
par.IFRmult = 1 + plusMinus(Theta.IFR, 0.5);      

par.IHR = par.IHRmult * parBase.IHR0;
par.IFR = par.IFRmult * parBase.IFR0;


% ---------------- Effect of antivirals -------------------
% Multiplier for the effect of antivirals on par.IHR and par.IFR.
% 0 - no effect; 1 - full protection from outcome
par.antiviralsEffectIHRmult = 0; 
par.antiviralsEffectIFRmult = 0.5 + plusMinus(Theta.aViralEffect, 0.1); 
if exist('scenPar', 'var')
        par.antiviralsEffectIHRmult = scenPar.antiviralsEffectIHRmult;
end


par.pTestTS = (par.pTestClin .* parBase.pClin + par.pTestSub .* (1 - parBase.pClin));


% --------------- Contact matrix adjustment ------------------
% Amount by which contract matrix relaxes back to Prem (0=not at all, 1=fully)
par.relaxAlpha = 0.4 + plusMinus(Theta.relaxAlpha, 0.4);             

Cw2 = (1-par.relaxAlpha) * parBase.Cw1 + ...
    par.relaxAlpha * triu(ones(parBase.nAgeBlocks));
Cw1 = parBase.Cw1 + triu(parBase.Cw1, 1)';                % make weights into symmetric matrices
Cw2 = Cw2 + triu(Cw2, 1)';
par.contactPar.weights = repelem(Cw1, parBase.ageBlockSizes, parBase.ageBlockSizes);        % expand blocks to create a 16x16 matrix that can multiplied elementwise with the contact matrix
par.contactPar.weightsChange = repelem(Cw2, parBase.ageBlockSizes, parBase.ageBlockSizes);

% Start date and time window for change of contact matrix - matrix will
% change linearly during the specified number of days followiong the start
% date
par.contactPar.changeDate = CtRampStarts(1);
par.contactPar.changeWindow = 70 + plusMinusInt(Theta.MRampDays, 20);



%---------------------------  VOC model ---------------------------------
% Coefficient determining what fraction of each post-infection susceptible
% compartment gets bumped down the immunity scale
par.vocWaneAmount = 0.4 + plusMinus(Theta.vocWane, 0.3) ;

par.VOC2active = false;

if exist('scenPar', 'var') && scenPar.VOC2active
    par.VOC2active = scenPar.VOC2active;
    par.vocWaneDate2 = datenum(scenPar.VOC2date); 
	par.vocWaneAmount2 = 0.25; 
end


%------------------------- Immunity parameters --------------------------
par.waneRateMult = 1 + plusMinus(Theta.waneRate, 0.5);    % Fitted multiplier on waning rate
par.waneRate_StoS = parBase.waneRateMean * par.waneRateMult;  % Rate of moving from one S compartment to the next one with lower immunity
par.waneRate_RtoS = parBase.waneRateMean * par.waneRateMult * parBase.relRate_RtoS; % Rate of moving R to S


%------------------------- Vaccination coverage --------------------------

% By default, the number of doses in each age group is as defined in getBasePar.m
par.nDoses1Smoothed = parBase.nDoses1Smoothed0;
par.nDoses2Smoothed = parBase.nDoses2Smoothed0;
par.nDoses3Smoothed = parBase.nDoses3Smoothed0;
par.nDoses4Smoothed = parBase.nDoses4Smoothed0;


if exist('scenPar', 'var')      % If running scenarios

    % Vaccination scenarios to run:
	% (1) as usual
	% (2) No vaccination
	% (3) No antivirals
	% (4) No vaccination and no antivirals
	% (5) No vaccination in U60s
	% (6) Vaccination at 90% (or x%) of what actually happened
	% (7) Vaccination rates same as 20-25 yos in all age groups (the China scenario)  
	% (8) Vaccination as for the Māori population (I have the data for this)
	% (9) Vaccination as for the European / other population


    % Create a 4x16 matrix of reduction multipliers to apply to the number
    % of vaccination doses in each age group.
    vaxCovRed = zeros(4, 16);
%     vaxCovRed(parBase.dosesToRed, parBase.agesToRed) = scenPar.vaxCovRedFactor;
    vaxCovRed(scenPar.dosesToRed{:}, scenPar.agesToRed{:}) = scenPar.vaxCovRedFactor; % get working with cell arrays
    
    % Apply reduction multipliers to the number of vax doses
    par.nDoses1Smoothed = parBase.nDoses1Smoothed0 .* (1 - vaxCovRed(1, :));
    par.nDoses2Smoothed = parBase.nDoses2Smoothed0 .* (1 - vaxCovRed(2, :));
    par.nDoses3Smoothed = parBase.nDoses3Smoothed0 .* (1 - vaxCovRed(3, :));
    par.nDoses4Smoothed = parBase.nDoses4Smoothed0 .* (1 - vaxCovRed(4, :));

    if scenPar.scenarioNumber == 7 % flag for China-like vaccination - copy 20-25 to all columns (scale by population size)
        par.nDoses1Smoothed(:, 4:16) = bsxfun(@times, repmat(parBase.nDoses1Smoothed0(:, 5)/parBase.popCount(5), 1, 13), parBase.popCount(4:16)'); 
        par.nDoses2Smoothed(:, 4:16) = bsxfun(@times, repmat(parBase.nDoses2Smoothed0(:, 5)/parBase.popCount(5), 1, 13), parBase.popCount(4:16)'); 
        par.nDoses3Smoothed(:, 4:16) = bsxfun(@times, repmat(parBase.nDoses3Smoothed0(:, 5)/parBase.popCount(5), 1, 13), parBase.popCount(4:16)'); 
        par.nDoses4Smoothed(:, 4:16) = bsxfun(@times, repmat(parBase.nDoses4Smoothed0(:, 5)/parBase.popCount(5), 1, 13), parBase.popCount(4:16)'); 
    elseif scenPar.scenarioNumber == 8 % flag for Maori vaccination rates - scale by population size
        [~, tempdoses1, tempdoses2, tempdoses3, tempdoses4plus] = ...
            getVaccineData('data/', "vaccine_data_Maori_2023-06-06", "reshaped_b2_projections_final_2022-07-13.csv", ...
            parBase.vaccImmDelay, parBase.date0, parBase.tEnd); % get vax doses from Maori uptake
        if parBase.sensitivity_flag == 1
            hsu_pop = readtable('data/popproj2018-21'); % read in population numbers (alternative, from Stats NZ)
        else
            hsu_pop = readtable('data/HSU_by_age_eth'); % read in population numbers
        end
        pop_scale_factor = [hsu_pop.Total(1:15); sum(hsu_pop.Total(16:19))]./[hsu_pop.Maori(1:15); sum(hsu_pop.Maori(16:19))]; % scaling Maori population up to total
        scaled_doses1 = bsxfun(@times, tempdoses1, pop_scale_factor'); % rescale doses up to population level
        scaled_doses2 = bsxfun(@times, tempdoses2, pop_scale_factor'); 
        scaled_doses3 = bsxfun(@times, tempdoses3, pop_scale_factor'); 
        scaled_doses4plus = bsxfun(@times, tempdoses4plus, pop_scale_factor'); 
        smoothWindow = 56;
        par.nDoses1Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses1, 'movmean', smoothWindow));]; % smooth out in same way
        par.nDoses2Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses2, 'movmean', smoothWindow))];
        par.nDoses3Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses3, 'movmean', smoothWindow))];
        par.nDoses4Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses4plus, 'movmean', smoothWindow))];
    elseif scenPar.scenarioNumber == 9 % flag for European / Other vaccination rates - scale by population size
        [~, tempdoses1, tempdoses2, tempdoses3, tempdoses4plus] = ...
            getVaccineData('data/', "vaccine_data_EuropeanOther_2023-06-06", "reshaped_b2_projections_final_2022-07-13.csv", ...
            parBase.vaccImmDelay, parBase.date0, parBase.tEnd); % get vax doses from Maori uptake
        if parBase.sensitivity_flag == 1
            hsu_pop = readtable('data/popproj2018-21'); % read in population numbers (alternative, from Stats NZ)
        else
            hsu_pop = readtable('data/HSU_by_age_eth'); % read in population numbers
        end
        pop_scale_factor = [hsu_pop.Total(1:15); sum(hsu_pop.Total(16:19))]./[hsu_pop.EuropeanorOther(1:15); sum(hsu_pop.EuropeanorOther(16:19))]; % scaling Maori population up to total
        scaled_doses1 = bsxfun(@times, tempdoses1, pop_scale_factor'); % rescale doses up to population level
        scaled_doses2 = bsxfun(@times, tempdoses2, pop_scale_factor'); 
        scaled_doses3 = bsxfun(@times, tempdoses3, pop_scale_factor'); 
        scaled_doses4plus = bsxfun(@times, tempdoses4plus, pop_scale_factor'); 
        smoothWindow = 56;
        par.nDoses1Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses1, 'movmean', smoothWindow));]; % smooth out in same way
        par.nDoses2Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses2, 'movmean', smoothWindow))];
        par.nDoses3Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses3, 'movmean', smoothWindow))];
        par.nDoses4Smoothed = [zeros(1, parBase.nAgeGroups); diff(smoothdata(scaled_doses4plus, 'movmean', smoothWindow))];        
    end

end


end

