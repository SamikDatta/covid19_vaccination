function par = getBasePar(tMax, myDataPath, dataFileNames)
% Function to create structure of base parameters that stay fixed for all
% simulations
% INPUT:
% - tMax: datenum variable corresponding to last datapoint
% - myDataPath: path to folder containing data
% - dataFileNames: structure with fields containing the filenames of the data files to be read in
% OUTPUT:
% - par: structure of parameters


% Simulation start date (start date set to 5Mar21 to include full
% vaccine rollout period and simulate correct waning dynamics)
par.date0 = datenum('05MAR2021');       

% Simulation end date, if no tMax specified run for 3 years
if isnan(tMax); par.tEnd = 3 * 360; else; par.tEnd = tMax - par.date0; end    

% Create time vector
par.tBase = par.date0:par.date0+par.tEnd;

%------------- SEIR parameters --------------
par.R0 = 3.25;   % 6
par.cSub = 0.5; % Relative infectiousness of subclinicals
par.tE = 1;     % Latent period
par.tI = 2.3;   % Infectious period


%------------- Specify Population Structure -------------
par.nAgeGroups = 16;

par.nSusComp = 14;
par.nVaxComp = 3;
par.nCaseComp = 3;
par.nHospComp = 5;
par.nDeathComp = 6;

% Load NZ population structure
fs = myDataPath + dataFileNames.popSizeFname; 
popSizeData = readmatrix(fs); 
popSizeDataBench = popSizeData;

% Fill entries with population distribution % Aggregate 75+ age-groups
par.popCount = [popSizeData(1:par.nAgeGroups-1, 2); sum(popSizeData(par.nAgeGroups:end, 2))];

par.totalPopSize = sum(par.popCount);
par.popDist = par.popCount/sum(par.popCount); 

% Northern region population size, only used to calculate hospitalisation 
% rates per 100,000 according to NR data
par.popSize_NR = 1937700;  


%--------------------- Seed Parameters --------------------
% Time window when community seed cases appear
par.seedDur = 7;
% Number of daily seed cases appearing in the community seeding window for 
% each age group
par.initialExp = 0.0001 * par.popCount;

% Date when border cases start getting seeded, approximately corresponding
% to the initial relaxation of border policies in 2022
par.borderTime = datenum("01-Mar-2022"); 
% Number of daily border seed cases appearing from borderTime. This is an
% approximation from the number of border arrivals in previous years
par.borderSeeds = 300;


%------------- Load Contact Matrix and Define NGM --------------
% Get Prem et al contact matrix from data folder
C = readmatrix(myDataPath + dataFileNames.CMdataFname); 

% Fill entries with population distribution % Aggregate 75+ age-groups
popCountBench = [popSizeDataBench(1:par.nAgeGroups-1, 2); ...
    sum(popSizeDataBench(par.nAgeGroups:end, 2))];
par.popDistBench = popCountBench/sum(popCountBench);

par.C_detBal = 0.5*(C + (par.popDistBench.')./(par.popDistBench) .* (C.'));


% ------------------ Population dynamics parameters ---------------------
% Demographic parameters birth, death, ageing - set all of these to zero to
% just have a static population
[Mu, b] = getDemogPars();
par.popnDeathRate = Mu;
par.popnBirthRate = b;
par.popnAgeingRate = 1/(5*365.25);


%------------- Disease Rate Data --------------

% Delta hazard ratios from Twohig et al & Fisman et al
HR_Hosp_Delta = 2.26;
OR_Death_Delta = 1;

% Omicron hazard ratios
HR_Hosp_Omi = 0.33;
HR_ICU_Omi = 0.3;
HR_Death_Omi = 0.3;

% Probability of developing symptoms (Fraser group)
par.pClin = [0.5440, 0.5550, 0.5770, 0.5985, 0.6195, 0.6395, 0.6585, ...
    0.6770, 0.6950, 0.7117, 0.7272, 0.7418, 0.7552, 0.7680, 0.7800, 0.8008]'; 

% Davies relative susceptibility
par.ui = [0.4000, 0.3950, 0.3850, 0.4825, 0.6875, 0.8075, 0.8425, ...
    0.8450, 0.8150, 0.8050, 0.8150, 0.8350, 0.8650, 0.8450, 0.7750, 0.7400]; 

% Get Herrera IHR and IFR and apply Delta and Omicron adjustments
[par.IHR0, ~, par.IFR0] = getHerreraRatesOmi(HR_Hosp_Delta, ...
    OR_Death_Delta, HR_Hosp_Omi, HR_ICU_Omi, HR_Death_Omi);


% ----- IHR adjustments -----

%%% Tested scaled IHR, scale factors obtained doing cumul actual daily
%%% hosp / cumul model daily hosp using the following model results file: 
% results/results_21MarPolicyChange_22MarRun_13FebVaxData/results_Uni_Filtered95_+0.0%_21Mar23_25-Feb-2023_fit.mat
% and with cumulative numbers between 25Jan22-18Mar23
IHRscale = [4.5450, 0.7335, 1.7772, 2.3990, 2.3895, 1.6952, ...
    1.2048, 0.8893, 0.6860, 0.5874, 0.4815, 0.4702, 0.4773, ...
    0.5401, 0.6539, 1.3063];

% Mean of the prior for the fitted IHR multiplier
IHRmultPriorMean = 0.5;

par.IHR0 = par.IHR0 .* IHRscale' .* IHRmultPriorMean;


% ----- IFR adjustments -----

% Ad hoc adjustment to IFR for oldest age group
par.IFR0(end) = 1.6 * par.IFR0(end);

% Mean of the prior for the fitted IFR multiplier
IFRmultPriorMean = 0.8;

par.IFR0 = par.IFR0 .* IFRmultPriorMean;


% --------------------- Testing and lag parameters -----------------------
par.tLatentToTest = 4;              % Days from onset of infectiousness to test
par.tTestToHosp = 1;                % days from test to hospital admission
par.tLOS = [2.0000 2.0000  2.0000  2.0000  2.0000  2.0000  2.6700 ...
    3.3400 4.0100 4.6800 5.3500 6.0200 6.6900 7.3600 8.0300 8.7000]';
par.tDeath = 14;                    % Days from admission to death 


% New way of determining CAR: setting fixed starting and ending CAR, linear
% % interpolation between the two, for 3 new age groups (0-30, 30-60, 60+)
par.pTest1_030 = 0.5; 
par.pTest1_3060 = 0.60;
par.pTest1_60p = 0.75;
par.pTest2_030 = 0.25; 
par.pTest2_3060 = 0.4;
par.pTest2_60p = 0.75;

% Dates between which the CAR decreases
CARchangeDates = datenum(["01MAY2022", "01JAN2023"]) - par.tBase(1);
CARchangeDays = CARchangeDates(2) - CARchangeDates(1) + 1;

par.pTestClin0 = ones(16, length(par.tBase));

% Before first date
par.pTestClin0(1:6, :) = par.pTest1_030;
par.pTestClin0(7:12, :) = par.pTest1_3060;
par.pTestClin0(13:16, :) = par.pTest1_60p;

% Linear decrease between dates
par.pTestClin0(1:6, CARchangeDates(1):CARchangeDates(2)) = repmat(linspace(par.pTest1_030, par.pTest2_030, CARchangeDays), 6, 1);
par.pTestClin0(7:12, CARchangeDates(1):CARchangeDates(2)) = repmat(linspace(par.pTest1_3060, par.pTest2_3060, CARchangeDays), 6, 1);
par.pTestClin0(13:16, CARchangeDates(1):CARchangeDates(2)) = repmat(linspace(par.pTest1_60p, par.pTest2_60p, CARchangeDays), 4, 1);

% After second date
par.pTestClin0(1:6, CARchangeDates(2):end) = par.pTest2_030;
par.pTestClin0(7:12, CARchangeDates(2):end) = par.pTest2_3060;
par.pTestClin0(13:16, CARchangeDates(2):end) = par.pTest2_60p;

% Scaling factor for testing probability of subclinical cases
par.subClinPtestMult = 0.4;

%-------------------------- Get vaccine data ----------------------------
par.vaccImmDelay = 14;  % delay in nb of days from vaccination to immunity
[~, par.doses1, par.doses2, par.doses3, par.doses4plus] = ...
    getVaccineData(myDataPath, dataFileNames.vaxDataFname, dataFileNames.vaxProjFname, par.vaccImmDelay, par.date0, par.tEnd);

smoothWindow = 56;
par.nDoses1Smoothed0 = [zeros(1, par.nAgeGroups); diff(smoothdata(par.doses1, 'movmean', smoothWindow));];
par.nDoses2Smoothed0 = [zeros(1, par.nAgeGroups); diff(smoothdata(par.doses2, 'movmean', smoothWindow))];
par.nDoses3Smoothed0 = [zeros(1, par.nAgeGroups); diff(smoothdata(par.doses3, 'movmean', smoothWindow))];
par.nDoses4Smoothed0 = [zeros(1, par.nAgeGroups); diff(smoothdata(par.doses4plus, 'movmean', smoothWindow))];


%---------------- Get antivirals data -------------------
therap_data = load(myDataPath + dataFileNames.AVdataFname);
movmean_period = 8 * 7; % Moving mean over 8 weeks
tailDays_toCut = 1 * 7; % Remove last 1 weeks of data to remove lagged entries
th_dates = therap_data.outTab.date(1:end-tailDays_toCut);

% Get number of daily treated, sum over the 75+, then smooth
daily_treated = therap_data.outTab.nTreated(1:end-tailDays_toCut, 1:15);
daily_treated(:, 16) = sum(therap_data.outTab.nTreated(1:end-tailDays_toCut, 16:end), 2);
daily_treated = movmean(daily_treated, movmean_period, 1);

% Get number of daily cases, sum over the 75+, then smooth
daily_cases = therap_data.outTab.nCases(1:end-tailDays_toCut, 1:15);
daily_cases(:, 16) = sum(therap_data.outTab.nCases(1:end-tailDays_toCut, 16:end), 2);
daily_cases = movmean(daily_cases, movmean_period, 1);

% Get ratio of smoothed daily treated and smoothed daily cases
daily_treatcaseratio = daily_treated ./ daily_cases;

% Add flat head and tail
par.daily_treatcaseratio = [zeros(datenum(th_dates(1))-par.date0, 16); 
    daily_treatcaseratio;
    repmat(daily_treatcaseratio(end, :), (par.date0+par.tEnd) - datenum(th_dates(end)), 1)];


%------------------------- Immunity parameters --------------------------
par.waneRateMean = 0.0045;          % Assumed mean daily waning rate
par.relRate_RtoS = 1.85;            % Relative rate of moving from R to S


% Run Khoury/Golding submodel to generate immunity parameters
kLog = 2.94/log(10);               % steepness of logistic relationship betweem log titre and VE -  2.94/log(10) -  Khoury Table S5
no50_sympt = log(0.2);         %  determines mapping from titre to VE symptoms - log(0.2)  match to NG 
no50_sev = log(0.03);            % determines mapping from titre to VE severe - log(0.04)  for match with NG results   - Khoury Table S5 is a bit more optimistic: no50_sev = log(0.03);           % Khoury offset (Table S5)
logTitreRatio = log(10);               % ratio of titre from one compartment to next (don't change this without also changing wanin rate and calibrating to Golding reulsts)

logTitre0_2 = log(0.2);                 % 0.2  strength of immunity (measured as initial titre) after 2 doses
logTitre0_3 = log(0.4);                 % 0.4  strength of immunity (measured as initial titre) after 3 doses
logTitre0_inf = log(0.8)+log(5);               % 0.8    strength of immunity (measured as initial titre) after 0/1 doses + infection
logTitre0_inf_plus2 = log(3)+log(5);           % 3       strength of immunity (measured as initial titre) after 2 doses + infection
logTitre0_inf_plus3 = log(7)+log(5);          % 7      strength of immunity (measured as initial titre) after 3 doses + infection

% Immunity to hospitalisation and death (from vaccine or prior infection) 
% cannot wane below this 
minVEsev = 0.5;  

% set titre levels for each susceptible compartment 
logTitreSequence = logTitreRatio*[0, -1, -2, -3];
logTitreLevels = [-inf -inf    logTitre0_2 + logTitreSequence   logTitre0_3 + logTitreSequence     logTitre0_inf_plus3 + logTitreSequence]; 

% convert titre levels to immunity agaist each outcome
par.VEi = 1./(1+ exp(-kLog*(logTitreLevels-no50_sympt))); % get Khoury VE
par.VEh = 1./(1+ exp(-kLog*(logTitreLevels-no50_sev))); % get Khoury VE
par.VEt = zeros(1, par.nSusComp);
par.VEs = par.VEi;
par.VEf = par.VEh;

% apply minimum immunity constraint to severe outcomes
par.VEh(3:end) = max(minVEsev, par.VEh(3:end));
par.VEf(3:end) = max(minVEsev, par.VEf(3:end));

% Calculate the proporton of post-infection indiciduala with 0/1 or 2 doses
% who go to each of the 4 post-infection susceptible compartments
% this is done by solving an ODE to calculate the proportion of a
% post-infection cohort who are in each compartment at time t such that
% their average titre has dropped by the specified amount
logTitreDrop =  [logTitre0_inf, logTitre0_inf_plus2] - logTitre0_inf_plus3;
Y0 = getImmPars(logTitreRatio, logTitreDrop);
Y0 = Y0./sum(Y0, 2);

% Proportion going into W1, W2, W3, W4 post recovery
postRecovDist_Unvaxed = Y0(1, :);
postRecovDist_Vaxed2 = Y0(2, :);
postRecovDist_Vaxed3 = [1 0 0 0];

[par.waneNet_StoS, par.waneNet_RtoS, par.vaxNet] = ...
    getATCmatrices(postRecovDist_Unvaxed, postRecovDist_Vaxed2, postRecovDist_Vaxed3);


% --------------- Contact matrix adjustment parameters ------------------
% divide the contact matrix up into blocks - this vector specifies the 
% number of 5-year age classes in each "block"
par.ageBlockSizes = [3 2 2 3 2 4];      
par.nAgeBlocks = length(par.ageBlockSizes);

% weighting matrix for initial contact matrix
par.Cw1 = [1.1 0.7  0.55 0.45 0.45 0.5;  
      0    1.2  0.7  0.5  0.5  0.3;
      0    0    1.1  0.5  0.5  0.5;
      0    0    0    0.15 0.15 0.45;
      0    0    0    0    0.15 0.45;
      0    0    0    0    0    0.15];

% ------------------------ VOC model ----------------------------------
% Date of arrival of the new variant        
par.vocWaneDate = datenum('20-Jun-2022');
% Time window over which new variant becomes predominant
par.vocWaneWindow = 2;    

% ------------------------ VOC VE model ----------------------------------
% To model reduced VE for VOC set VOC_titreDrop < 1 
par.VOC_logTitreDrop = log(0.4);
logTitre0_2_VOC = par.VOC_logTitreDrop + logTitre0_2;                 
logTitre0_3_VOC = par.VOC_logTitreDrop + logTitre0_3;
logTitreLevels_VOC = [-inf -inf logTitre0_2_VOC+logTitreSequence ...
    logTitre0_3_VOC+logTitreSequence logTitre0_inf_plus3+logTitreSequence]; 

% Get Khoury VE
par.VEi_VOC = 1./(1+ exp(-kLog*(logTitreLevels_VOC-no50_sympt)));
par.VEh_VOC = 1./(1+ exp(-kLog*(logTitreLevels_VOC-no50_sev)));

par.VEt_VOC = zeros(1, par.nSusComp);
par.VEs_VOC = par.VEi_VOC;
par.VEf_VOC = par.VEh_VOC;

% apply minimum immunity constraint to severe outcomes
par.VEh_VOC(3:end) = max(minVEsev, par.VEh_VOC(3:end));
par.VEf_VOC(3:end) = max(minVEsev, par.VEf_VOC(3:end));


%---------------- Scenarios set-up -------------------

% Date of policy change
policyDate = '15-May-2023';

% Transmission (Ct) increase multiplier for policy changes
% Each column is a policy, rows are the lower bound, middle value, and 
% upper bound of the transmission increase for each policy
% par.cRampDeltaPolicy = [1.0742, 1.0126, 1.0013, 1, 0.9859, 0.9678, 1.05, 1.1;
%                     1.1478, 1.0252, 1.0027, 1, 0.9906, 0.9785, 1.075, 1.125;
%                     1.2213, 1.0377, 1.0040, 1, 0.9953, 0.9892, 1.1, 1.15];
% % For merged results files, only used if cRampDeltaPolicy has multiple lines (levels):
% policyNames = ["Policy1", "Policy2", "Policy3", "Policy4", "Policy5", ...
%     "Policy6", "Policy7", "Policy8"]; 

% % Example with two policies giving two transmission increases each:
% par.cRampDeltaPolicy = [1.1, 1.2];
% % Example with two policies and 2 levels (lower-upper) each
% par.cRampDeltaPolicy = [1.0, 1.1; 1.2, 1.3];
% policyNames = ["Policy1", "Policy2"];

par.cRampDeltaPolicy = [1];


% Dates when VOC becomes predominant
VOC2date = '15-Nov-2022';
% Boolean to switch on (=1) or off (=0) second variant arriving at date specified above
VOC2active = [1];

% Seasonality multiplier applied to Ct from 1Apr23
seasonMultAmp = [0];

% Effect of antivirals on IHR (0-no effect; 1-full protection)
antiviralsEffectIHRmult = [0];

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

% Reduction multiplier applied to vaccination coverage (0 - no reduction, 1 - all doses removed)


redFactor = [0, 1, 0, 1, 1, 0.1, 0, 0, 0];
% par.dosesToRed = 3:4; % Eg 3:4 will apply the reduction factor to 3+ doses
% par.agesToRed = 1:16; % Eg 1:16 will apply the reduction factor to all age groups
dosesToRed = repmat({1:4}, numel(redFactor), 1); % Eg 3:4 will apply the reduction factor to 3+ doses (cell array for each scenario)
agesToRed = {1:16; 1:16; 1:16; 1:16; 1:12; 1:16; 1:16; 1:16; 1:16}; % Eg 1:16 will apply the reduction factor to all age groups  (cell array for each scenario)

% Table with combination of scenarios
par.scenarios = array2table(combvec(reshape(par.cRampDeltaPolicy',1,[]), ...
    VOC2active, seasonMultAmp, antiviralsEffectIHRmult, redFactor)');
par.scenarios.Properties.VariableNames = {'cRampDeltaPolicy', 'VOC2active', ...
    'seasonMultAmp', 'antiviralsEffectIHRmult', 'vaxCovRedFactor'};
par.nScenarios = size(par.scenarios, 1);
par.scenarios.dosesToRed = dosesToRed;
par.scenarios.agesToRed = agesToRed;


% Scenario names for output files
par.scenario_names = strings(par.nScenarios, 1);
for is = 1:par.nScenarios
%     scenario_names(is) = sprintf("%.1faviralEff", ...
%         par.scenarios(is, :).antiviralsEffectIFRmult);
    par.scenario_names(is) = sprintf("scenario%i", is);
end

par.scenarios.policyDate = repmat(policyDate, par.nScenarios, 1);
par.scenarios.VOC2date = repmat(VOC2date, par.nScenarios, 1);
par.scenarios.scenarioNumber = (1:par.nScenarios)';

% flag for sensitivity run

if dataFileNames.popSizeFname == 'popproj_national2018-21.xlsx'
    par.sensitivity_flag = 1; % make 1
else
    par.sensitivity_flag = 0; % normally 0
end

end
