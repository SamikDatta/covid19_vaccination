function epiVarsCompact = extractEpiVarsCompact_vaxSplitInf(t, Y, par)

% Interpolate between original and VOC VE parameters
Alpha = normcdf( t'-par.vocWaneDate(1) , 0, par.vocWaneWindow);
VEi = (1-Alpha).*par.VEi + Alpha.*par.VEi_VOC;

% Take a matrix representing the ODE solution over time and restructure as a set of matrices representing epidemiological
% variables
% The number of rows in the input matrix Y is the number of time steps
% The function returns a structure with fields which are matrices each of which has 16 columns corresponding to 16 age groups and 
% rows representing time points. The fields are:
% N - total population size in each group
% V1, V2, V3 - number with at least 1, 2, 3 doses
% S - number in a susceptible compartment
% Sw  - weighted susceptibility (weighted by vulnerability to infection
% 1-VEi in each susceptibility class)
% E1 - number in an exposed compartment for 1st infection (S-class 1-10)
% E2 - number in an exposed compartment for re-infectoin (S-class 11-14)
% I - number infectious & clinical
% A - number infection & subclinical
% Ci - cumulative cases in immunity categorry i (i = 0, 1, 2, 3 doses, r = reinfection)
% Hiocc - current hospital occupancy
% Hidis - cumulative hospital discharges
% Fi - cumulative fatalities


[nDays, ~] = size(Y);

iCount = 1;
nC = par.nAgeGroups;
epiVarsCompact.N = Y(:, iCount : iCount+nC-1);

iCount = iCount+nC;
epiVarsCompact.V1 = Y(:, iCount : iCount+nC-1 );
iCount = iCount+nC;
epiVarsCompact.V2 = Y(:, iCount : iCount+nC-1 );
iCount = iCount+nC;
epiVarsCompact.V3 = Y(:, iCount : iCount+nC-1 );

iCount = iCount+nC;
nC = par.nAgeGroups * par.nSusComp;
Stemp = Y(:, iCount : iCount+nC-1);         % matrix of nDays x 224 whose columns are the 16 x 14 = 224 susceptible compartments for the 16 age groups and 14 susceptible comppartments
epiVarsCompact.S = condenseSusLevels(Stemp, par);       % sum over susceptible compartments to get a nDays x 16 matrix of total susceptible in each age groups
epiVarsCompact.Sw = condenseSusLevels(Stemp.*repelem(1-VEi, 1, par.nAgeGroups) , par);      % weighted by (1-immunity to infection) to get average susceptibility to infection in each age group

iCount = iCount+nC;
E1temp = Y(:, iCount:iCount+(par.nAgeGroups*10)-1 );      % first infections exposed class
E2temp = Y(:,  (iCount+par.nAgeGroups*10):iCount+nC-1);     % reinfections exposed class
Ev01temp = Y(:, iCount:iCount+(par.nAgeGroups*2)-1 );   % first infection exposed unvax or 1dose
Ev2temp = Y(:, (iCount+par.nAgeGroups*2):(iCount+par.nAgeGroups*6)-1 );   % first infection exposed 2dose
Ev3temp = Y(:, (iCount+par.nAgeGroups*6):(iCount+par.nAgeGroups*10)-1 );   % first infection exposed 3dose
epiVarsCompact.E1 = condenseSusLevels(E1temp, par);
epiVarsCompact.E2 = condenseSusLevels(E2temp, par);
epiVarsCompact.Ev01 = condenseSusLevels(Ev01temp, par);
epiVarsCompact.Ev2 = condenseSusLevels(Ev2temp, par);
epiVarsCompact.Ev3 = condenseSusLevels(Ev3temp, par);

iCount = iCount+nC;
Itemp = Y(:, iCount : iCount+nC-1);
epiVarsCompact.I = condenseSusLevels(Itemp, par);

iCount = iCount+nC;
Atemp = Y(:, iCount : iCount+nC-1);
epiVarsCompact.A = condenseSusLevels(Atemp, par);

iCount = iCount+nC;
nC = par.nAgeGroups * (par.nSusComp-1);
% skipping R

iCount = iCount+nC;
nC = par.nAgeGroups * par.nCaseComp;
epiVarsCompact.C0 = Y(:, (iCount + (par.nCaseComp-1)*par.nAgeGroups):(iCount+nC-1) );       % skipping two unobserved case compartments
iCount = iCount+nC;
epiVarsCompact.C1 = Y(:, (iCount + (par.nCaseComp-1)*par.nAgeGroups):(iCount+nC-1) );       % skipping two unobserved case compartments
iCount = iCount+nC;
epiVarsCompact.C2 = Y(:, (iCount + (par.nCaseComp-1)*par.nAgeGroups):(iCount+nC-1) );       % skipping two unobserved case compartments
iCount = iCount+nC;
epiVarsCompact.C3 = Y(:, (iCount + (par.nCaseComp-1)*par.nAgeGroups):(iCount+nC-1) );       % skipping two unobserved case compartments
iCount = iCount+nC;
epiVarsCompact.Cr = Y(:, (iCount + (par.nCaseComp-1)*par.nAgeGroups):(iCount+nC-1) );       % skipping two unobserved case compartments

iCount = iCount+nC;
nC = par.nAgeGroups * par.nHospComp;
epiVarsCompact.H0occ = Y(:, (iCount + (par.nHospComp-2)*par.nAgeGroups):(iCount+(par.nHospComp-1)*par.nAgeGroups-1));    % skipping three unobserved hosp compartments
epiVarsCompact.H0dis = Y(:, (iCount + (par.nHospComp-1)*par.nAgeGroups):(iCount+nC-1));    
iCount = iCount+nC;
epiVarsCompact.H1occ = Y(:, (iCount + (par.nHospComp-2)*par.nAgeGroups):(iCount+(par.nHospComp-1)*par.nAgeGroups-1));    % skipping three unobserved hosp compartments
epiVarsCompact.H1dis = Y(:, (iCount + (par.nHospComp-1)*par.nAgeGroups):(iCount+nC-1));    
iCount = iCount+nC;
epiVarsCompact.H2occ = Y(:, (iCount + (par.nHospComp-2)*par.nAgeGroups):(iCount+(par.nHospComp-1)*par.nAgeGroups-1));    % skipping three unobserved hosp compartments
epiVarsCompact.H2dis = Y(:, (iCount + (par.nHospComp-1)*par.nAgeGroups):(iCount+nC-1));    
iCount = iCount+nC;
epiVarsCompact.H3occ = Y(:, (iCount + (par.nHospComp-2)*par.nAgeGroups):(iCount+(par.nHospComp-1)*par.nAgeGroups-1));    % skipping three unobserved hosp compartments
epiVarsCompact.H3dis = Y(:, (iCount + (par.nHospComp-1)*par.nAgeGroups):(iCount+nC-1));    
iCount = iCount+nC;
epiVarsCompact.Hrocc = Y(:, (iCount + (par.nHospComp-2)*par.nAgeGroups):(iCount+(par.nHospComp-1)*par.nAgeGroups-1));    % skipping three unobserved hosp compartments
epiVarsCompact.Hrdis = Y(:, (iCount + (par.nHospComp-1)*par.nAgeGroups):(iCount+nC-1));    

iCount = iCount+nC;
nC = par.nAgeGroups * par.nDeathComp;
epiVarsCompact.F0 = Y(:, (iCount + (par.nDeathComp-1)*par.nAgeGroups):(iCount+nC-1));       % skipping 1 unobserved death compartment
iCount = iCount+nC;
epiVarsCompact.F1 = Y(:, (iCount + (par.nDeathComp-1)*par.nAgeGroups):(iCount+nC-1));       % skipping 1 unobserved death compartment
iCount = iCount+nC;
epiVarsCompact.F2 = Y(:, (iCount + (par.nDeathComp-1)*par.nAgeGroups):(iCount+nC-1));       % skipping 1 unobserved death compartment
iCount = iCount+nC;
epiVarsCompact.F3 = Y(:, (iCount + (par.nDeathComp-1)*par.nAgeGroups):(iCount+nC-1));       % skipping 1 unobserved death compartment
iCount = iCount+nC;
epiVarsCompact.Fr = Y(:, (iCount + (par.nDeathComp-1)*par.nAgeGroups):(iCount+nC-1));       % skipping 1 unobserved death compartment



