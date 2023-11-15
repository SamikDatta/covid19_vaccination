function [ N, V, S, E, I, A, R, C0, C1, C2, C3, Cr, H0, H1, H2, H3, Hr, F0, F1, F2, F3, Fr] = extractEpiVarsFull(y, par)

% Take a column vector representing the ODE state at a single time point
% and restructure as a set of matrices representing epidemiological
% variables
% Each matrix has 16 rows representing 16 age groups
% Infection states (S, E, I, A, R) have 14 columns corresponding to the 14
% susceptible compartments
% N is a single column (total population size in each age group)
% V has three columns for number of people with at least 1, 2, 3 doses
% Ci, Hi and Fi have a specified number of columns corresponding to the number
% of unobserved+observed states for each of cases, hospitalisations and
% fatalities


iCount = 1;
nC = par.nAgeGroups;
N = reshape( y(iCount : iCount+nC-1), par.nAgeGroups, 1);

iCount = iCount+nC;
nC = par.nAgeGroups * par.nVaxComp;
V = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nVaxComp);      % column j is the number in each age group who have had at least j doses

iCount = iCount+nC;
nC = par.nAgeGroups * par.nSusComp;
S = reshape( y(        iCount : iCount+nC-1), par.nAgeGroups, par.nSusComp);

iCount = iCount+nC;
E = reshape( y(  iCount : iCount+nC-1), par.nAgeGroups, par.nSusComp);

iCount = iCount+nC;
I = reshape( y( iCount : iCount+nC-1), par.nAgeGroups, par.nSusComp);

iCount = iCount+nC;
A = reshape( y( iCount : iCount+nC-1), par.nAgeGroups, par.nSusComp);

iCount = iCount+nC;
nC = par.nAgeGroups * (par.nSusComp-1);
Rpart = reshape( y( iCount : iCount+nC-1), par.nAgeGroups, par.nSusComp-1);
R = [Rpart, max(0, N - sum(S+E+I+A, 2) - sum(Rpart, 2) )];

iCount = iCount+nC;
nC = par.nAgeGroups * par.nCaseComp;
C0 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nCaseComp);
iCount = iCount+nC;
C1 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nCaseComp);
iCount = iCount+nC;
C2 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nCaseComp);
iCount = iCount+nC;
C3 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nCaseComp);
iCount = iCount+nC;
Cr = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nCaseComp);


iCount = iCount+nC;
nC = par.nAgeGroups * par.nHospComp;
H0 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nHospComp);
iCount = iCount+nC;
H1 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nHospComp);
iCount = iCount+nC;
H2 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nHospComp);
iCount = iCount+nC;
H3 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nHospComp);
iCount = iCount+nC;
Hr = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nHospComp);


iCount = iCount+nC;
nC = par.nAgeGroups * par.nDeathComp;
F0 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nDeathComp);
iCount = iCount+nC;
F1 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nDeathComp);
iCount = iCount+nC;
F2 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nDeathComp);
iCount = iCount+nC;
F3 = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nDeathComp);
iCount = iCount+nC;
Fr = reshape( y( iCount : iCount+nC-1 ), par.nAgeGroups, par.nDeathComp);


