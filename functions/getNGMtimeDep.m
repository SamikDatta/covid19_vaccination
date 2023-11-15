function [NGM, NGMclin] = getNGMtimeDep(t, par)

% Returns unitary NGM (dominant eigenvalue = 1), multiply by REI for actual
% NGM in a fully susceptible population

% Alpha is 0 before the change in contact matrix, 1 after the change has
% completed and increased from 0 to 1 linearly with time during the change window 
Alpha = max(0, min(1, (t-par.contactPar.changeDate)/par.contactPar.changeWindow));
% Current value of contact matrix weights
weightsNow = (1-Alpha)*par.contactPar.weights + Alpha*par.contactPar.weightsChange;

% Apply weights to contact matrix
C_detBalNow = par.C_detBal .* weightsNow;

NGM0 = diag(par.ui)*C_detBalNow'*diag(par.pClin + par.cSub*(1-par.pClin)); % Construct unnormalised NGM"
u = 1/max(abs(eig(NGM0)));                                                  % Normalising factor to make dominant eigenvalue = 1

% Re-adjust for actual population distribution being modelling as per Prem
% et al
C_popAdj = (par.popDist./par.popDistBench).' .* C_detBalNow;

NGM =     u * diag(par.ui) * C_popAdj' * diag(par.pClin + par.cSub*(1-par.pClin));    % Set final NGM
NGMclin = u * diag(par.ui) * C_popAdj';                                           % NGM for clinical individuals (as used in simulation models)





