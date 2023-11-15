function posterior = getPosteriorStats(pUniFiltered100, dUniFiltered100, parsToFit, index)
% Author: GV, Dec 2022
% Function that outputs structure of parameter values and corresponding 
% distances given the normalised vector of fitted parameters and the 
% indexes of interest
% INPUTS:
% - pUniFiltered100: matrix of 1% best parameter sets (normalised)
% - dUniFiltered100: vector of distances corresponding to each pUni set
% - parsToFit: cell vector of fitted parameter names
% - indexes: indexes of pUniFiltered set of interest
% OUTPUT:
% - posterior: actual values of retained parameters and corresp. distances

% Uniform perturbation between +/- r, z is a random deviate [0,1]
plusMinus = @(z, r) (2 * r * z - r);
% Uniform perturbation on integers between +/- r, z is a random deviate [0,1]
plusMinusInt = @(z, r) (floor((2 * r + 1) * z) - r);

posterior = table();

% Table with best fit value and quartile values within retained range
ThetaAll = array2table([pUniFiltered100(index, :); quantile(pUniFiltered100, ...
    [0.25, 0.75], 1)], 'VariableNames', parsToFit);

for i = 1:size(ThetaAll, 1)

    Theta = ThetaAll(i, :);

    posterior.dateSeed(i) = datetime(datenum('19JAN2022') + ...
        plusMinusInt(Theta.dateSeed, 3), 'ConvertFrom', 'datenum');

    posterior.Ct(i) = 0.68 + plusMinus(Theta.Cstart, 0.1);

    posterior.CtRampStart1(i) = datetime(datenum('10-Mar-2022') + ...
        plusMinusInt(Theta.rampStart, 5), 'ConvertFrom', 'datenum');
    posterior.CtRampDays1(i) = 55 + plusMinusInt(Theta.rampDays, 20);
    posterior.CtRamp1(i) = 1.1 + plusMinus(Theta.Cramp, 0.21);

    posterior.CtRampStart2(i) = datetime(datenum('15-Sep-2022') + ...
        plusMinusInt(Theta.ramp2Start, 5), 'ConvertFrom', 'datenum');
    posterior.CtRampDays2(i) = 10 + plusMinusInt(Theta.ramp2Days, 9);
    posterior.CtRamp2(i) = (1.2 + plusMinus(Theta.Cramp2, 0.1));

    posterior.pTestMult(i) = 1 + plusMinus(Theta.pTestMult, 0.2);

    posterior.IFRmult(i) = 0.8 * (1 + plusMinus(Theta.IFR, 0.5));
    posterior.IHRmult(i) = 0.5 * (1 + plusMinus(Theta.IHR, 0.5));

    posterior.relaxContactsAlpha(i) = 0.4 + plusMinus(Theta.relaxAlpha, 0.4);
    posterior.changeWindow(i) = 70 + plusMinusInt(Theta.MRampDays, 20);

    posterior.vocWaneAmount(i) = 0.4 + plusMinus(Theta.vocWane, 0.3) ;
    posterior.vocWaneAmount2(i) = 0.5 + plusMinus(Theta.vocWane, 0.3) ;
    posterior.waneRateMult(i) = (1 + plusMinus(Theta.waneRate, 0.5)) * 0.009*0.5;               

end


end