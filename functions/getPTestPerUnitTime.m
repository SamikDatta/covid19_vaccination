function [pTestClinNow, pTestSubNow] = getPTestPerUnitTime(t, date0, ...
    tEnd, pTestClin, pTestSub)
% Function that takes the timeseries vector of age-split testing
% probabilities and returns the current value
% INPUTS:
% - t: current time
% - date0: start of simulation, defined in getPar
% - tEnd: end of simulation, defined in getPar
% - pTestClin: timeseries age-split vector of pTest for clinical cases
% - pTestSub: timeseries age-split vector of pTest for subclinical cases
% OUTPUTS:
% - pTestClinNow: current age-split vector of pTest for clinical cases
% - pTestSubNow: current age-split vector of pTest for subclinical cases

ta = date0:date0+tEnd;
[pTestClinNow, pTestSubNow] = deal(zeros(16, 1));
for ag = 1:16
    pTestClinNow(ag) = interp1(ta, pTestClin(ag, :), t)';
    pTestSubNow(ag) = interp1(ta, pTestSub(ag, :), t)';
end

end
