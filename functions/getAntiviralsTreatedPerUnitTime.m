function [antiviralsTreatedPropNow] = getAntiviralsTreatedPerUnitTime(t, date0, ...
    tEnd, dailyTreatedProp)
% Function that takes the timeseries vector of age-split testing
% probabilities and returns the current value
% INPUTS:
% - t: current time
% - date0: start of simulation, defined in getPar
% - tEnd: end of simulation, defined in getPar
% - dailyTreatedProp: timeseries age-split vector of proportion of daily 
%                     treated cases
% OUTPUTS:
% - antiviralsTreatedNow: current age-split vector of treated cases

ta = date0:date0+tEnd;
antiviralsTreatedPropNow = zeros(16, 1);
for ag = 1:16
    antiviralsTreatedPropNow(ag) = interp1(ta, dailyTreatedProp(:, ag), t)';
end

end
