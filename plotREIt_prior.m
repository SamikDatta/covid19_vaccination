clear
close all

nSample = 10000;
nToPlot = 10;

par.date0 = datenum('01JAN2022'); 
par.tEnd = datenum('30JUN2023')-par.date0;

par.R0 = 3.25;


% Functions to add a zero-centred random variable to selected parameters
% according to random deviates [0,1] specified in input parameter Theta
plusMinus = @(z, r)(2*r*z - r);                     % uniform perturbation between +/- r, z is a random deviate [0,1]
plusMinusInt = @(z, r)(floor((2*r+1)*z) - r );      % uniform perturbation on integers between +/- r, z is a random deviate [0,1]

t = par.date0:par.date0+par.tEnd;

REI = nan(nSample, par.tEnd+1);

for iSample = 1:nSample
    
    Theta.Cstart = rand;
    Theta.rampStart = rand;
    Theta.ramp2Start = rand;
    Theta.rampDays = rand;
    Theta.ramp2Days = rand;
    Theta.Cramp = rand;
    Theta.Cramp2 = rand;
    
    
    
    % Starting value of Ct (fitted) - between 0.58-0.78
    par.Ct = (0.68 + plusMinus(Theta.Cstart, 0.1)) .* ones(1, par.tEnd+1);
    
    CtRampStarts = [datenum('10-Mar-2022') + plusMinusInt(Theta.rampStart, 5), ...
        datenum('15-Sep-2022') + plusMinusInt(Theta.ramp2Start, 5) ];
    CtRampDays = [55 + plusMinusInt(Theta.rampDays, 20), ...
        10 + plusMinusInt(Theta.ramp2Days, 9) ];
    
    % Ct after 1st ramp up (0.89-1.31) 
    % and after 2nd ramp up (0.89-1.31 * 1.1-1.3)
    CtRamp = [1.1 + plusMinus(Theta.Cramp, 0.21), ...
        (1.1 + plusMinus(Theta.Cramp, 0.21)) * (1.2 + plusMinus(Theta.Cramp2, 0.1))];
    
    
    % Adding Ct ramp-ups at each date:
    for pci = 1:length(CtRampStarts)
        ti = datenum(CtRampStarts(pci)) - par.date0;
        par.Ct(ti:ti+CtRampDays(pci)-1) = linspace(par.Ct(ti), CtRamp(pci), CtRampDays(pci));
        par.Ct(ti+CtRampDays(pci):end) = CtRamp(pci);
    %    par.Ct = par.Ct(1:parBase.tEnd+1); % Making sure size stays the same
    end
    
    
    REI(iSample, :) = par.R0*par.Ct;

end

Q = quantile(REI, [0.025, 0.975]);
Rlower = Q(1, :);
Rupper = Q(2, :);

dates = datetime(t, 'ConvertFrom', 'datenum');
darkGrey = [0.4 0.4 0.4];
lightGrey = [0.9 0.9 0.9];
figure(1);
fill( [dates, fliplr(dates) ], [Rlower, fliplr(Rupper)], lightGrey, 'FaceAlpha', 0.5, 'linestyle', 'none' );  
hold on
plot(dates, REI(1:nToPlot, :), 'color', darkGrey);
ylabel('R_{EI}(t)')
ylim([0, 5.5])
grid on
