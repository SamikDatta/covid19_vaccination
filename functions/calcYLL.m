function YLL = calcYLL(age, exM, exF, nDeathsActualM, nDeathsActualF, ageModel, nDeathsModel)

% Calculate CIs for YLL from life expectancy tables (as returned by importLifeTables()), ensemble of model outputs on deaths
% in 5-year age bands (nDeathsModel), and empirical age distribution of actual deaths (which is used to distribute model deaths within each age band)
%
% age is as returned by importLifeTables()
% exM and exF are the middle column (central estimate) of the corresponding matrices returned by importLifeTables()
% nDeathsActualM and nDeathsActualF should be vectors of the same length as 'age' containing actual male and female deaths in one year age bands
% ageModel is a vector of length n containng the lower age breaks of the n model age groups (n=16 usually)
% nDeathsModel should be n x m matrix whose m columns are the number of deaths in each age group for each of m model realisations 
%
% Returns YLL as a 1 x m vector containing the total YLL for each of the m model realisations

nAges = length(age);
[nAgesModel, nSims] = size(nDeathsModel);

nDeathsDistM = zeros(nAges, nSims);      % vector for model deaths, distributed to one year age bands
nDeathsDistF = zeros(nAges, nSims);      % vector for model deaths, distributed to one year age bands
for iAge = 1:nAgesModel
    % indices for actual deaths that correspond to this model age band
    if iAge < nAgesModel
        jInBand = age >= ageModel(iAge) & age < ageModel(iAge+1);
    else
        jInBand = age >= ageModel(iAge);
    end

    % Within this model age band, calculate proportion of actual deaths in
    % each one year age group and sex
   if sum(nDeathsActualM(jInBand)+nDeathsActualF(jInBand)) > 0
        pM = nDeathsActualM(jInBand)/sum(nDeathsActualM(jInBand)+nDeathsActualF(jInBand));
        pF = nDeathsActualF(jInBand)/sum(nDeathsActualM(jInBand)+nDeathsActualF(jInBand));
   else    % if no actual deaths in this model age band, distribute equally
        pM = ones(sum(jInBand), 1)/(2*sum(jInBand));
        pF = pM;
    end
    
    nDeathsDistM(jInBand, :) = nDeathsModel(iAge, :).*pM;
    nDeathsDistF(jInBand, :) = nDeathsModel(iAge, :).*pF;
end

YLLbyAgeM = nDeathsDistM.*exM;
YLLbyAgeF = nDeathsDistF.*exF;

% figure;
% plot(age, quantile(YLLbyAgeM, [0.025, 0.975], 2), age, quantile(YLLbyAgeF, [0.025, 0.975], 2))
% legend('M lower', 'M upper', 'F lower', 'F upper')

YLL = sum(YLLbyAgeM+YLLbyAgeF);       



