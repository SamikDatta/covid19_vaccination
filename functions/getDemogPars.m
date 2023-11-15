function [Mu, b] = getDemogPars()

      
% StatsNZ infoshare death rate per 1000 for 2019 for 0-75 year age groups
% Death rate for 75+ chosen to giver a similar equilibrium distribution to
% current population
Mu = 1/1000/365.25 * [1.07 0.08 0.17 0.41 0.6 0.56 0.73 0.83 1.21 1.95 3.07 4.45 6.49 10.27 16.69 136 ]';


% StatsNZ infoshare total births 2019
b = 59637/365.25;     

