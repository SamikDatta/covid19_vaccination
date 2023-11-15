function M = condenseSusLevels(A, par)

% Takes a n x m matrix whose n rows are time points and whose m = k x 16 columns correspond
% to the k categories (e.g. different susceptible compartments or groups of
% susceptibility classes) and 16 age groups.
% Columns are in order: [class1, age1], [class1, age2], ... [class1, age16], [class2, age1], ...  
% Returns a n x 16 matrix whose columns are summed over the k categories 

[nDays, nCols] = size(A);
nToSum = nCols/par.nAgeGroups;
 
% counter = 0:par.nAgeGroups:(par.nAgeGroups*(nToSum-1)); % vector of indices to pick out columns correspondin to the 14 susceptible compartments for the 1st age group  
% 
% M = zeros(nDays, par.nAgeGroups);
% for iAge = 1:par.nAgeGroups
%   M(:, iAge) = sum(A(:, iAge+counter), 2); 
% end


M = sum(reshape(A, nDays, par.nAgeGroups, nToSum), 3);


