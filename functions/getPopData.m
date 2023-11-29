function [popData, ethNames] = getPopData(fName);

opts = detectImportOptions(fName);
opts = setvartype(opts, {'Age', 'Asian', 'EuropeanorOther', 'Maori', 'PacificPeoples', 'Various' , 'Total'}, 'double');
ethNames = string(opts.VariableNames(2:end));

popData = readtable(fName, opts);

% exclude 'Total' row from data table
popData = popData(1:end-1, :);
