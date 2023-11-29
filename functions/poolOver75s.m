function HSU16 = poolOver75s(HSU, ethNames)

lastAgeGroup = 16;

HSU16 = HSU(1:lastAgeGroup, :);
nEth = length(ethNames);
for iEth = 1:nEth
    HSU16.(ethNames(iEth))(lastAgeGroup) = sum(HSU.(ethNames(iEth))(lastAgeGroup:end) );
end

