function Y = combineBands(X)


Y  = X(:, 1:2:end, :);
Y2 = X(:, 2:2:end, :);
nCols = size(Y2, 2);
Y(:, 1:nCols, :) = Y(:, 1:nCols, :) + Y2;
