function dYdt = calcPopnDynamics(Y, b, Mu, r)

[~, nCols] = size(Y);

dYdt = [b; r*Y(1:end-1, :)] - [r*Y(1:end-1, :); zeros(1, nCols)] - Mu.*Y;
