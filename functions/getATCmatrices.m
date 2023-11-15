function [waneNet_StoS, waneNet_RtoS, vaxNet] = getATCmatrices(postRecovDist_Unvaxed, postRecovDist_Vaxed2, postRecovDist_Vaxed3)

% Calculate "ATC" matrices:
% - waneNet_StoS specifies fluxes between susceptible compartments
% as a result of waning
% - waneNet_RtoS specifies which S compartments individuals move to from which R compartments
% - vaxNet specified fluxes between susceptible compartments as a result of
% receiving a vaccine dose

waneNet_StoS = [0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                0  0 -1  1  0  0  0  0  0  0  0  0  0  0;
                0  0  0 -1  1  0  0  0  0  0  0  0  0  0;
                0  0  0  0 -1  1  0  0  0  0  0  0  0  0;
                0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                0  0  0  0  0  0 -1  1  0  0  0  0  0  0;
                0  0  0  0  0  0  0 -1  1  0  0  0  0  0;
                0  0  0  0  0  0  0  0 -1  1  0  0  0  0;
                0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                0  0  0  0  0  0  0  0  0  0 -1  1  0  0;
                0  0  0  0  0  0  0  0  0  0  0 -1  1  0;
                0  0  0  0  0  0  0  0  0  0  0  0 -1  1;
                0  0  0  0  0  0  0  0  0  0  0  0  0  0];
                
                
waneNet_RtoS = [zeros(14, 10), [repmat(postRecovDist_Unvaxed, 2, 1); repmat(postRecovDist_Vaxed2, 4, 1); repmat(postRecovDist_Vaxed3, 4, 1); [ones(4, 1), zeros(4, 3)]] ];


vaxNet = zeros(14);  
vaxNet(1, 1) = -1;  vaxNet(1, 2) = 1;  
vaxNet(2, 2) = -1;  vaxNet(2, 3) = 1;  
vaxNet(3:6, 3:6) = -eye(4);  vaxNet(3:6, 7) = 1;  
vaxNet(8:10, 8:10) = -eye(3);  vaxNet(8:10, 7) = 1;  
vaxNet(12:14, 11) = 1; vaxNet(12:14, 12:14) = -eye(3);                   

end
