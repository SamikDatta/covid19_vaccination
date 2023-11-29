function writeTableLatex(results, fNameOut)

nScenarios = height(results);
tableOrder = [1 3 2 4 6 5 7 9 8];           % Order in which scenarios in 'results' will appear in formateted table
tableOrderSupp =   [1 2 6 5 7 9 3 4 9 8];   % Order in which scenarios in 'results' will appear in formateted supplementary table (Deltas)
baseForDeltaFlag = [1 0 0 0 0 0 1 0 1 0];   % Flag indicting which scnearios are comparators in the Supp Table against which Deltas are calculated for the subsequent scenarios
scenarioNumbers = ["0", "1", "0a", "1a", "3", "2", "4", "6", "5"]; 


scaleInf = 1e6;
scaleHosp = 1000;
scaleDeaths = 1;
scaleYLL = 1000;
scaleOcc = 1;
formInf = '%.2f';
formHosp = '%.1f';
formDeaths = '%.0f';
formYLL = '%.1f';
formOcc = '%.0f';

scenarioNames = string(results.scenario);
scenarioNames = strrep(scenarioNames, '%' , '\%');
scenarioNames = strrep(scenarioNames, 'Maori' , 'M\=aori');

fid = fopen(fNameOut, 'w');

nRows = length(tableOrder);
% Main table of results
for iRow = 1:nRows
    iScenario = tableOrder(iRow);
    infString = makeCI(results.nInfTot(iScenario, :), scaleInf, formInf  );
    hospString = makeCI(results.nHospTot(iScenario, :), scaleHosp, formHosp  );
    deathsString = makeCI(results.nDeathsTot(iScenario, :), scaleDeaths, formDeaths  );
    YLLString = makeCI(results.YLL(iScenario, :), scaleYLL, formYLL  );
    occString = makeCI(results.peakOcc(iScenario, :), scaleOcc, formOcc  );

    fprintf(fid, '(%s) %s & %s & %s & %s & %s & %s  \\\\ \n', scenarioNumbers(iScenario), scenarioNames(iScenario), infString, hospString, deathsString, YLLString, occString);
end

fprintf(fid, '\n\n');

% Supplementary table of Deltas

nRows = length(tableOrderSupp);
for iRow = 1:nRows
    iScenario = tableOrderSupp(iRow);
    if baseForDeltaFlag(iRow)
        fprintf(fid, '\\hline \n');
        fprintf(fid, '(%s) %s & - & - & - & -   \\\\ \n', scenarioNumbers(iScenario), scenarioNames(iScenario));
    else
        infString = makeCI(results.dInf(iScenario, :), scaleInf, formInf  );
        hospString = makeCI(results.dHosp(iScenario, :), scaleHosp, formHosp  );
        deathsString = makeCI(results.dDeaths(iScenario, :), scaleDeaths, formDeaths  );
        YLLString = makeCI(results.dYLL(iScenario, :), scaleYLL, formYLL  );
        fprintf(fid, '(%s) %s & %s & %s & %s & %s   \\\\ \n', scenarioNumbers(iScenario), scenarioNames(iScenario), infString, hospString, deathsString, YLLString);
    end
end



fclose(fid);

