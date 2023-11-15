function [dates, doses1, doses2, doses3, doses4plus] = ...
    getVaccineData(myDataPath, vaxDataFname, vaxProjFname, vaccImmDelay, date0, tEnd)
% Function that reads in vaccination data and projections (if used), and
% outputs arrays of cumulative doses for each age group, shifted by a 
% delay to account for the time it takes for the vaccine to increase 
% immunity, and cropped to match simulation's time array. 
% If projection not used, a flat tail will be added at the end of
% available data.
% INPUTS:
% - myDataPath: path to folder where vaccination data is stored
% - vaxDataFname: filename of vaccination data
% - vaxProjFname: filename of vaccination projections data
% - vaccImmDelay: number of days after which vax dose is assumed to take
%                 effect (as defined in getBasePar.m)
% - date0: datenum corresponding to simulation start date (as defined in
%          getBasePar.m)
% - tEnd: number of simulated days (as defined in getBasePar.m)
% OUTPUTS:
% - dates: array of datenum dates corresponding to delayed vax doses
% - doses1-4plus: arrays of cumulative doses for each shifted date

fName = myDataPath + vaxDataFname;
vaxData = readtable(fName);

% Set use_projections to 1 if wanting to use projected vaccination data, if
% set to 0 a flat tail is added after the end of the available vax data
use_projections = 0;

if use_projections == 1
    fNameProj = myDataPath + vaxProjFname; 
    vaxProj = readtable(fNameProj);

    % Join actuals and projections into one table
    vaxMerged = outerjoin(vaxData, vaxProj, 'LeftKeys', 'dates', 'RightKeys', 'date', 'MergeKeys', true);

    % Replacing 4th dose for 50+ with projections once data has run out
    fourth_dose_data = 60:65; % selecting columns for dose 4 50-54, ..., 75+
    disp(vaxMerged.Properties.VariableNames(fourth_dose_data)); % check names
    fourth_dose_proj = 66:71; % columns for fourth dose projections
    disp(vaxMerged.Properties.VariableNames(fourth_dose_proj)); % check names
    smoothing_period = 7; % number of days to smooth between data and projections
    for j = 1:length(fourth_dose_data)
        ind = find(isnan(vaxMerged{:, fourth_dose_data(j)}), 1, 'first'); % get first date without data
        if vaxMerged{ind-1, fourth_dose_data(j)} > vaxMerged{ind-1, fourth_dose_proj(j)} % if data lies above projection, smooth to projection
            vaxMerged(ind+smoothing_period:end, fourth_dose_data(j)) = vaxMerged(ind+smoothing_period:end, fourth_dose_proj(j)); % replace after smoothing period
            vaxMerged{ind:ind+smoothing_period-1, fourth_dose_data(j)} = ...
                round(interp1(datenum(vaxMerged.dates_date), vaxMerged{:, fourth_dose_data(j)}, ...
                datenum(vaxMerged.dates_date(ind:ind+smoothing_period-1)), 'spline')); % smooth for smoothing period
        else % if data below projection, pull projection down to match last data point
            scalar = vaxMerged{ind-1, fourth_dose_data(j)}/vaxMerged{ind-1, fourth_dose_proj(j)}; % scalar
            vaxMerged{:, fourth_dose_proj(j)} = round(scalar*vaxMerged{:, fourth_dose_proj(j)}); % scale down whole curve
            vaxMerged(ind:end, fourth_dose_data(j)) = vaxMerged(ind:end, fourth_dose_proj(j)); % replace after data ends
        end
    end

    % Overwrite NaNs by filling the last non-nan value down the column (skipping 1st column which is dates)
    [~, nCols] = size(vaxMerged);
    for iCol = 2:nCols
        ind = find(~isnan(vaxMerged{:, iCol}), 1, 'last');
        vaxMerged(ind+1:end, iCol) = vaxMerged(ind, iCol);
    end

    dates = datenum(vaxMerged.dates_date)';
    doses1 =     [vaxMerged.doses1_1 vaxMerged.doses1_2  vaxMerged.doses1_3  vaxMerged.doses1_4  vaxMerged.doses1_5  vaxMerged.doses1_6  vaxMerged.doses1_7  vaxMerged.doses1_8  vaxMerged.doses1_9  vaxMerged.doses1_10  vaxMerged.doses1_11  vaxMerged.doses1_12  vaxMerged.doses1_13  vaxMerged.doses1_14  vaxMerged.doses1_15  vaxMerged.doses1_16 ];
    doses2 =     [vaxMerged.doses2_1 vaxMerged.doses2_2  vaxMerged.doses2_3  vaxMerged.doses2_4  vaxMerged.doses2_5  vaxMerged.doses2_6  vaxMerged.doses2_7  vaxMerged.doses2_8  vaxMerged.doses2_9  vaxMerged.doses2_10  vaxMerged.doses2_11  vaxMerged.doses2_12  vaxMerged.doses2_13  vaxMerged.doses2_14  vaxMerged.doses2_15  vaxMerged.doses2_16 ];
    doses3 =     [vaxMerged.doses3_1 vaxMerged.doses3_2  vaxMerged.doses3_3  vaxMerged.doses3_4  vaxMerged.doses3_5  vaxMerged.doses3_6  vaxMerged.doses3_7  vaxMerged.doses3_8  vaxMerged.doses3_9  vaxMerged.doses3_10  vaxMerged.doses3_11  vaxMerged.doses3_12  vaxMerged.doses3_13  vaxMerged.doses3_14  vaxMerged.doses3_15  vaxMerged.doses3_16 ];
    doses4plus = [vaxMerged.doses4_1 vaxMerged.doses4_2  vaxMerged.doses4_3  vaxMerged.doses4_4  vaxMerged.doses4_5  vaxMerged.doses4_6  vaxMerged.doses4_7  vaxMerged.doses4_8  vaxMerged.doses4_9  vaxMerged.doses4_10  vaxMerged.doses4_11  vaxMerged.doses4_12  vaxMerged.doses4_13  vaxMerged.doses4_14  vaxMerged.doses4_15  vaxMerged.doses4_16 ];

else
    dates = datenum(vaxData.dates)';
    doses1 = table2array(vaxData(:, 2:17));
    doses2 = table2array(vaxData(:, 18:33));
    doses3 = table2array(vaxData(:, 34:49));
    doses4plus = table2array(vaxData(:, 50:65));

end


% %% DIAGNOSTIC PLOTS to run when changing vax data file
% Note: the 4+ cumulative doses will cross the 1-3 curves since Apr23
% because of the rollout of the new booster. To check 4th and 5th doses
% separately, use utils4_5doses.m and check_4th_5th_doses.m in the
% processVaxData folder
plotVax = 0;
if plotVax == 1
    f = figure;
    f.Position = [100 100 800 400];
    titles = {'0-4', '5-9', '10-14', '15-19', '20-24', '25-29', ...
        '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', ...
        '65-69', '70-74', '75+'};
    dtimedates = datetime(dates,'ConvertFrom','datenum');
    f.Position = [100 100 600 400];
    tiledlayout(2, 2);
    for ag = 13:16
        nexttile
        hold on
        title(titles(ag))
        plot(dtimedates, doses1(:, ag), dtimedates, doses2(:, ag), ...
            dtimedates, doses3(:, ag), dtimedates, doses4plus(:, ag))
        hold off
        if ag == 14
            l = legend({'1+ dose', '2+ doses', '3+ doses', '4+ doses'});
            l.Layout.Tile = 'East';
        end
    end
end



% Shift vaccination dates to allow for delay in taking effect
dates = dates + vaccImmDelay;     

% Find index of vax data line that corresponds to first day of simulation
iDate = find(datenum(dates) == date0);

% Find out how many more (or less) dates there are in the vax data compared
% to the simulations's time array
nPad = tEnd - (length(dates) - iDate);
if nPad > 0
    % If the simulation has more dates than the vax data, add a flat tail
    % at the end of the vax data
    dates = [dates, dates(end)+1:dates(end)+nPad];
    doses1 = [doses1(1:end, :); repmat(doses1(end, :), nPad, 1)];
    doses2 = [doses2(1:end, :); repmat(doses2(end, :), nPad, 1)];
    doses3 = [doses3(1:end, :); repmat(doses3(end, :), nPad, 1)];
    doses4plus = [doses4plus(1:end, :); repmat(doses4plus(end, :), nPad, 1)];
elseif nPad < 0
    % If the simulation has less dates than the vax data, cut vaccination
    % data down to match
    dates = dates(1:end+nPad);
    doses1 = doses1(1:end+nPad, :);
    doses2 = doses2(1:end+nPad, :);
    doses3 = doses3(1:end+nPad, :);
    doses4plus = doses4plus(1:end+nPad, :);
end

% Crop out excess dates from beginning of vax data, if any
dates = dates(iDate:end);
doses1 = doses1(iDate:end, :);
doses2 = doses2(iDate:end, :);
doses3 = doses3(iDate:end, :);
doses4plus = doses4plus(iDate:end, :);




end


