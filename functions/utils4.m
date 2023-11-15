
%% Process and combine previous and bookings

clear all
close all


regionDHBs{1} = categorical(["Northland", "Waitemata", "Auckland", "Counties Manukau"]);
regionDHBs{2} = categorical(["Waikato", "Bay of Plenty", "Lakes", "Tairawhiti", "Taranaki"]);
regionDHBs{3} = categorical(["Whanganui", "Hawkes Bay", "Wairarapa", "Hutt Valley", "Capital and Coast", "MidCentral"]);
regionDHBs{4} = categorical(["Nelson Marlborough", "West Coast", "Canterbury", "South Canterbury", "Southern"]);

datelbl_administered = "2023-03-06";
%datelbl_forecast = "2022-01-26";
nAgeGroups = 16;

%F0 = readtable(sprintf("vaccine_bookings_%s_withage.xlsx", datelbl_forecast));
T0 = readtable(sprintf("TPM_vaccine_%s.csv", datelbl_administered));
%%
T0.DHB_OF_RESIDENCE = categorical(T0.DHBOFRESIDENCE);
T0.ageBand = min(16, 1+floor(double(string(T0.AGEATVACCINATIONEVENT))/5) );
validAgeFlag = (T0.ageBand >= 1 & T0.ageBand <= 16);
fprintf('Discarding elements of T0 with invalid age:\n')
T0(~validAgeFlag, :)

T0 = T0(validAgeFlag, :);
% Make sure T0 has an entry for each age band (with count = 0)
nRows = height(T0);
T0 = [T0; T0(1:nAgeGroups, :)];
T0.count = [ones(nRows, 1); zeros(nAgeGroups, 1)];
T0.ageBand(nRows+1:nRows+nAgeGroups) = 1:nAgeGroups;
Uc = unstack(T0, 'count', 'ageBand', 'AggregationFunction', @sum, 'GroupingVariables', {'ACTIVITYDATE', 'DOSENUMBER', 'DHB_OF_RESIDENCE'});

% F0.ageBand = min(16, 1+double(extractBefore(string(F0.AGE_BAND_5), '-'))/5);
% F0.ORGANISATION_NAME = categorical(extractBefore(string(F0.ORGANISATION_NAME), " District Health Board"));
% F0.DOSE_NUMBER = double(string(F0.DOSE_NUMBER));
% F0 = table(F0.APPOINTMENT_DATE, F0.ORGANISATION_NAME, F0.DOSE_NUMBER, F0.count, F0.ageBand, 'VariableNames', {'ACTIVITYDATE', 'DHB_OF_RESIDENCE', 'DOSENUMBER', 'count', 'ageBand'});
% nRows = height(F0);
% F0 = [F0; F0(1:nAgeGroups, :)];
% F0.count(nRows+1:nRows+nAgeGroups) = 0;
% F0.ageBand(nRows+1:nRows+nAgeGroups) = 1:nAgeGroups;
% W0 = unstack(F0, 'count', 'ageBand', 'AggregationFunction', @sum, 'GroupingVariables', {'ACTIVITYDATE', 'DOSENUMBER', 'DHB_OF_RESIDENCE'});

%Uc = [U0; W0(W0.ACTIVITYDATE > max(U0.ACTIVITYDATE), :)];


%%


st = min(Uc.ACTIVITYDATE); en = max(Uc.ACTIVITYDATE); 
dts = st:en;
% writematrix(dts, "processed_for_modelling/vaccine_dates.csv");

change1 = ["national", "northern" , "midland", "central", "southern", string(regionDHBs{1}), string(regionDHBs{2}), string(regionDHBs{3}), string(regionDHBs{4})];


for iChange1 = 1:length(change1)

    REGION = change1(iChange1)

    % Load Data

    if REGION == "national"
        inRegionFlag = ones(height(Uc), 1);
    elseif REGION == "northern"
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, regionDHBs{1});
    elseif REGION == "midland"
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, regionDHBs{2});
    elseif REGION == "central"
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, regionDHBs{3});
    elseif REGION == "southern"
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, regionDHBs{4});
    elseif REGION == "auckland_metro"
        DHBs = categorical(["Waitemata", "Auckland", "Counties Manukau"]);
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, DHBs);
    elseif ismember(REGION, Uc.DHB_OF_RESIDENCE)
        inRegionFlag = ismember(Uc.DHB_OF_RESIDENCE, REGION);
    else
        error("Please choose a valid REGION");
    end

    dose1 = zeros(nAgeGroups, length(dts));
    dose2 = zeros(nAgeGroups, length(dts));
    dose3 = zeros(nAgeGroups, length(dts));
    dose4plus = zeros(nAgeGroups, length(dts));
    dose1(:, 1) = nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(1) & Uc.DOSENUMBER == 1, 4:end)), 1 )';
    dose2(:, 1) = nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(1) & Uc.DOSENUMBER == 2, 4:end)), 1 )';
    dose3(:, 1) = nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(1) & Uc.DOSENUMBER == 3, 4:end)), 1 )';
    dose4plus(:, 1)     = nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(1) & Uc.DOSENUMBER >= 4, 4:end)), 1 )';
    for ii = 2:length(dts)
        dt = dts(ii);
        dose1(:, ii) = dose1(:, ii-1) + nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(ii) & Uc.DOSENUMBER == 1, 4:end)), 1 )';
        dose2(:, ii) = dose2(:, ii-1) + nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(ii) & Uc.DOSENUMBER == 2, 4:end)), 1 )';
        dose3(:, ii) = dose3(:, ii-1) + nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(ii) & Uc.DOSENUMBER == 3, 4:end)), 1 )';
        dose4plus(:, ii) = dose4plus(:, ii-1)         + nansum(table2array(Uc(inRegionFlag & Uc.ACTIVITYDATE == dts(ii) & Uc.DOSENUMBER >= 4, 4:end)), 1 )';
    end

    V1 = max(dose1-dose2, 0);
    V2 = max(dose2-dose3, 0);
    V3 = max(dose3-dose4plus, 0);
    V4 = dose4plus;
    
     figure
     subplot(2, 2, 1)
    plot(dts, dose1); title("Number of Each Age Group With At Least One Dose"); xline(datetime("today"));
    subplot(2, 2, 2)
    plot(dts, dose2); title("Number of Each Age Group With At Least Two Dose"); xline(datetime("today"));   
    subplot(2, 2, 3)
    plot(dts, dose3); title("Number of Each Age Group With At Least Three Dose"); xline(datetime("today"));
    subplot(2, 2, 4)
    plot(dts, dose4plus); title("Number of Each Age Group With At Least Four Dose"); xline(datetime("today"));
    
%      writematrix(V1, sprintf("processed_for_modelling/firstDoseCount_%s.csv", REGION ));
%      writematrix(V2, sprintf("processed_for_modelling/secondDoseCount_%s.csv", REGION ));
%      writematrix(V3, sprintf("processed_for_modelling/thirdDoseCount_%s.csv", REGION ));
%      writematrix(V4, sprintf("processed_for_modelling/fourthDoseCount_%s.csv", REGION ));
     
     outTab = [];
     outTab.dates = dts';    
     outTab.doses1 = dose1';
     outTab.doses2 = dose2';
     outTab.doses3 = dose3';
     outTab.doses4 = dose4plus';
     outTab = struct2table(outTab);
     
     
     writetable(outTab, sprintf("processed_for_modelling/vaccine_data_%s_%s.csv", REGION, datelbl_administered ));
     
end



