function [myDataPath, dataFileNames] = getDataFileNames(dateLbl)


myDataPath = 'data/';

dataFileNames.epiDataFname =  "epidata_by_age_and_vax_" + dateLbl + ".mat";    % Line data
dataFileNames.vaxDataFname = "vaccine_data_national_2023-06-06";              % Vax data
dataFileNames.vaxProjFname = "reshaped_b2_projections_final_2022-07-13.csv";  % Vax projections (currently not used)
dataFileNames.AVdataFname = "therapeutics_by_age_14-Aug-2023.mat";            % Antiviral data
dataFileNames.hospOccFname = "covid-cases-in-hospital-counts-location-16-Aug-2023.xlsx";           % Only used for plotting
dataFileNames.popSizeFname = "popsize_national.xlsx";                         % NZ population structure (HSU)
% popSizeFname = "popproj_national2018-21.xlsx";                % NZ population structure (Stats NZ population projections)
dataFileNames.CMdataFname = "nzcontmatrix.xlsx";                              % Prem contact matrix
dataFileNames.borderIncFname = "border_incidence.xlsx";
