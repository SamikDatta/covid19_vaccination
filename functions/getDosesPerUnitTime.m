function [nDose1, nDose2, nDose3, nDose4] = getDosesPerUnitTime(t, par);

ta = par.date0:par.date0+par.tEnd;

nDose1 = interp1(ta, par.nDoses1Smoothed, t)';
nDose2 = interp1(ta, par.nDoses2Smoothed, t)';
nDose3 = interp1(ta, par.nDoses3Smoothed, t)';
nDose4 = interp1(ta, par.nDoses4Smoothed, t)';

