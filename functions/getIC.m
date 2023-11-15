function IC = getIC(par)
% Get initial conditions for ODE

N0 = par.popCount;
V0 = [par.doses1(1, :)', par.doses2(1, :)', par.doses3(1, :)'];
S0 = [par.popCount-par.doses1(1, :)', par.doses1(1, :)'-par.doses2(1, :)', ...
    par.doses2(1, :)'-par.doses3(1, :)', zeros(par.nAgeGroups, 3), ...
    par.doses3(1, :)', zeros(par.nAgeGroups, 7)];
E0 = zeros(size(S0));
I0 = zeros(size(S0));
A0 = zeros(size(S0));
R0 = zeros(par.nAgeGroups, par.nSusComp-1);
C00 = zeros(par.nAgeGroups, 3);
C10 = zeros(par.nAgeGroups, 3);
C20 = zeros(par.nAgeGroups, 3);
C30 = zeros(par.nAgeGroups, 3);
Cr0 = zeros(par.nAgeGroups, 3);
H00 = zeros(par.nAgeGroups, 5);
H10 = zeros(par.nAgeGroups, 5);
H20 = zeros(par.nAgeGroups, 5);
H30 = zeros(par.nAgeGroups, 5);
Hr0 = zeros(par.nAgeGroups, 5);
F00 = zeros(par.nAgeGroups, 6);
F10 = zeros(par.nAgeGroups, 6);
F20 = zeros(par.nAgeGroups, 6);
F30 = zeros(par.nAgeGroups, 6);
Fr0 = zeros(par.nAgeGroups, 6);
IC = [ N0(:); V0(:); S0(:); E0(:); I0(:); A0(:); R0(:); C00(:); C10(:); ...
    C20(:); C30(:); Cr0(:); H00(:); H10(:); H20(:); H30(:); Hr0(:); ...
    F00(:); F10(:); F20(:); F30(:); Fr0(:)];

end