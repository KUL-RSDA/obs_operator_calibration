function [param_0] = WCMparam_ini(pol, theta, r_A_B_to_37, model) 

% ALL initial parameter & obs error values (incl. those not calibrated)
% depend on polarization, incidence angle, and model SSM
% Gabrielle De Lannoy - 13 March 2023
%==========================================================================

% Initial values of parameters A, B, C, D, and obs error
% e.g. Modanesi et al. (2022), De Lannoy et al. (JAMES, 2024)
param_0.m.A = 0.;       %[linear scale]
param_0.m.B = 0.;       %[linear scale]
param_0.m.C = -20.;     %[dB]
param_0.m.D = 40.;      %[db*m3/m3]
param_0.m.s_0 = 1.;     %[dB] backscatter error std in obs space

if strcmp(pol,'VH')
    param_0.m.C = -30.; % [dB] (values of VH backscatter are lower)
end

% Range of values of parameters A, B, C, D, and obs error
% e.g. Modanesi et al. (2022), De Lannoy et al. (JAMES, 2024)
param_0.r.A = [0. 1.0];   %[linear scale]
param_0.r.B = [0. 0.5];   %[linear scale]
param_0.r.C = [-35 -10];  %[dB]
param_0.r.D = [15 80];    %[dB*m3/m3]
param_0.r.s_0 = [0.1 10.]; %prior range of backscatter obs error

if (theta ~= 37 && r_A_B_to_37)
    param_0.r.A = [0. 0.4/cosd(37)]; %[linear scale]
    param_0.r.B = [0. 0.4*cosd(37)]; %[linear scale] 
end

if strcmp(model,'AquaCrop')    % different range in SM
    % de Roos et al. (RSE, 2023)	
    param_0.r.C = [-35. -10.]; %[dB]
    param_0.r.D = [10. 80.];   %[dB*m3/m3]
    % de Roos et al. (JGR, Biogeosciences, 2024)
    %param_0.r.C = [-35. -5.]; %[dB]
    %param_0.r.D = [5. 80.];   %[dB*m3/m3]
end

end
%----------------------------------EOF-------------------------------------
