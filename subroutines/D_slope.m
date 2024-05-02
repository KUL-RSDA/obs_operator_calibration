function [param_0] = D_slope(param_0, DD, ...
    LIS_ssm, LIS_veg, LIS_obs, N_min)

% Use a change detection approach to constrain D, and subsequently C. 
% Michel Bechtold, Gabrielle De Lannoy - 17 Mar 2023
%=========================================================================

% Optional settings
%
% Narrow prior parameter range around estimated D [D-D_plusmin D+D_plusmin]
D_plusmin = 5;
% Quantile to split the year into low and high vegetation density
% 0.5 will split it into two parts of six months
% D will be estimated only for low vegetation period
vegetation_quantile = 0.5;

%--------------------------------------------------------------------------

% Constrain D with change detection approach only 
% for period when both S1 A and B are available 
% (i.e. repeat time of 6 days)

period_S1AB = DD>datetime(2016,10,1) & DD<datetime(2022,12,30);
ssm_S1AB = squeeze(LIS_ssm(period_S1AB));
veg_S1AB = squeeze(LIS_veg(period_S1AB));
obs_S1AB = squeeze(LIS_obs(period_S1AB));
% determine veg threshold for low vegetation based on defined vegetation_quantile
veg_threshold = quantile(veg_S1AB,vegetation_quantile);
% calculate 6 day delta
delta_ssm = ssm_S1AB(7:end)-ssm_S1AB(1:end-6);
delta_obs = obs_S1AB(7:end)-obs_S1AB(1:end-6);
cond1 = squeeze(~isnan(delta_obs));
cond2 = squeeze(veg_S1AB(4:end-3)<veg_threshold);
delta_ssm_nonan = delta_ssm(cond1&cond2);
delta_obs_nonan = delta_obs(cond1&cond2);
% constrain D if there are enough delta points (~50)
if length(find(~isnan(delta_obs_nonan)))>N_min
    pf = polyfit(delta_ssm_nonan, delta_obs_nonan, 1)';
    param_0.m.D=pf(1);%[db*m3/m3]
    % reset estimated D if outside range
    param_0.m.D = max([param_0.r.D(1) param_0.m.D]);
    param_0.m.D = min([param_0.r.D(2) param_0.m.D]);
    % set new narrowed range within defined narrower range
    param_0.r.D=[max([param_0.r.D(1) param_0.m.D-D_plusmin]) min([param_0.r.D(2) param_0.m.D+D_plusmin])]; %[dB*m3/m3]
    % calculate initial value for C based on estimated D and low vegetation season
    condveg = veg_S1AB<veg_threshold;
    param_0.m.C = nanmean(obs_S1AB(condveg)) - param_0.m.D * mean(ssm_S1AB(condveg));
end

    
end    
%---------------------------EOF--------------------------------------------
