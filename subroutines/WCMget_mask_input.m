function [fwd_in, obs, dates] = WCMget_mask_input(...
        DD, ssm, veg, obs, theta)

% Get LIS output and obs, and crossmask all
% Gabrielle De Lannoy - 13 March 2023
%==========================================================================
 
M = [DD' ssm veg obs];
% remove all lines/times with any missing data
M(isnan(mean(M,2)),:)=[];

fwd_in.sm   = M(:,2);
fwd_in.veg  = M(:,3);
fwd_in.theta= theta+zeros(length(M(:,1)),1);
obs         = M(:,4);
dates       = M(:,1);

clear M     
            
end
%-------------------------------EOF----------------------------------------
