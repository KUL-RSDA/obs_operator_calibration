%
% Set parameters for OF
%
% Input: choice of objective function
%
% Last update:
% 13 Dec 2019: Gabrielle De Lannoy
% 14 Mar 2023: Gabrielle De Lannoy, 
%                 clean up to only work for SSE for now
%==========================================================================

function [OF_param] = Initialize_OFparam_cali(OF_option,cali_option)
    
disp('Initializing OF parameters');

if (strcmp(OF_option,'SSE'))

    %W = Weights on each term in the OF:
    OF_param.W_p   = 1;             %Weight on parameter penalty

    %Prior is treated separately in DREAMZS and MCMC
    if (strcmp(cali_option,'DREAMZS') || strcmp(cali_option,'MCMC'))
      OF_param.W_p   = 0;           %Weight on parameter penalty 
    end

    OF_param.OF_type = 'SSE';

else

  error('this OF option is not ready for WCM calibration')

end


%=============================EOF==========================================


