function [OF, OF_sse, OF_p] = CalWCM_SSE(fwd_in,param_0_tmp,obs,...
            p_constraint,OF_param,select_par)

% Gabrielle De Lannoy - Bayesian OF approximation
% Gabrielle De Lannoy - 13 March 2023
%==========================================================================
            
fn = fieldnames(param_0_tmp);

A=param_0_tmp.A;
B=param_0_tmp.B;
C=param_0_tmp.C;
D=param_0_tmp.D;

sm_arr   =fwd_in.sm;
veg_arr  =fwd_in.veg;
theta_arr=fwd_in.theta;
obs_arr  =obs;

si=param_0_tmp.s_0; %this is the obs sigma which goes into the SSE

% Water Cloud Model
sigma_tot = WCM(A,B,C,D,sm_arr,veg_arr,theta_arr); %linear scale
sigma_tot = 10.*log10(sigma_tot); %sigma_tot converted in dB

%sum of square errors for each pixel
SSE = ((obs_arr-sigma_tot).^2)./(2*si^2);
OF_sse=sum(SSE); 

clear sigma_tot obs_arr veg_arr sm_arr 

% OF = SSE (Sum of Square Errors) + Sum of squared parameter deviations
% limited by the variance of a uniform prior parameter distribution

par_size  = length(select_par);
OF_p  = 0;

if OF_param.W_p ~= 0

  for p=1:par_size

    OF_p = OF_p + ...
        (p_constraint.pen.(fn{select_par(p)})...
         -param_0_tmp.(fn{select_par(p)}))^2/...
        (2*p_constraint.var.(fn{select_par(p)})); 
    % second term is the variance of uniform distribution

  end

  OF_p=OF_param.W_p*(OF_p); %parameter penalty term

end

OF=OF_sse+OF_p; %sum of square error weighted by a parameter penalty term

end

%================================EOF=======================================
