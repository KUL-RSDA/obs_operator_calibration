function [sigma_tot] = WCM(A,B,C,D,sm_arr,veg_arr,theta_arr)

% WCM returning backscatter (sigma, gamma,...) in 
% linear scale

sigma_s   =C+D.*sm_arr;       %sigma_bare calculated in dB
sigma_bare=10.^(sigma_s./10); %sigma bare goes into the WCM in linear scale

co=cos(theta_arr.*pi/180);
tt=exp(-2*B.*veg_arr./co);
sigma_can=(A.*veg_arr).*co.*(1-tt);
sigma_soil=tt.*sigma_bare;
sigma_tot=sigma_can+sigma_soil;

end


