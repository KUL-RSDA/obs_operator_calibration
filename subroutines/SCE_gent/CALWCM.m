function [obj]=CALWCM(X,mat_cell);

A=X(1);
B=X(2);
C=X(3);
D=X(4);
 
mat=cell2mat(mat_cell);
sm=mat(1,:)';
veg=mat(2,:)';
obs=mat(3,:)';

n_el=length(sm);

% Water Cloud Model
y0soil=C+D.*sm;
y0soil_lin=10.^(y0soil./10);
tt=exp(-2*B.*veg);
y0veg=(A.*veg);
y0tot_lin=(1-tt).*y0veg+tt.*y0soil_lin;
y0tot=10.*log10(y0tot_lin);
% Objective function

% Parameter deviation
PD_A=((0-A)^2)/(1.0^2/12);
PD_B=((0-B)^2)/(0.5^2/12);
PD_C=((-15-C)^2)/(25^2/12);
PD_D=((40-D)^2)/(65^2/12);
PD=(1/4)*(PD_A+PD_B+PD_C+PD_D);

% RMSE
RMSE=sqrt((1/n_el)*sum((y0tot-obs).^2));

obj=RMSE+0.1*PD;
%obj=RMSE+0*PD;
end
