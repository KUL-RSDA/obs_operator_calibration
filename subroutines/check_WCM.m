%Quick checking of WCM input, output, parameters
%GDL, June 2021
%------------------------------------------------------------------------

function check_WCM(dates,fwd_in,obs,param_0,param_opt,cali_option,f_out)

%Run WCM with initial parameters
par = param_0.m;
sigma_0 = WCM(par.A,par.B,par.C,par.D,fwd_in.sm,fwd_in.veg,fwd_in.theta);
sigma_0 = 10.*log10(sigma_0);
%disp(sigma_0);

%Run WCM with optimized parameters
par = param_opt;
sigma_opt = WCM(par.A,par.B,par.C,par.D,fwd_in.sm,fwd_in.veg,fwd_in.theta);
sigma_opt = 10.*log10(sigma_opt);
%disp(sigma_opt);

%Plotting
figure('units','centimeters','position',[10 10 40 15]); hold on
subplot(3,1,1)
yyaxis left
plot(dates,fwd_in.sm,'Color',[0.3,0.3,1],'LineWidth',2); hold on;
ylabel('Soil moisture [m³/m³]');
yyaxis right
plot(dates,obs,'ro','Markersize',2,'MarkerFaceColor','r'); hold on;
plot(dates,sigma_0,'Color',[0.25, 0.25, 0.25]); hold on;
plot(dates,sigma_opt,'k.','LineWidth',2); hold on;
ylabel('\gamma_0 [dB]');
axis tight
ax = gca;
ax.YAxis(1).Color = [0.3 0.3 1];
ax.YAxis(2).Color = [1 0 0];
legend('LIS SM','obs','\gamma\_ini','\gamma\_opt','Location','South','Orientation','horiz');

subplot(3,1,2);
yyaxis left
plot(dates,fwd_in.veg,'Color',[0 0.8 0],'LineWidth',2);hold on;
ylabel('Vegetation, LAI [-]');
yyaxis right
plot(dates,obs,'ro','Markersize',2,'MarkerFaceColor','r'); hold on;
plot(dates,sigma_opt,'k.','LineWidth',2); hold on;
ylabel('\gamma_0 [dB]');
axis tight
ax = gca;
ax.YAxis(1).Color = [0 0.8 0];
ax.YAxis(2).Color = [1 0 0];
legend('LIS Veg','obs','\gamma\_opt','Location','South','Orientation','horiz');

subplot(3,2,5);
scatter(obs,sigma_0,2,repmat([0.6 0.6 0.6]',1,length(obs))','filled');hold on;
scatter(obs,sigma_opt,2,repmat([0 0 0]',1,length(obs))','filled');hold on;
xlabel('obs \gamma_0 [dB]');ylabel('sim [dB]');
ax=gca;
mi=min(min(ax.XLim),min(ax.YLim));
ma=max(max(ax.XLim),max(ax.YLim));
plot([mi ma],[mi ma]);
legend('\gamma\_ini','\gamma\_opt','Location','NorthWest');
title(cali_option);

skill_text = {['R_0=',num2str(corr(obs,sigma_0),2),'; R\_opt=',num2str(corr(obs,sigma_opt),2),' [-]'],...
    ['rmsd_0=',num2str(sqrt(mean((obs-sigma_0).^2)),3),'; rmsd\_opt=',num2str(sqrt(mean((obs-sigma_opt).^2)),3),' [dB]']};
if exist('par.s_0')>0
  param_text = ['A=',num2str(par.A,2),'; B=',num2str(par.B,2),'; C=',num2str(par.C,2),'; D=',num2str(par.D,2),'; s_0=',num2str(par.s_0,2)];
else
  param_text = ['A=',num2str(par.A,2),'; B=',num2str(par.B,2),'; C=',num2str(par.C,2),'; D=',num2str(par.D,2)];    
end

annotation('textbox',[0.48 0.25 0.46 0.05],'String',param_text)
annotation('textbox',[0.48 0.15 0.46 0.1],'String',skill_text)

end
%=========EOF===============================
