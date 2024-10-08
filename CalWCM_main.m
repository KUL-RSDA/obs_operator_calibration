%
% A. Initialize: set calibration options
% B. Datasets: get input and observation data from LIS output
% C. Calibrate
% D. Write parameter output
%
% Gabrielle De Lannoy - 13 March 2023: Major clean up for 4DMED
%
%% ==============A. INITIALIZE=============================================

clc;
clear all;

debug = 1;  %0/1 to check and get some display

path_subr = ['/data/leuven/314/vsc31402/obs_operator_calibration/',...
             'subroutines/'];
addpath(path_subr); 

%--------------------------------------------------------------------------                       

% Calibration algorithms   
cali_option  = 'DREAMZS'; %'DREAMZS'; %'SCE_gent'; %'pso'; 'sceua'
                          %'sceua' WARNING: with D via change detection 
                          %    does not converge well
% Objective function to optimize  
OF_option    = 'SSE';           % 'SSE' (sum(backscatter dev) + sum(prior dev))
                                % other options not ready
% Underlying LSM or crop model
model        = 'LSM';           % 'LSM'; 'AquaCrop';

% Calibration time period                        
start_cal_date = '15-apr-2015'; %'28-dec-2021'; %'01-jan-2018'
end_cal_date   = '15-apr-2023'; %'30-dec-2021'; % 31<-> 30-dec // 31 dec not available 
                                % LIS SURFACEMODEL saves previous day
exclude_months = [];            % e.g. [11,12]; none: []
N_min          = 50;            % minimum number of data pairs 
                                % needed for calibration; 
                                % need ~50 if D is a priori constrained
                                % before calibration;
                                % insufficient data -> set param to NaN

% Initial parameter & obs error values 
% depend on polarization and incidence angle
pol         = 'VV';             %or 'VH';        
theta       = 0;                %or 37, for now assumed constant
r_A_B_to_37 = 1;                %rescale UB of A and B to theta=37 
                                %==> 1 for Sara; 0 for Shannon

                        
% Input/Output specifications
experiment      = 'po_ol_hymap_irr'; %'po_ol_hymap_noirr';
path_LIS_output = '/staging/leuven/stg_00024/OUTPUT/projects/4DMED/output_2015-2023/';
path_cal_output = '/scratch/leuven/314/vsc31402/4DMED/';

% experiment      = 'po_ol_no_irr'; %'po_ol_irr';
% path_LIS_output = '/staging/leuven/stg_00024/OUTPUT/projects/4DMED/';
% path_cal_output = '/scratch/leuven/314/vsc31402/4DMED/';

% experiment      = 'GERMANY'; 
% path_LIS_output = '/staging/leuven/stg_00024/OUTPUT/louiseb/WCM_calib/';
% path_cal_output = '/scratch/leuven/314/vsc31402/Louise/';
 
% Domain (e.g. subdomain)
min_lat  = NaN; %44.2350; %NaN;
max_lat  = NaN; %44.4850; %NaN;
min_lon  = NaN; %11.5050; %NaN;
max_lon  = NaN; %11.9750; %NaN;

% Select which parameters to calibrate
% 'A', 'B', 'C', 'D', 's_0' = param_0, = all parameters (below)
par_names = {'A', 'B', 'C', 'D', 's_0'};

% Use change detection approach to constrain D 
D_change_detect = 1;                       

%--------------------------------------------------------------------------                       

path_LIS_output = [path_LIS_output,'/',experiment,'/'];
path_cal_output = [path_cal_output,'/',experiment,'/'];

% Initialize all parameters
[param_0] = WCMparam_ini(pol, theta, r_A_B_to_37, model); 
fn        = fieldnames(param_0.m);
ftag      = '';
problem   = D_change_detect;
for p=1:length(par_names)
    select_par(p) = find(strcmp(par_names{p},fn));
    ftag = [ftag,fn{p}];
    if strcmp(par_names{p},'D') && D_change_detect == 1
        problem = 0;
    end
end
ftag = [ftag,num2str(D_change_detect)];
if problem
    error('A priori slope detection is only useful when D is calibrated');
end

% Set calibration options
par_opt  = Initialize_Algparam_cali(cali_option,select_par,path_subr);
OF_param = Initialize_OFparam_cali(OF_option,cali_option);

%% ============B. DATA SPECS===============================================

% Load data from LIS output, where simulations and observations are aligned
% and perfectly masked (e.g. flags for soil temperature, snow, etc.).
% Limit data to calibration period and domain only.

path_cal_output = [path_cal_output,'/',start_cal_date,'_',end_cal_date,'/']; 
if ~exist(path_cal_output, 'dir')
     mkdir(path_cal_output);
end

% Optionally make a domain subdirectory for the output
if (~isnan(min_lon) && ~isnan(min_lat) && ...
    ~isnan(min_lat) && ~isnan(max_lat))
    subdir = [num2str(round(100*min_lon),'%d'),'_',...
              num2str(round(100*max_lon),'%d'),'_',...
              num2str(round(100*min_lat),'%d'),'_',...
              num2str(round(100*max_lat),'%d')];
    path_cal_output = [path_cal_output,'/',subdir,'/']; 
    if ~exist(path_cal_output, 'dir')
        mkdir(path_cal_output);
    end
end

file_out_mat = [path_cal_output,'/',experiment,'.mat'];

if ~exist(file_out_mat,'file')
    
    % Script is hardwired for daily LIS output
    [DD, lat_all, lon_all, lat, lon, domain, LIS_ssm, LIS_veg, LIS_obs]...
        = read_LIS_output(...
        path_LIS_output, start_cal_date, end_cal_date, exclude_months,...
        min_lat, max_lat, min_lon, max_lon, '.a01', 1, model);
    save (file_out_mat,...
        'DD', 'lat_all', 'lon_all', 'lat', 'lon', 'domain',...
        'LIS_ssm', 'LIS_veg', 'LIS_obs' ,'-v7.3');
    
else
    disp(['Loading matlabfile ',file_out_mat]);
    load(file_out_mat);
end

%--------------------------------------------------------------------------                       

% Initialize calibrated parameters
par_size= length(select_par);
for p=1:par_size
    eval([fn{select_par(p)},'_cal = NaN(length(domain),1);']);
    eval([fn{select_par(p)},'_std = NaN(length(domain),1);']);
end 

%% ============C. CALIBRATION LOOP=========================================
% % test for some individual locations
% l_lon = [11,    12];
% l_lat = [44.75, 44.75];
% for t=1:length(l_lon)
%     dist = abs(l_lat(t)-lat_all) + abs(l_lon(t)-lon_all);
%     ll(t) = find(dist == min(dist));
% end
%
% % test for a small domain
% ll = find(lat_all >= 45.003 & lat_all <= 45.004 & lon_all >= 10.9895 & lon_all <= 10.9899);
% ll = find(lat_all >= 44.697 & lat_all <= 44.698 & lon_all >= 9.0006 & lon_all <= 9.0008);
% 
% %for  i= ll
for  i= 1:length(domain) 

    % progr = i/(end_i-start_i) * 100;
    % disp(sprintf('%0.2f / 100',progr))
    
    % Possibly update the prior param_0 per pixel i, if the D parameter
    % is constrained based on the available obs.
    param_0_i = param_0;
    if D_change_detect
        [param_0_i] = D_slope(param_0_i, DD, squeeze(LIS_ssm(:,i)),...
            squeeze(LIS_veg(:,i)), squeeze(LIS_obs(:,i)), N_min);
    end
        
    [fwd_in, obs, dates] = WCMget_mask_input(...
        datenum(DD), LIS_ssm(:,i), LIS_veg(:,i), LIS_obs(:,i), theta);

    % Calibration for this location i
    if (isempty(fwd_in.sm) || isempty(fwd_in.veg) || isempty(obs) ||...
            length(obs) < N_min)
            
       disp(['NOT calibrating lat-lon: ', num2str(lat(i)),'-', num2str(lon(i))]);
       
    else
       disp(['Calibrating lat-lon: ', num2str(lat(i)),'-', num2str(lon(i))]);
    
       [param_opt, param_std] = WCMcalibration(i, cali_option, par_opt,...
            param_0_i, select_par, ...
            fwd_in, ...
            obs, OF_param,debug);
       
       for p=1:par_size
            eval([fn{select_par(p)},'_cal(i) = param_opt.',...
                fn{select_par(p)},';']);
            eval([fn{select_par(p)},'_std(i) = param_std.',...
                fn{select_par(p)},';']);
       end 
        
       if debug 
            disp(['Calibrated --',...
                  ' A ',num2str(A_cal(i)),' B ',num2str(B_cal(i)),...
                  ' C ',num2str(C_cal(i)),' D ',num2str(D_cal(i))]);
            f_out = [path_cal_output,'/check_',cali_option,pol,num2str(i)];     
            check_WCM(datetime(dates,'ConvertFrom','datenum'),...
                fwd_in,obs,param_0_i,param_opt,cali_option,f_out);
       end
        
    end
    
end    

%% ============D. WRITE OUTPUT=============================================

filename = [path_cal_output,'Par_s',pol,'_',cali_option,'_',ftag,'_'];

for p=1:par_size
    fname = [filename,'_',fn{select_par(p)},'.txt'];
    save (fname, [fn{select_par(p)},'_cal'], '-ascii')
    fname = [filename,'_',fn{select_par(p)},'_sd.txt'];
    save (fname, [fn{select_par(p)},'_std'], '-ascii')
end 
save ([filename,'_lat.txt'], 'lat', '-ascii')
save ([filename,'_lon.txt'], 'lon', '-ascii')

%if debug
for p=1:par_size
    eval([fn{select_par(p)},'_cal_all = NaN(length(lon_all),1);']);
    eval([fn{select_par(p)},'_std_all = NaN(length(lon_all),1);']);
    eval([fn{select_par(p)},'_cal_all(domain)=',fn{select_par(p)},'_cal;']);
    eval([fn{select_par(p)},'_std_all(domain)=',fn{select_par(p)},'_cal;']);
end 
save ([filename,'.mat'],'lat_all', 'lon_all',...
      'A_cal_all','B_cal_all','C_cal_all','D_cal_all',...
      'A_std_all','B_std_all','C_std_all','D_std_all','-v7.3');
%end

%=============================EOF=========================================
