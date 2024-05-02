function [dates, lat_all, lon_all, lat, lon, domain, LIS_ssm_sd, LIS_veg_swe, LIS_obs] = ...
    read_LIS_output(path, start_cal_date, end_cal_date, exclude_months,...
    min_lat, max_lat, min_lon, max_lon, DA_tag, read_obs, pol)

% Sara Modanesi & Gabrielle De Lannoy - 13 March 2023 
% Gabrielle De Lannoy - 15 March 2023:
%    - Major overhaul, time loop restricted to calibration period
%    - For now fixed for daily avg output data (change dtstep for other)
% Gabrielle De Lannoy - 18 March 2023:
%    - Added domain subsetting inside LIS-reader to save memory
% Gabrielle De Lannoy - 05 May 2023:
%    - Added optional argument 'pol': 
%       if not given, then DAOBS filenames without specification of the pol
%          are read, otherwise, the filenames with VV or VH spec are read.
%==========================================================================

dtstep   = 86400;
file_tag = '.d01';

%--------------------------------------------------------------------------

% DAOBS file specs
precision= 'float32';
F  = 1;
NaNvalue = -9999.;
if read_obs == 0
    LIS_obs = NaNvalue;
end
pol_tag  = '';
if exist('pol')
  if strcmp(pol, 'VV')
    pol_tag  = '_VV_';
  elseif strcmp(pol, 'VH')
    pol_tag  = '_VH_';
  end
end

%--------------------------------------------------------------------------

Ds = datevec(datenum(start_cal_date));
start_time.year = Ds(:,1);
start_time.month= Ds(:,2);
start_time.day  = Ds(:,3);
start_time.hour = Ds(:,4);
start_time.min  = Ds(:,5);
start_time.sec  = Ds(:,6);
De = datevec(datenum(end_cal_date));
end_time.year = De(:,1);
end_time.month= De(:,2);
end_time.day  = De(:,3);
end_time.hour = De(:,4);
end_time.min  = De(:,5);
end_time.sec  = De(:,6);

date_time = start_time;
t_ind     = 0;
stop      = 0;

while (stop~=1)
    
  if (date_time.year == end_time.year   & ...
      date_time.month== end_time.month  & ...
      date_time.day  == end_time.day    & ...
      date_time.hour == end_time.hour   & ...
      date_time.min  == end_time.min    & ...
      date_time.sec  == end_time.sec          )
    
    stop = 1; %on next loop; first finish this day
    
  end
  
  % Get simulations
  % Daily simulations of today are averaged in a file with the 
  % time stamp of next day
  date_time_tmp = augment_date_time(dtstep, date_time);
  
  YYYYMM = [ num2str(date_time_tmp.year,  '%4.4d'),     ...
             num2str(date_time_tmp.month, '%2.2d') ];    
  DDHHMM = [ num2str(date_time_tmp.day,    '%2.2d'),     ...
             num2str(date_time_tmp.hour,   '%2.2d'),     ...
             num2str(date_time_tmp.min,    '%2.2d') ]; 
        
  fname = [ path, '/SURFACEMODEL/',...
        YYYYMM, '/LIS_HIST_',    ...
        YYYYMM, DDHHMM, file_tag, '.nc'];
    
  if t_ind == 0
    lat_all=ncread(fname,'lat');
    lon_all=ncread(fname,'lon');
    %newer output versions
    if any(size(lat_all)==1) && any(size(lon_all)==1)
        [lon_all,lat_all]=meshgrid(lon_all,lat_all);
        lon_all = lon_all';
        lat_all = lat_all';
    end
    [nx, ny] = size(lat_all); 
    lat_all=lat_all(:);
    lon_all=lon_all(:);
    
    % Potential domain subset
    [domain, lat, lon] = domain_subset(lat_all, lon_all, ...
    min_lat, max_lat, min_lon, max_lon);    

  end
  t_ind = t_ind + 1;

  if strcmp(DA_tag, '.a01')
      S=ncread(fname,'SoilMoist_tavg');
      S=S(:,:,1); %first layer
      LIS_ssm_sd(t_ind,:)=S(domain);
      S=ncread(fname,'LAI_tavg');
      LIS_veg_swe(t_ind,:)=S(domain);
  elseif strcmp(DA_tag, '.a02')
      S=ncread(fname,'SnowDepth_inst'); %m
      LIS_ssm_sd(t_ind,:)=S(domain);
      S=ncread(fname,'SWE_inst'); %kg m-2
      LIS_veg_swe(t_ind,:)=S(domain);
  end
  %check the data
  %scatter(lon,lat,30,nanmean(LIS_ssm,1),'s','filled')

  if read_obs>0
      % Get (mean) observations within the day
      % Time stamp is at the time of assimilation today.
      YYYYMM = [ num2str(date_time.year,   '%4.4d'),     ...
                 num2str(date_time.month,  '%2.2d') ];     
      tmp = NaN(length(domain), 24);  
      for hh=1:24   
          fname  = [ path, '/DAOBS/',...
                 YYYYMM, '/LISDAOBS_',    ...
                 pol_tag, YYYYMM, num2str(date_time.day,'%2.2d'),...
                 num2str(hh,'%2.2d'),'00', DA_tag, file_tag, '.1gs4r'];
          if (exist(fname)==2)
            %The rows/colums and their order are a bit weird for S1 gamma data,
            %but it works out when everything is put in a vector
            [data] = read_LISout1gdr(fname,precision,nx,ny,F,NaNvalue);
            data(data==NaNvalue) = nan;
            tmp(:,hh) = squeeze(data(domain));
          end
      end
      LIS_obs(t_ind,:) = nanmean(tmp,2);
      %check the data
      %scatter(lon,lat,30,nanmean(LIS_obs,1),'s','filled')
  end
  
  YY(t_ind) = date_time.year;
  MM(t_ind) = date_time.month;
  DD(t_ind) = date_time.day;
  date_time = augment_date_time(dtstep, date_time);
    
end

dates = datetime(YY,MM,DD);

% Exclude months
logical_member=ismember(MM',exclude_months);
dates(logical_member)     = [];
LIS_ssm_sd(logical_member,:) = [];  
LIS_veg_swe(logical_member,:) = []; 
LIS_obs(logical_member,:) = [];

end

%============================EOF===========================================
