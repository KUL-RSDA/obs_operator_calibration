%
%Merge calibration output from multiple CPUS and multiple domains
%into 1 global file 
%
%Gabrielle De Lannoy - 21/aug/2013
%aGbrielle De Lannoy - 18/mar/2023: updated for WCM calibration
%                      WARNING: hardwired assumption that A was calibrated
%====================================================================

clear all

addpath('/data/leuven/314/vsc31402/process_code/MATLAB_basics/');

%---------------------------------------------------------------------
runs = {'po_ol_hymap_irr','poE5_ol_hymap_irr',...
        'poE5_ol_hymap_noirr','po_ol_hymap_noirr'};
ncpu = [31 31 30 30];

for rr=1:length(runs)
run = runs{rr};

par_file_path = ['/scratch/leuven/314/vsc31402/4DMED/',run,'/15-apr-2015_15-apr-2023/'];
fname_in      = {'Par_sVV_DREAMZS_ABCDs_01__lat','Par_sVV_DREAMZS_ABCDs_01__lon','Par_sVV_DREAMZS_ABCDs_01__s_0'};
num           = [1:ncpu(rr)]; 

% par_file_path = '/staging/leuven/stg_00024/OUTPUT/projects/4DMED/WCM_calibration/po_ol_no_irr/01-apr-2015_30-dec-2021/';
% fname_in      = {'Par_sVV_DREAMZS_ABCDs_01__lat','Par_sVV_DREAMZS_ABCDs_01__lon','Par_sVV_DREAMZS_ABCDs_01__s_0'};  %'Par_sVV_sceua_ABCD1_'; %'Par_sVV_sceua_ABCDs_01_'; %'Par_sVV_sceua_ABCD0_'; %'Par_sVV_sceua_ABCD1_';
% num           = [1:30]; %[1:30] %[1:9]; %

for fin=1:length(fname_in)
    
c=0;
for i=1:length(num) 
    c=c+1;
    filenames{c} = ['/',num2str(i),'_',num2str(length(num)),'/',fname_in{fin},'.txt'];
end

fname_out     = [par_file_path, '/', fname_in{fin},'_all.txt'];

%--------------------------------------------------------------------

% Simply concatenate

for f=1:length(filenames)
 
    fid = fopen([par_file_path,filenames{f}]);
    if fid>=0
        disp([par_file_path,filenames{f}]);
        tmp = textscan(fid,'%f');
        data = tmp{1};
        fclose(fid);
    end
    
    if f==1
        M = data;
    else
        M = [M; data];
    end
    
end

if ~isempty(M)
  disp(['Saving file...']);
  save (fname_out, 'M', '-ascii');
end

%--------------------------------------------------------------------------  

%Check output
if fin==3
fname     = [par_file_path, '/Par_sVV_DREAMZS_ABCDs_01__lon_all.txt'];
fid = fopen(fname);
tmp = textscan(fid,'%f');
fclose(fid);
lon_all = tmp{1};
fname     = [par_file_path, '/Par_sVV_DREAMZS_ABCDs_01__lat_all.txt'];
fid = fopen(fname);
tmp = textscan(fid,'%f');
fclose(fid);
lat_all = tmp{1};


ps = 1; %pixel size
figure; 
MM = nanmean(M);
scatter(lon_all,lat_all,ps,M,'s','filled');axis tight;
title([run,' s0 mean ',num2str(MM,'%2.2d')]);
colorbar; caxis([0.0 2]);
%
print('-djpeg',[par_file_path, '/figs/', run,fname_in{fin}],'-r300')
end
%
end

clear filenames

end
%-------------------EOF----------------------------------------------------

