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
fname_in      = {'Par_sVV_DREAMZS_ABCDs_01_'};  %'Par_sVV_sceua_ABCD1_'; %'Par_sVV_sceua_ABCDs_01_'; %'Par_sVV_sceua_ABCD0_'; %'Par_sVV_sceua_ABCD1_';
num           = [1:ncpu(rr)]; 

% par_file_path = '/staging/leuven/stg_00024/OUTPUT/projects/4DMED/WCM_calibration/po_ol_irr/01-apr-2015_30-dec-2021/';
% fname_in      = {'Par_sVV_DREAMZS_ABCD0_','Par_sVV_DREAMZS_ABCD1_','Par_sVV_DREAMZS_ABCDs_01_'};  %'Par_sVV_sceua_ABCD1_'; 
% num           = [1:30]; 

for fin=1:length(fname_in)
    
c=0;
for i=1:length(num) 
    c=c+1;
    filenames{c} = ['/',num2str(i),'_',num2str(length(num)),'/',fname_in{fin},'.mat'];
end

fname_out     = [par_file_path, '/', fname_in{fin},'all.mat'];

%--------------------------------------------------------------------

%Get structure from 'first' param file

disp([par_file_path,filenames{1}]);

all_fields = load([par_file_path,filenames{1}]);

fn_all  = fieldnames(all_fields);

for i=1:length(fn_all)
    cmd = [fn_all{i},' = all_fields.(fn_all{i});'];
    eval(cmd);
end

%--------------------------------------------------------------------

%Keep everything from the 'first'/previous param file, and replace entries
%where the new file has calibrated parameters while the 'first'/previous does not.

for f=2:length(filenames)

    disp([par_file_path,filenames{f}]);
 
    TMP=load([par_file_path,filenames{f}]);
    if (any(TMP.lat_all-lat_all) || any(TMP.lon_all-lon_all) || ...
            ~any(strcmp(fn_all,fieldnames(TMP))) )
        error('different files have different lat-lon or fields');
    end

    %Find calibrated parameters
    %!!CAREFUL, change this if A happened to not be calibrated!
    ind = find(~isnan(TMP.A_cal_all));

    for ii = 1:length(fn_all)

        new_param = TMP.(fn_all{ii});
        cmd = [fn_all{ii},'(ind,:) = new_param(ind,:);'];
        eval(cmd);

    end

end

clear filenames

disp(['Saving file...']);

%Copied from CalWCM_main.m
save (fname_out,'lat_all', 'lon_all',...
      'A_cal_all','B_cal_all','C_cal_all','D_cal_all',...
      'A_std_all','B_std_all','C_std_all','D_std_all','-v7.3');    

%--------------------------------------------------------------------------  
  
%Check output
if 1
ps = 2; %pixel size
figure; 
subplot(2,2,1)
scatter(lon_all,lat_all,ps,A_cal_all,'s','filled');axis tight;
title(['WCM(A) ',num2str(nanmean(A_cal_all),'%2.2d')]);
colorbar; caxis([0.0 1.0]);
subplot(2,2,2)
scatter(lon_all,lat_all,ps,B_cal_all,'s','filled');axis tight;
title(['WCM(B) ',num2str(nanmean(B_cal_all),'%2.2d')]);
colorbar; caxis([0.0 0.5]);
subplot(2,2,3)
scatter(lon_all,lat_all,ps,C_cal_all,'s','filled');axis tight;
title(['WCM(C) ',num2str(nanmean(C_cal_all),'%2.2d')]);
colorbar; caxis([-35 -10]);
subplot(2,2,4)
scatter(lon_all,lat_all,ps,D_cal_all,'s','filled');axis tight;
title(['WCM(D) ',num2str(nanmean(D_cal_all),'%2.2d')]);
colorbar; caxis([15 80]);
%
print('-djpeg',[par_file_path, '/figs/', fname_in{fin}],'-r300')
end
%
end
end
%-------------------EOF----------------------------------------------------

