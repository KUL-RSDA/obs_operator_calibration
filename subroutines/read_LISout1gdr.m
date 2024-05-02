% Read and plot a LIS direct/sequential access-file
%
% EXAMPLE:
% fname = '/scratch/leuven/320/vsc32041/LIS_chaco/NLDAS2-a/input/PARAMETERS/UMD/1KM/landcover_IGBP_NCEP.1gd4r';
% nx = 36000;  %number of columns
% ny = 15000;  %number of rows
% F  = 1;      %which field to read 
% NaNvalue = -9999;
% precision= 'int32';
% 
% Gabrielle De Lannoy - 19/12/2017
% Michel Bechtold - 02/04/2021  % adapted for 1gs4r, there was one integer value more in the binary files that needs to be skipped 
%=========================================================================

function [data] = read_LISout1gdr(fname,precision,nx,ny,F,NaNvalue) 

tol = 0.0001;

disp(['reading from ', fname])

ifp = fopen( fname, 'r', 'b' );

int_precision   = 'int32';      % precision of fortran tag

for i=1:F
  % if 1gsd4r, skip one integer value
  if strcmp(fname(end-4:end),'1gs4r')
      tmp = fread( ifp,  1, int_precision );
  end
  data = fread( ifp, [nx ny], precision );
end
 
fclose(ifp);

%=========================================================================








