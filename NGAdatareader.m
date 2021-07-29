% This script reads in data from files generated by NGA
function data = NGAdatareader(filename)

% clear
% clc
% close all
PRODRATES_FLAG = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User modified variables

%filename = '~/Research/H2_Example/1D_Flame_Detailed/phi0_4/thermal_diffusion_comparison/thermal_diffusion/optdata_1.000E-01'

Lx      = 0.048;
Ly      = 1;
Lz      = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid     = fopen(filename,'r');
dims    = fread(fid,4,'integer*4','ieee-le'); % read only the first 4 elements from the file

% Dimensions
nx      = dims(1);
ny      = dims(2);
nz      = dims(3);
nvar    = dims(4);

nsize   = nx*ny*nz;

deltax= Lx/nx;      % mesh size

% The following lines are necessary to move the file pointer to the correct
% location befory you start reading in the variables
fread(fid,1,'real*8','ieee-le'); % skipping over 'dt'
fread(fid,1,'real*8','ieee-le'); % skipping over 'time'
varnames = [];
for var=1:nvar
    varnames = [varnames; fread(fid,8,'*char','ieee-le')]; % skipping over the stored names
end

varnames;

% There are usually 5 variables (or nvar number of variables) stored in the data
% file: U, V, W, P, ZMIX. Each fread call makes the file pointer shift to the end
% of the number of values read. Hence each subsequent fread call will start after
% where the previous fread call stopped. Add on extra variables if necessary.
% %

% if isempty(strfind(filename,'opt'))
%     dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
%     U       = reshape(dummy,nx,ny,nz);                % now turning the column vector into a 3D matrix
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     V       = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     W       = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     P       = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{1}    = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{2}    = reshape(dummy,nx,ny,nz);
%     
%     data = Y;
% else
%     dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
%     SRC_PROG       = reshape(dummy,nx,ny,nz);                % SRC_PROG
%     
%     data = SRC_PROG;
% end
if nvar == 4
    dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
    U       = reshape(dummy,nx,ny,nz);                % now turning the column vector into a 3D matrix
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    V       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    W       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    P       = reshape(dummy,nx,ny,nz);
    data = {U,V,W,P};
elseif isempty(strfind(filename,'opt')) && nvar > 10 && isempty(strfind(filename,'table'))
    dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
    U       = reshape(dummy,nx,ny,nz);                % now turning the column vector into a 3D matrix
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    V       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    W       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    P       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    RHO       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    dRHO       = reshape(dummy,nx,ny,nz);
    for i=1:nvar-6
        dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
        Y{i}       = reshape(dummy,nx,ny,nz);                % N2
    end
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{2}       = reshape(dummy,nx,ny,nz);                % H
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{3}       = reshape(dummy,nx,ny,nz);                % O2
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{4}       = reshape(dummy,nx,ny,nz);                % O
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{5}       = reshape(dummy,nx,ny,nz);                % OH
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{6}       = reshape(dummy,nx,ny,nz);                % H2
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{7}       = reshape(dummy,nx,ny,nz);                % H2O
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{8}       = reshape(dummy,nx,ny,nz);                % HO2
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{9}       = reshape(dummy,nx,ny,nz);                % H2O2
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{10}       = reshape(dummy,nx,ny,nz);                % T
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{11}       = reshape(dummy,nx,ny,nz);                % ZMIX
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     Y{12}       = reshape(dummy,nx,ny,nz);                % VELGRAD
    
    if PRODRATES_FLAG == 0
        for i=1:nvar-6
            Y{i}(Y{i}<0) = 0;
        end
    end
    data = {U,V,W,P,RHO,dRHO,Y};
    
elseif strfind(filename,'optdata')% Read OPTDATA
    for i=1:nvar
        dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
        data{i}       = reshape(dummy,nx,ny,nz);                % N2
    end
    
elseif isempty(strfind(filename,'opt'))
    dummy   = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
    U       = reshape(dummy,nx,ny,nz);                % now turning the column vector into a 3D matrix
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    V       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    W       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    P       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    RHO       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    dRHO       = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    ZMIX    = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    PROG    = reshape(dummy,nx,ny,nz);
    data = {U,V,W,P,RHO,dRHO,ZMIX,PROG};
%     data = {U,V,W,P,RHO,dRHO,PROG};
    
else
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     SRC_PROG    = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     T    = reshape(dummy,nx,ny,nz);
%     dummy   = fread(fid,nsize,'real*8','ieee-le');
%     DIFFPROG    = reshape(dummy,nx,ny,nz);
%     data = {SRC_PROG,T,DIFFPROG};
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    srcTdiff    = reshape(dummy,nx,ny,nz);
    dummy   = fread(fid,nsize,'real*8','ieee-le');
    srcTenth    = reshape(dummy,nx,ny,nz);
    data = {srcTdiff,srcTenth};
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now use U, V, W, P, ZMIX for computations you need
%end