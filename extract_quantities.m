%%% Plots instantaneous snapshots
%% setup
clear
close all

data_path = '/nfs/attic/niemeyek/ProCI_Vorticity_paper';

MA_filename = fullfile(data_path, 'data_MA_Vort_2.020E-02');
MC_filename = fullfile(data_path, 'data_MC_Vort_2.020E-02');
configfile = fullfile(data_path, 'config');

addpath 'colormaps';

[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(configfile);
[X,Y,Z] = meshgrid(ym,xm,zm);
L = y(end);

%% MA calculations
filename = MA_filename;
U = NGAdatareader_large(filename,1);
V = NGAdatareader_large(filename,2);
W = NGAdatareader_large(filename,3);

[omegaX,omegaY,omegaZ] = curl(X,Y,Z,U,V,W);

T = NGAdatareader_large(filename,16);
H = NGAdatareader_large(filename,12);

save('vorticity_MA.mat', 'omegaX', 'omegaY', 'omegaZ', 'T', 'H', 'W')


%% MC calculations
filename = MC_filename;
U = NGAdatareader_large(filename,1);
V = NGAdatareader_large(filename,2);
W = NGAdatareader_large(filename,3);

[omegaX,omegaY,omegaZ] = curl(X,Y,Z,U,V,W);

T = NGAdatareader_large(filename,16);
H = NGAdatareader_large(filename,12);

save('vorticity_MC.mat', 'omegaX', 'omegaY', 'omegaZ', 'T', 'H', 'W')

