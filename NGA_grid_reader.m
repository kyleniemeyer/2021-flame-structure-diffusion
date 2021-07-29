function [x_t,y_t,z_t,xm_t,ym_t,zm_t,nx_t,ny_t,nz_t,dx,dy,dz] = NGA_grid_reader(srcfilename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User modified variables

% srcfilename = '~/Research/H2_Example/2D_Flame/phi0_4/24mm_domain_refined_long/config';
% Also modify the variables to be written later on in the file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the mesh
fid = fopen(srcfilename,'r');

simulation = fread(fid, 64, '*char', 'ieee-le');

icyl = fread(fid, 1, 'integer*4', 'ieee-le');

xper = fread(fid, 1, 'integer*4', 'ieee-le');
yper = fread(fid, 1, 'integer*4', 'ieee-le');
zper = fread(fid, 1, 'integer*4', 'ieee-le');

nx_t = fread(fid, 1, 'integer*4', 'ieee-le');
ny_t = fread(fid, 1, 'integer*4', 'ieee-le');
nz_t = fread(fid, 1, 'integer*4', 'ieee-le');

x_t = fread(fid, nx_t+1,'real*8', 'ieee-le');
y_t = fread(fid, ny_t+1,'real*8', 'ieee-le');
z_t = fread(fid, nz_t+1,'real*8', 'ieee-le');

mask = fread(fid, nx_t*ny_t,'integer*4', 'ieee-le');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx, dy, and dz in constant region

% dx = min(diff(x_t));
% dy = min(diff(y_t));
% dz = min(diff(z_t));

if nx_t==1
    dx=x_t(end)-x_t(1);
else
	dx=x_t(ceil(nx_t/2))-x_t(ceil(nx_t/2)-1);
end

if ny_t==1
    dy=y_t(end)-y_t(1);
else
    dy=y_t(ceil(ny_t/2))-y_t(ceil(ny_t/2)-1);
end


if nz_t==1
    dz=z_t(end)-z_t(1);
else
    dz=z_t(ceil(nz_t/2))-z_t(ceil(nz_t/2)-1);
end

fclose(fid);

for i=1:nx_t
    xm_t(i) = (x_t(i+1)+x_t(i))/2;
end
xm_t = xm_t';

for j=1:ny_t
    ym_t(j) = (y_t(j+1)+y_t(j))/2;
end
ym_t = ym_t';

for k=1:nz_t    
    zm_t(k) = (z_t(k+1)+z_t(k))/2;
end
zm_t = zm_t';

end
