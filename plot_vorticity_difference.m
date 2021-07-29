%%% Plots instantaneous snapshots
%% setup
clear
close all

configfile = 'config';
addpath 'colormaps';

[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(configfile);
[X,Y,Z] = meshgrid(ym,xm,zm);
L = y(end);

%% load data

load vorticity_MA.mat;
vorticity_MA = sqrt(omegaX.^2 + omegaY.^2 + omegaZ.^2);

load vorticity_MC.mat;
vorticity_MC = sqrt(omegaX.^2 + omegaY.^2 + omegaZ.^2);

vorticity_diff = log10(abs(vorticity_MA - vorticity_MC));

%% Plot

[xg,yg,zg] = meshgrid(x(1:end-1),y(1:end-1),z(1:end-1));

for i=1:nx
    xm(i) = (x(i+1)+x(i))/2;
end
for j=1:ny
    ym(j) = (y(j+1)+y(j))/2;
end
for k=1:nz
    zm(k) = (z(k+1)+z(k))/2;
end

startval = 1520/4; endval = 3*1520/4;

Tpeak = 1185;
Tu = 400;
T_slice_1(:,:) = T(:,:,120);
isoline_1 = contourc(T_slice_1',[Tpeak+300 Tpeak+300]);
checkedlength = 0;
indexcheck = 1;
while checkedlength < length(isoline_1)
   isolength = isoline_1(2, indexcheck);
   isoline_1(:, indexcheck) = []; % Remove column with isovalue and # points
   checkedlength = checkedlength + isolength;
   indexcheck = checkedlength + 1;
end
isoline_2 = contourc(T_slice_1',[Tpeak-300 Tpeak-300]);
checkedlength = 0;
indexcheck = 1;
while checkedlength < length(isoline_2)
   isolength = isoline_2(2, indexcheck);
   isoline_2(:, indexcheck) = []; % Remove column with isovalue and # points
   checkedlength = checkedlength + isolength;
   indexcheck = checkedlength + 1;
end


h = figure();
%set(gcf,'units','centimeters','position',[0,0,6.7*4,6.7*2]);
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L, vorticity_diff(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image;
scatter(isoline_1(1,:).*dx/L, isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L, isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);
xlabel('$x/L$','interpreter','latex','FontSize',16);
ylabel('$y/L$','interpreter','latex','FontSize',18);
ax=gca;
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); 
set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',18);
shading interp
colorbar

%% Save

exportgraphics(h, 'vorticity_difference.png', 'Resolution', 300)
%exportgraphics(h, 'vorticity_difference.pdf', 'ContentType', 'vector')