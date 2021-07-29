%%% Plots instantaneous snapshots
%% setup
clear
close all

configfile = 'config';
addpath 'colormaps';

[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(configfile);
[X,Y,Z] = meshgrid(ym,xm,zm);
L = y(end);

% check for presence of large files
if ~isfile('vorticity_MA.mat') || ~isfile('vorticity_MC.mat')
    disp('Error: vorticity_MA.mat and/or vorticiy_MC.mat are missing.')
    disp('Check the README for details on obtaining')
    return
end

%% MA data

load vorticity_MA.mat;
vorticity = log10(sqrt(omegaX.^2 + omegaY.^2 + omegaZ.^2));

%% Plots
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

t = tiledlayout(4,2);

ax = nexttile;
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L,T(:,:,120),256,'EdgeColor','none'),hold on ; 
colormap viridis; 
axis image; %startval:endval
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%'fontweight','bold','fontname','Times New Roman')
ylabel('$y/L$','Interpreter','latex','FontSize',18);
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar

title('(a) MA temperature [K]', 'Interpreter','latex');

ax = nexttile(3);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L,H(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image; 
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
shading interp
set(gcf,'color','w')
set(gcf,'units','centimeters')
cbh=colorbar;
set(cbh,'YTick',[1*10^-3 5*10^-3 10*10^-3]);

title('(c) MA $Y_{\rm{H}_2}$', 'Interpreter','latex');


ax = nexttile(5);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L, W(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image;
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar

title('(e) MA $u_z$ [m/s]', 'Interpreter','latex');

ax = nexttile(7);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L, vorticity(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image;
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
xlabel('$x/L$','Interpreter','latex');
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar
cbh=colorbar;
set(cbh,'YTick',[0 1 2 3 4 5])

title('(g) MA $\log_{10}(\omega^2)$', 'Interpreter','latex');

%% MC calculations

load 'vorticity_MC.mat';
vorticity = log10(sqrt(omegaX.^2 + omegaY.^2 + omegaZ.^2));

%% Plots

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


ax = nexttile(2);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L,T(:,:,120),256,'EdgeColor','none'),hold on ; 
colormap viridis; 
axis image; %startval:endval
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar

title('(b) MC temperature [K]', 'Interpreter','latex');

ax = nexttile(4);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L,H(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image; 
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
shading interp
set(gcf,'color','w')
set(gcf,'units','centimeters')
cbh=colorbar;
set(cbh,'YTick',[1*10^-3 5*10^-3 10*10^-3])

title('(d) MC $Y_{\rm{H}_2}$', 'Interpreter','latex');

ax = nexttile(6);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L, W(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image;
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar

title('(f) MA $u_z$ [m/s]', 'Interpreter','latex');

ax = nexttile(8);
%set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7]);
% c1 = subplot(2,1,1,'units','centimeters')
contourf(xg(:,:,1)'/L,yg(:,:,1)'/L, vorticity(:,:,120),256,'EdgeColor','none');
hold on; 
colormap viridis; 
axis image;
scatter(isoline_1(1,:).*dx/L,isoline_1(2,:).*dx/L,1,'r.')
scatter(isoline_2(1,:).*dx/L,isoline_2(2,:).*dx/L,1,'w.')
axis([180*dx/L 890*dx/L 0 190*dx/L])
set(gca,'FontSize',16,'linewidth',1);%,'fontweight','bold','fontname','Times New Roman')
xlabel('$x/L$','Interpreter','latex');
ylabel('$y/L$','Interpreter','latex','FontSize',18);%,'fontweight','bold');
set(ax,'XTick',[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]); set(ax,'YTick',[0 0.5 1]);
pos=get(gca,'pos');
set(gcf,'color','w')
set(gcf,'units','centimeters')
set(gca,'FontSize',16);
shading interp
colorbar
cbh=colorbar;
set(cbh,'YTick',[0 1 2 3 4 5])

title('(h) MC $\log_{10}(\omega^2)$', 'Interpreter','latex');

%% save

exportgraphics(t, 'instantaneous_fields.png', 'Resolution', 300)
%exportgraphics(t, 'instantaneous_fields.pdf', 'ContentType', 'vector')