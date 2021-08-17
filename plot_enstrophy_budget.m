%%% Plots enstrophy budget

%% Setup

clear
close all

configfile = 'config';

[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(configfile);

base_path = 'enstrophy_terms';

%% Mixture averaged
num_files = 23;

baroclinic_MA_avg = 0;
dilatation_MA_avg = 0;
viscous_MA_avg = 0;
forcing_MA_avg = 0;
stretch_MA_avg = 0;
enstrophy_MA_avg = 0;

for i = 2 : num_files
    %filename = fullfile(base_path, join(['avg_', int2str(i), '__MA.mat']));
    filename = fullfile(base_path, join(['enstrophy_', int2str(i), '_MA.mat']));
    load(filename);
    baroclinic_MA_avg = baroclinic_MA_avg + baroclinic;
    dilatation_MA_avg = dilatation_MA_avg + dilatation;
    viscous_MA_avg = viscous_MA_avg + viscous_effects;
    forcing_MA_avg = forcing_MA_avg + forcing;
    stretch_MA_avg = stretch_MA_avg + stretch;
    enstrophy_MA_avg = enstrophy_MA_avg + enstrophy;
end

baroclinic_MA_avg = baroclinic_MA_avg / num_files;
dilatation_MA_avg = dilatation_MA_avg / num_files;
viscous_MA_avg = viscous_MA_avg / num_files;
forcing_MA_avg = forcing_MA_avg / num_files;
stretch_MA_avg = stretch_MA_avg / num_files;
enstrophy_MA_avg = enstrophy_MA_avg / num_files;

%% Multicomponent
num_files = 25;

baroclinic_MC_avg = 0;
dilatation_MC_avg = 0;
viscous_MC_avg = 0;
forcing_MC_avg = 0;
stretch_MC_avg = 0;
enstrophy_MC_avg = 0;

for i = 2 : num_files
    %filename = fullfile(base_path, join(['avg_', int2str(i), '__MC.mat']));
    filename = fullfile(base_path, join(['enstrophy_', int2str(i), '_MC.mat']));
    load(filename);
    baroclinic_MC_avg = baroclinic_MC_avg + baroclinic;
    dilatation_MC_avg = dilatation_MC_avg + dilatation;
    viscous_MC_avg = viscous_MC_avg + viscous_effects;
    forcing_MC_avg = forcing_MC_avg + forcing;
    stretch_MC_avg = stretch_MC_avg + stretch;
    enstrophy_MC_avg = enstrophy_MC_avg + enstrophy;
end

baroclinic_MC_avg = baroclinic_MC_avg / num_files;
dilatation_MC_avg = dilatation_MC_avg / num_files;
viscous_MC_avg = viscous_MC_avg / num_files;
forcing_MC_avg = forcing_MC_avg / num_files;
stretch_MC_avg = stretch_MC_avg / num_files;
enstrophy_MC_avg = enstrophy_MC_avg / num_files;

space = ones(length(xm),1).*xm;

%% Plot 

h = figure;
viscous_MC_avg_x = movmean(viscous_MC_avg,30);
viscous_MA_avg_x = movmean(viscous_MA_avg,30);
H = plot(xm./0.00643,viscous_MA_avg_x,'--','linewidth',2);
Colors(:,1) = get(H,'Color');
hold on
plot(xm./0.00631,viscous_MC_avg_x,'linewidth',2,'Color',Colors(:,1))


baroclinic_MC_avg_x = movmean(baroclinic_MC_avg,30);
baroclinic_MA_avg_x = movmean(baroclinic_MA_avg,30);
H = plot(xm./0.00643,baroclinic_MA_avg_x,'--','linewidth',2);
Colors(:,2) = get(H,'Color');
hold on
plot(xm./0.00631,baroclinic_MC_avg_x,'linewidth',2,'Color',Colors(:,2))

dilatation_MC_avg_x = movmean(dilatation_MC_avg,30);
dilatation_MA_avg_x = movmean(dilatation_MA_avg,30);
H = plot(xm./0.00643,dilatation_MA_avg_x,'--','linewidth',2);
Colors(:,3) = get(H,'Color');
hold on;
plot(xm./0.00631,dilatation_MC_avg_x,'linewidth',2,'Color',Colors(:,3))

forcing_MC_avg_x = movmean(forcing_MC_avg,30);
forcing_MA_avg_x = movmean(forcing_MA_avg,30);%*200000
H = plot(xm./0.00643,forcing_MA_avg_x,'k--','linewidth',2);
hold on
plot(xm./0.00631,forcing_MC_avg_x,'k','linewidth',2)
% Colors(:,1) = get(H,'Color');

stretch_MC_avg_x = movmean(stretch_MC_avg,30);
stretch_MA_avg_x = movmean(stretch_MA_avg,30);
H = plot(xm./0.00643,stretch_MA_avg_x,'r--','linewidth',2);
hold on
plot(xm./0.00631,stretch_MC_avg_x,'r','linewidth',2)

%LHS_MC_avg_x = movmean(viscous_MC_avg+baroclinic_MC_avg-dilatation_MC_avg+forcing_MC_avg+stretch_MC_avg,30); %movmean(LHS_MC_avg,30);
%LHS_MA_avg_x = movmean(viscous_MA_avg+baroclinic_MA_avg-dilatation_MA_avg+forcing_MA_avg+stretch_MA_avg,30); %movmean(LHS_MA_avg,30);-Dissipation_MA_avg+Baroclinic_MA_avg+Dilatation_MA_avg+Forcing_MA_avg+Stretch_MA_avg;
%plot(xm./0.00631,LHS_MC_avg_x,'m','linewidth',2),hold on
%plot(xm./0.00643,LHS_MA_avg_x,'m--','linewidth',2)

%LHS_MC_avg_x = movmean(enstrophy_MC_avg,30);
%LHS_MA_avg_x = movmean(enstrophy_MA_avg,30);
%H = plot(xm./0.00643,LHS_MA_avg_x,'r--','linewidth',2);
%hold on
%plot(xm./0.00631,LHS_MC_avg_x,'r','linewidth',2)

set(gca,'FontSize',14,'linewidth',1)
xlabel('$x/\delta_L$','interpreter','latex','fontsize',16)
ylabel('Enstrophy $[1/s^3]$','interpreter','latex','fontsize',16)
axis square
axis([0 10 -1.75e13 1.75e13])
set(gcf,'color','w')
% legend('MC H_{2}','MA H_{2}','MC \itn\rm-C_{7}H_{16}','MA \itn\rm-C_{7}H_{16}','MC A_{1}CH_{3}','MA A_{1}CH_{3}','location','ne')
% set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
legend('Viscous MA','Viscous MC','Baroclinic MA','Baroclinic MC',...
    'Dilatation MA','Dilatation MC','Forcing MA','Forcing MC',...
    'Stretch MA','Stretch MC','interpreter','latex',...
    'location','ne','fontsize',16);
legend boxoff
set(gcf,'color','w')
set(gcf,'units','centimeters')
pos = get(gcf,'position');
AR = pos(3)/pos(4);
set(gcf,'units','centimeters','position',[0,0,6.7*3,6.7*3]);
grid on

exportgraphics(h, 'enstrophy_budget_dim.pdf', 'ContentType', 'vector')

%% enstrophy

h = figure;

enstrophy_MC_avg_x = movmean(enstrophy_MC_avg,30);
enstrophy_MA_avg_x = movmean(enstrophy_MA_avg,30);
H = plot(xm./0.00643,enstrophy_MA_avg_x,'r--','linewidth',2,'Color','b');
hold on
plot(xm./0.00631,enstrophy_MC_avg_x,'r','linewidth',2,'Color','b')

set(gca,'FontSize',14,'linewidth',1)
xlabel('$x/\delta_L$','interpreter','latex','fontsize',16)
ylabel('Enstrophy $[1/s^3]$','interpreter','latex','fontsize',16)
axis square
axis([0 10 0 2.5e9])
set(gcf,'color','w')
legend('MA','MC','interpreter','latex','location','ne','fontsize',16);
legend boxoff
set(gcf,'color','w')
set(gcf,'units','centimeters')
pos = get(gcf,'position');
AR = pos(3)/pos(4);
set(gcf,'units','centimeters','position',[0,0,6.7*3,6.7*3]);
grid on

exportgraphics(h, 'enstrophy.pdf', 'ContentType', 'vector')

%% nondimensional plot

RHO_u = 1.0224;
RHO_b = 0.2333;
RHO_avg = 0.4261;
gamma = (RHO_u-RHO_b)/RHO_avg;


S_L_MC = 0.223;
S_L_MA = 0.230;

l_f_MC = 0.000631;
l_f_MA = 0.000643;

% tau_MC = 1.816e-5;
% tau_MA = 1.876e-5;

tau_f_MC = l_f_MC/S_L_MC;
tau_f_MA = l_f_MA/S_L_MA;

tau_eta_MC = (tau_f_MC/151);
tau_eta_MA = (tau_f_MA/149);

u_p_MC = 18.6*S_L_MC;
u_p_MA = 18*S_L_MA;

viscous_MC_scale = .1/tau_eta_MC^3;
viscous_MA_scale = .1/tau_eta_MA^3;

stretch_MC_scale =  .1/tau_eta_MC^3;
stretch_MA_scale =  .1/tau_eta_MA^3;
 
dilatation_MC_scale = gamma*(10*S_L_MC/(l_f_MC*tau_eta_MC^2));
dilatation_MA_scale = gamma*(10*S_L_MA/(l_f_MA*tau_eta_MA^2));

baroclinic_MC_scale = gamma*(u_p_MC/(l_f_MC*tau_eta_MC^2));
baroclinic_MA_scale = gamma*(u_p_MA/(l_f_MA*tau_eta_MA^2));

forcing_MC_scale = 1/((1/(2*973.05))*tau_eta_MC^2);
forcing_MA_scale = 1/((1/(2*973.05))*tau_eta_MA^2);

h = figure;
viscous_MC_avg_x = movmean(viscous_MC_avg,30);
viscous_MA_avg_x = movmean(viscous_MA_avg,30);
baroclinic_MC_avg_x = movmean(baroclinic_MC_avg,30);
baroclinic_MA_avg_x = movmean(baroclinic_MA_avg,30);
dilatation_MC_avg_x = movmean(dilatation_MC_avg,30);
dilatation_MA_avg_x = movmean(dilatation_MA_avg,30);
forcing_MC_avg_x = movmean(forcing_MC_avg,30);%
forcing_MA_avg_x = movmean(forcing_MA_avg,30);%*5000
stretch_MC_avg_x = movmean(stretch_MC_avg,30);
stretch_MA_avg_x = movmean(stretch_MA_avg,30);

plot(xm./0.00643,viscous_MA_avg_x./viscous_MA_scale,'--','Color',Colors(:,1),'linewidth',2); hold on
plot(xm./0.00631,viscous_MC_avg_x./viscous_MC_scale,'Color',Colors(:,1),'linewidth',2);

plot(xm./0.00643,baroclinic_MA_avg_x./baroclinic_MA_scale,'--','Color',Colors(:,2),'linewidth',2);hold on
plot(xm./0.00631,baroclinic_MC_avg_x./baroclinic_MC_scale,'Color',Colors(:,2),'linewidth',2)

plot(xm./0.00643,dilatation_MA_avg_x./dilatation_MA_scale,'--','Color',Colors(:,3),'linewidth',2),hold on
plot(xm./0.00631,dilatation_MC_avg_x./dilatation_MC_scale,'Color',Colors(:,3),'linewidth',2)

plot(xm./0.00643,forcing_MA_avg_x./forcing_MA_scale,'k--','linewidth',2),hold on
plot(xm./0.00631,forcing_MC_avg_x./forcing_MC_scale,'k','linewidth',2)

plot(xm./0.00643,stretch_MA_avg_x./stretch_MA_scale,'r--','linewidth',2),hold on
plot(xm./0.00631,stretch_MC_avg_x./stretch_MC_scale,'r','linewidth',2)


set(gca,'FontSize',14,'linewidth',1)
xlabel('$x/\delta_L$','interpreter','latex','fontsize',16)
ylabel('Normalized enstrophy budget','interpreter','latex','fontsize',16)
axis square
xlim = get(gca,'xlim');
axis([0 10 -0.7 1.05])
yticks([-0.6 -0.3 0 0.3 0.6 0.9])
set(gcf,'color','w')
legend('Viscous MA','Viscous MC','Baroclinic MA','Baroclinic MC',...
    'Dilatation MA','Dilatation MC','Forcing MA','Forcing MC',...
    'Stretch MA','Stretch MC','location','ne','interpreter','latex','fontsize',16)
legend boxoff
set(gcf,'color','w')
set(gcf,'units','centimeters')
grid on
pos = get(gcf,'position');
AR = pos(3)/pos(4);
set(gcf,'units','centimeters','position',[0,0,6.7*3,6.7*3]);

exportgraphics(h, 'enstrophy_budget_nondim.pdf', 'ContentType', 'vector')