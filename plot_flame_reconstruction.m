%%% Plots flame width and reconstruction

%% Setup

clear, clc
close all

base_path = 'conditional_means';

ref_MC = importdata(fullfile(base_path, 'Condmean_Test_MC'));
ref_lam_MC = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC', 'Condmean_Test'));
ref_MA = importdata(fullfile(base_path, 'Condmean_Test_MA'));
ref_lam_MA = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA', 'Condmean_Test'));

ref_MC_norm = ref_MC;
ref_lam_MC_norm = ref_lam_MC;
ref_MA_norm = ref_MA;
ref_lam_MA_norm = ref_lam_MA;

ref_MC_norm(:,2) = ref_MC(:,2) / ref_MC(1,2);
ref_lam_MC_norm(:,2) = ref_lam_MC(:,2) / ref_lam_MC(1,2);
ref_MA_norm(:,2) = ref_MA(:,2) / ref_MA(1,2);
ref_lam_MA_norm(:,2) = ref_lam_MA(:,2) / ref_lam_MA(1,2);

delta_L_MC = 0.000631;
delta_L_MA = 0.000643;

gradN2_on_T     = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradN2_on_T'));
gradH_on_T      = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradH_on_T'));
gradO2_on_T     = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradO2_on_T'));
gradO_on_T      = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradO_on_T'));
gradOH_on_T     = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradOH_on_T'));
gradH2_on_T     = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradH2_on_T'));
gradH2O_on_T    = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradH2O_on_T'));
gradHO2_on_T    = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradHO2_on_T'));
gradH2O2_on_T   = importdata(fullfile(base_path, 'Condmean_on_T', 'Condmean_gradH2O2_on_T'));

gradT_on_T_MC     = importdata(fullfile(base_path, 'Condmean_on_H2', 'Condmean_gradT_on_T'));
gradT_on_T_lam_MC = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC', 'Condmean_gradT_on_T'));
gradT_on_T_MA     = importdata(fullfile(base_path, 'Condmean_on_H2_MA', 'Condmean_gradT_on_T'));
gradT_on_T_lam_MA = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA', 'Condmean_gradT_on_T'));

gradN2_on_N2      = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradN2_on_N2'));
gradH_on_H        = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradH_on_H'));
gradO2_on_O       = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradO2_on_O2'));
gradO_on_O        = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradO_on_O'));
gradOH_on_OH      = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradOH_on_OH'));
gradH2_on_H2      = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradH2_on_H2'));
gradH2O_on_H2O    = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradH2O_on_H2O'));
gradHO2_on_HO2    = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradHO2_on_HO2'));
gradH2O2_on_H2O2  = importdata(fullfile(base_path, 'Condmean_on_Species', 'Condmean_gradH2O2_on_H2O2'));

inv_GradN2_on_N2     = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradN2_on_N2'));
inv_GradH_on_H       = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradH_on_H'));
inv_GradO2_on_O2     = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradO2_on_O2'));
inv_GradO_on_O       = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradO_on_O'));
inv_GradOH_on_OH     = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradOH_on_OH'));
inv_GradH2_on_H2     = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradH2_on_H2'));
inv_GradH2O_on_H2O   = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradH2O_on_H2O'));
inv_GradHO2_on_HO2   = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradHO2_on_HO2'));
inv_GradH2O2_on_H2O2 = importdata(fullfile(base_path, 'Condmean_on_Species','Condmean_inv_gradH2O2_on_H2O2'));

gradN2_on_H2      = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradN2_on_H2'));
gradH_on_H2_MC    = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradH_on_H2'));
gradO2_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradO2_on_H2'));
gradO_on_H2_MC    = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradO_on_H2'));
gradOH_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradOH_on_H2'));
gradH2_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradH2_on_H2'));
gradH2O_on_H2_MC  = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradH2O_on_H2'));
gradHO2_on_H2_MC  = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradHO2_on_H2'));
gradH2O2_on_H2_MC = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_gradH2O2_on_H2'));

gradN2_on_lam_H2      = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradN2_on_H2'));
gradH_on_H2_lam_MC    = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradH_on_H2'));
gradO2_on_H2_lam_MC   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradO2_on_H2'));
gradO_on_H2_lam_MC    = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradO_on_H2'));
gradOH_on_H2_lam_MC   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradOH_on_H2'));
gradH2_on_H2_lam_MC   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradH2_on_H2'));
gradH2O_on_H2_lam_MC  = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradH2O_on_H2'));
gradHO2_on_H2_lam_MC  = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradHO2_on_H2'));
gradH2O2_on_H2_lam_MC = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_gradH2O2_on_H2'));

gradN2_on_lam_H2      = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradN2_on_H2'));
gradH_on_H2_lam_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradH_on_H2'));
gradO2_on_H2_lam_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradO2_on_H2'));
gradO_on_H2_lam_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradO_on_H2'));
gradOH_on_H2_lam_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradOH_on_H2'));
gradH2_on_H2_lam_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradH2_on_H2'));
gradH2O_on_H2_lam_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradH2O_on_H2'));
gradHO2_on_H2_lam_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradHO2_on_H2'));
gradH2O2_on_H2_lam_MA = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_gradH2O2_on_H2'));

gradH_on_H2_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradH_on_H2'));
gradO2_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradO2_on_H2'));
gradO_on_H2_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradO_on_H2'));
gradOH_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradOH_on_H2'));
gradH2_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradH2_on_H2'));
gradH2O_on_H2_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradH2O_on_H2'));
gradHO2_on_H2_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradHO2_on_H2'));
gradH2O2_on_H2_MA = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_gradH2O2_on_H2'));

inv_GradN2_on_H2      = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradN2_on_H2'));
inv_GradH_on_H2_MC    = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradH_on_H2'));
inv_GradO2_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradO2_on_H2'));
inv_GradO_on_H2_MC    = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradO_on_H2'));
inv_GradOH_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradOH_on_H2'));
inv_GradH2_on_H2_MC   = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradH2_on_H2'));
inv_GradH2O_on_H2_MC  = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradH2O_on_H2'));
inv_GradHO2_on_H2_MC  = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradHO2_on_H2'));
inv_GradH2O2_on_H2_MC = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradH2O2_on_H2'));

inv_GradH_on_H2_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradH_on_H2'));
inv_GradO2_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradO2_on_H2'));
inv_GradO_on_H2_MA    = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradO_on_H2'));
inv_GradOH_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradOH_on_H2'));
inv_GradH2_on_H2_MA   = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradH2_on_H2'));
inv_GradH2O_on_H2_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradH2O_on_H2'));
inv_GradHO2_on_H2_MA  = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradHO2_on_H2'));
inv_GradH2O2_on_H2_MA = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradH2O2_on_H2'));

inv_GradT_on_T_MC     = importdata(fullfile(base_path, 'Condmean_on_H2','Condmean_inv_gradT_on_T'));
inv_GradT_on_T_lam_MC = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MC','Condmean_inv_gradT_on_T'));
inv_GradT_on_T_MA     = importdata(fullfile(base_path, 'Condmean_on_H2_MA','Condmean_inv_gradT_on_T'));
inv_GradT_on_T_lam_MA = importdata(fullfile(base_path, 'Condmean_on_H2_lam_MA','Condmean_inv_gradT_on_T'));

gradN2_on_C     = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradN2_on_C'));
gradH_on_C      = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradH_on_C'));
gradO2_on_C     = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradO2_on_C'));
gradO_on_C      = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradO_on_C'));
gradOH_on_C     = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradOH_on_C'));
gradH2_on_C     = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradH2_on_C'));
gradH2O_on_C    = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradH2O_on_C'));
gradHO2_on_C    = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradHO2_on_C'));
gradH2O2_on_C   = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_gradH2O2_on_C'));

inv_GradN2_on_C   = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradN2_on_C'));
inv_GradH_on_C    = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradH_on_C'));
inv_GradO2_on_C   = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradO2_on_C'));
inv_GradO_on_C    = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradO_on_C'));
inv_GradOH_on_C   = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradOH_on_C'));
inv_GradH2_on_C   = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradH2_on_C'));
inv_GradH2O_on_C  = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradH2O_on_C'));
inv_GradHO2_on_C  = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradHO2_on_C'));
inv_GradH2O2_on_C = importdata(fullfile(base_path, 'Condmean_on_C','Condmean_inv_gradH2O2_on_C'));

%% read laminar data files

temp_MC_lam = read_data(fullfile(base_path, 'data_MC_1.500E-01'));
temp_MA_lam = read_data(fullfile(base_path, 'data_MA_1.500E-01'));
[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(fullfile(base_path, 'config'));

X = xm;
Y = y(2:end);
Z = z(2);

%% plot flame width

Trans_T = 1185-300;

h = figure;
plot(inv_GradT_on_T_MA(2:end,1),inv_GradT_on_T_MA(2:end,2)/(delta_L_MC/(max(temp_MA_lam)-min(temp_MC_lam))),'k--','linewidth',2)%delta_L_MA
hold on
plot(inv_GradT_on_T_MC(2:end,1),inv_GradT_on_T_MC(2:end,2)/(delta_L_MC/(max(temp_MC_lam)-min(temp_MC_lam))),'k','linewidth',2)%delta_L_MC 1/max(gradient(T_MC_lam,X))
plot(inv_GradT_on_T_lam_MA(:,1),inv_GradT_on_T_lam_MA(:,2)/(delta_L_MC/(max(temp_MA_lam)-min(temp_MC_lam))),'r--','linewidth',2)%delta_L_MA
plot(inv_GradT_on_T_lam_MC(:,1),inv_GradT_on_T_lam_MC(:,2)/(delta_L_MC/(max(temp_MC_lam)-min(temp_MC_lam))),'r','linewidth',2)%delta_L_MC
plot([Trans_T Trans_T],[-2 0.7],'k-.','linewidth',1)
text(Trans_T-350,-1,'$\leftarrow$ Preheat','interpreter','latex','FontSize',14)
text(Trans_T+35,-1,'Reaction $\rightarrow$','interpreter','latex','FontSize',14)
set(gca,'FontSize',16,'linewidth',1,'fontname','Times New Roman')
xlabel('$T$ [K]','interpreter','latex','FontSize',16);
ylabel('$\langle \delta_{t}|T \rangle/ \delta_{L}$','interpreter','latex','FontSize',16);
legend('MA','MC','MA Laminar','MC Laminar','location','southeast','interpreter','latex')
legend boxoff
set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7*2]);
axis_p = get(gca,'position');
ylim = get(gca,'ylim');
axis([300 1700 -2 25])
axis square
axes('Position',[.30 .5 .35 .35]);
box on
axis square
plot(gradT_on_T_MA(:,1),gradT_on_T_MA(:,2)*(delta_L_MA/(max(temp_MA_lam)-min(temp_MA_lam))),'k--','linewidth',2)
hold on
plot(gradT_on_T_MC(:,1),gradT_on_T_MC(:,2)*(delta_L_MC/(max(temp_MC_lam)-min(temp_MC_lam))),'k','linewidth',2)
plot(gradT_on_T_lam_MA(:,1),gradT_on_T_lam_MA(:,2)*(delta_L_MA/(max(temp_MA_lam)-min(temp_MA_lam))),'r--','linewidth',2)
plot(gradT_on_T_lam_MC(:,1),gradT_on_T_lam_MC(:,2)*(delta_L_MC/(max(temp_MC_lam)-min(temp_MC_lam))),'r','linewidth',2)
set(gca,'FontSize',14,'linewidth',1,'fontname','Times New Roman')
xlabel('$T$ [K]','interpreter','latex','FontSize',14);
ylabel('$\langle \tilde{\chi}|T \rangle \delta_{L}$','interpreter','latex','FontSize',14);
axis([300 1700 0 2])
set(gcf,'color','w')
set(gcf,'units','centimeters')
pos = get(gcf,'position');
AR = pos(3)/pos(4);
set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7*2]);

exportgraphics(h, 'flame_width_conditionalmean.pdf', 'ContentType', 'vector')

%% plot flame reconstruction

T_MC = linspace(298,1700,5000);
T_lam_MC = linspace(interp1(ref_lam_MC_norm(:,2),ref_lam_MC_norm(:,1),.9999,'spline','extrap'),interp1(ref_lam_MC_norm(:,2),ref_lam_MC_norm(:,1),.001,'spline','extrap'),5000);
T_MA = linspace(298,1700,5000);
T_lam_MA = linspace(interp1(ref_lam_MA_norm(:,2),ref_lam_MA_norm(:,1),.9999,'spline','extrap'),interp1(ref_lam_MA_norm(:,2),ref_lam_MA_norm(:,1),.001,'spline','extrap'),5000);

inv_Grad_MC = interp1(inv_GradT_on_T_MC(2:end,1), inv_GradT_on_T_MC(2:end,2),T_MC,'spline','extrap');
grad_MC = interp1(gradT_on_T_MC(2:end,1), gradT_on_T_MC(2:end,2),T_MC,'spline','extrap');
inv_Grad_MA = interp1(inv_GradT_on_T_MA(2:end,1), inv_GradT_on_T_MA(2:end,2),T_MA,'spline','extrap');
grad_MA = interp1(gradT_on_T_MA(2:end,1), gradT_on_T_MA(2:end,2),T_MA,'spline','extrap');

inv_Grad_lam_MC = interp1(inv_GradT_on_T_lam_MC(2:44,1), inv_GradT_on_T_lam_MC(2:44,2),T_lam_MC,'spline','extrap');
grad_lam_MC = interp1(gradT_on_T_lam_MC(2:44,1),gradT_on_T_lam_MC(2:44,2),T_lam_MC,'spline','extrap');
inv_Grad_lam_MA = interp1(inv_GradT_on_T_lam_MA(2:44,1),inv_GradT_on_T_lam_MA(2:44,2),T_lam_MA,'spline','extrap');
grad_lam_MA = interp1(gradT_on_T_lam_MA(2:44,1),gradT_on_T_lam_MA(2:44,2),T_lam_MA,'spline','extrap');

temp_MC = cumtrapz(T_MC, inv_Grad_MC);
struc_MC = temp_MC - interp1(T_MC,temp_MC,interp1(grad_MC,T_MC,max(grad_MC),'spline','extrap'),'spline','extrap');

temp_MC = cumtrapz(T_lam_MC, inv_Grad_lam_MC);
struc_lam_MC = temp_MC - interp1(T_lam_MC,temp_MC,interp1(grad_lam_MC,T_lam_MC,max(grad_lam_MC),'spline','extrap'),'spline','extrap');

temp_MA = cumtrapz(T_MA, inv_Grad_MA);
struc_MA = temp_MA - interp1(T_MA, temp_MA, interp1(grad_MA,T_MA,max(grad_MA),'spline','extrap'),'spline','extrap');
temp_MA = cumtrapz(T_lam_MA, inv_Grad_lam_MA);
struc_lam_MA = temp_MC - interp1(T_lam_MA, temp_MA, interp1(grad_lam_MA, T_lam_MA,max(grad_lam_MA),'spline','extrap'),'spline','extrap');

h = figure;
plot(struc_MA/delta_L_MA,T_MA,'k--','linewidth',2), hold on
plot(struc_MC/delta_L_MC,T_MC,'k','linewidth',2)
plot([-1 0.1],[Trans_T Trans_T],'k-.','linewidth',1)
set(gca,'FontSize',14,'linewidth',1);%'fontname','Times New Roman')
xlabel('$\langle n|T \rangle/ \delta_{L}$','interpreter','latex','FontSize',16);
ylabel('$T$ [K]','interpreter','latex','FontSize',16);
axis([-1 12 0 1700])
axis square
axes('Position',[.38 .2 .5 .5]);
box on
axis square
plot(struc_MA/delta_L_MA,T_MA,'k--','linewidth',2), hold on
plot(struc_MC/delta_L_MC,T_MC,'k','linewidth',2)
plot(struc_lam_MA./delta_L_MA,T_lam_MA,'r--','linewidth',2)
plot(struc_lam_MC./delta_L_MC,T_lam_MC,'r','linewidth',2)
plot([-1 0.1],[Trans_T Trans_T],'k-.','linewidth',1)
text(-.9,Trans_T-60,'$\downarrow$ Preheat','interpreter','latex','FontSize',12)
text(-.9,Trans_T+60,'$\uparrow$ Reaction','interpreter','latex','FontSize',12)
set(gca,'FontSize',14,'linewidth',1);%,'fontname','Times New Roman')
axis([-1 2 0 1400])
set(gcf,'color','w')
set(gcf,'units','centimeters')
pos = get(gcf,'position');
AR = pos(3)/pos(4);
set(gcf,'units','centimeters','position',[0,0,6.7*2,6.7*2]);

legend('MA','MC','MA Laminar','MC Laminar','location','SE','interpreter','latex')
legend boxoff

exportgraphics(h, 'reconstructed_flame.pdf', 'ContentType', 'vector')

%% functions

function T = read_data(filename)
% reads temperature only from NGA data file

fid     = fopen(filename, 'r');

dims    = fread(fid,4,'integer*4','ieee-le'); % read only the first 4 elements from the file
% Dimensions 
nx      = dims(1);
ny      = dims(2);
nz      = dims(3);
nvar    = dims(4);
% dt      = dims(5);
nsize   = nx*ny*nz;

% The following lines are necessary to move the file pointer to the correct
% location befory you start reading in the variables
dt      = fread(fid,1,'real*8','ieee-le'); % skipping over 'dt'
Time    = fread(fid,1,'real*8','ieee-le'); % skipping over 'time'
nt      = Time/dt;

varnames = [];
for var=1:nvar
    varnames = [varnames; fread(fid,8,'*char','ieee-le')]; % skipping over the stored names
end

% There are usually 5 variables (or nvar number of variables) stored in the data
% file: U, V, W, P, ZMIX. Each fread call makes the file pointer shift to the end 
% of the number of values read. Hence each subsequent fread call will start after 
% where the previous fread call stopped. Add on extra variables if necessary.

dummy       = fread(fid,nsize,'real*8','ieee-le');    % will produce a column vector
dummy           = reshape(dummy,nx,ny,nz);            % now turning the column vector into a 3D matrix
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy         = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy        = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy          = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy          = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy          = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy          = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy         = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy         = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy        = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
T           = reshape(dummy,nx,ny,nz);
dummy       = fread(fid,nsize,'real*8','ieee-le');
dummy        = reshape(dummy,nx,ny,nz);

fclose(fid);

end