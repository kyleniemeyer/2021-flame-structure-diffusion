%%% Calculates terms in enstrophy budget

%% Setup

clear, clc
close all

tic

transport_model = 'MA';
%transport_model = 'MC';

data_path = '/nfs/attic/niemeyek/ProCI_Vorticity_paper';
configfile = fullfile(data_path, 'config');

filenames_MA = ["data.init"; ...
    "data_MA_Vort_2.020E-02"; ...
    "data_MA_Vort_2.040E-02"; ...
    "data_MA_Vort_2.044E-02"; ...
    "data_MA_Vort_2.060E-02"; ...
    "data_MA_Vort_2.080E-02"; ...
    "data_MA_Vort_2.100E-02"; ...
    "data_MA_Vort_2.120E-02"; ...
    "data_MA_Vort_2.140E-02"; ...
    "data_MA_Vort_2.160E-02"; ...
    "data_MA_Vort_2.180E-02"; ...
    "data_MA_Vort_2.200E-02"; ...
    "data_MA_Vort_2.220E-02"; ...
    "data_MA_Vort_2.240E-02"; ...
    "data_MA_Vort_2.260E-02"; ...
    "data_MA_Vort_2.280E-02"; ...
    "data_MA_Vort_2.300E-02"; ...
    "data_MA_Vort_2.320E-02"; ...
    "data_MA_Vort_2.340E-02"; ...
    "data_MA_Vort_2.360E-02"; ...
    "data_MA_Vort_2.380E-02"; ...
    "data_MA_Vort_2.400E-02"; ...
    "data_MA_Vort_2.420E-02"];

filenames_MC = ["data.init"; ...
    "data_MC_Vort_2.020E-02"; ...
    "data_MC_Vort_2.040E-02"; ...
    "data_MC_Vort_2.044E-02"; ...
    "data_MC_Vort_2.060E-02"; ...
    "data_MC_Vort_2.080E-02"; ...
    "data_MC_Vort_2.100E-02"; ...
    "data_MC_Vort_2.120E-02"; ...
    "data_MC_Vort_2.140E-02"; ...
    "data_MC_Vort_2.160E-02"; ...
    "data_MC_Vort_2.180E-02"; ...
    "data_MC_Vort_2.200E-02"; ...
    "data_MC_Vort_2.220E-02"; ...
    "data_MC_Vort_2.240E-02"; ...
    "data_MC_Vort_2.256E-02"; ...
    "data_MC_Vort_2.260E-02"; ...
    "data_MC_Vort_2.280E-02"; ...
    "data_MC_Vort_2.300E-02"; ...
    "data_MC_Vort_2.320E-02"; ...
    "data_MC_Vort_2.340E-02"; ...
    "data_MC_Vort_2.360E-02"; ...
    "data_MC_Vort_2.380E-02"; ...
    "data_MC_Vort_2.400E-02"; ...
    "data_MC_Vort_2.420E-02"; ...
    "data_MC_Vort_2.440E-02"];

filenames = filenames_MA;

% these files are not actually needed
optdata_filenames = ["optdata.init"];
numfiles = length(filenames);

[x,y,z,xm,ym,zm,nx,ny,nz,dx,dy,dz] = NGA_grid_reader(configfile);
[X,Y,Z] = meshgrid(ym,xm,zm);

%% Read data

for n = 1:numfiles
    U = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 1);
    V = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 2);
    W = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 3);
    
    % calculate enstrophy: 0.5 * omega^2
    [omegaX,omegaY,omegaZ] = curl(X,Y,Z,U,V,W);
    TMP1 = 0.5*(omegaX.*omegaX + omegaY.*omegaY + omegaZ.*omegaZ);
    LHS = mean(mean(TMP1,3),2);

    clear TMP1
    
     % Viscous effects part 1
    [tmpXx,tmpXy,tmpXz] = gradient(U,dx,dy,dz);
    [tmpYx,tmpYy,tmpYz] = gradient(V,dx,dy,dz);
    [tmpZx,tmpZy,tmpZz] = gradient(W,dx,dy,dz);
    
    % Stretch term
    tmpX = 0.5*((tmpXx+tmpXx).*omegaX + (tmpXy+tmpYx).*omegaY + (tmpXz+tmpZx).*omegaZ);
    tmpY = 0.5*((tmpYx+tmpXy).*omegaX + (tmpYy+tmpYy).*omegaY + (tmpYz+tmpZy).*omegaZ);
    tmpZ = 0.5*((tmpZx+tmpXz).*omegaX + (tmpZy+tmpYz).*omegaY + (tmpZz+tmpZz).*omegaZ);
    TMP1 = omegaX.*tmpX + omegaY.*tmpY + omegaZ.*tmpZ;
    stretch = mean(mean(TMP1,3),2);

    clear TMP1 tmpX tmpY tmpZ
    
    % Viscous effects
    % $\omega \times \left( \frac{1}{\rho} \nabla \cdot \tau \right)$
    % where $\tau = 2 \mu \left( S - \frac{1}{3} (\nabla \cdot u) I \right)$
    % and $S = \frac{1}{2} \left( \nabla u + \nabla u^{\top} \right)$

    viscosity = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 17);
    rho = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 5);
    rho_inv = 1.0 ./ rho;
    
    % this is div(u)
    TMP1 = divergence(X,Y,Z,U,V,W);
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                % this is local grad(u)
                A = [tmpXx(i,j,k), tmpXy(i,j,k), tmpXz(i,j,k);
                     tmpYx(i,j,k), tmpYy(i,j,k), tmpYz(i,j,k);
                     tmpZx(i,j,k), tmpZy(i,j,k), tmpZz(i,j,k)];
                % this is local tau, the viscous stress tensor
                B = viscosity(i,j,k).*((A+A') - (2./3.).*TMP1(i,j,k)*eye(3));
                
                % this overwrites the local grad(u) entries, which are not needed again
                tmpXx(i,j,k) = B(1,1);
                tmpXy(i,j,k) = B(1,2);
                tmpXz(i,j,k) = B(1,3);
                tmpYx(i,j,k) = B(2,1);
                tmpYy(i,j,k) = B(2,2);
                tmpYz(i,j,k) = B(2,3);
                tmpZx(i,j,k) = B(3,1);
                tmpZy(i,j,k) = B(3,2);
                tmpZz(i,j,k) = B(3,3);
            end
        end
    end
    % divergence of tau
    tmpX = divergence(X,Y,Z,tmpXx,tmpXy,tmpXz);
    tmpY = divergence(X,Y,Z,tmpYx,tmpYy,tmpYz);
    tmpZ = divergence(X,Y,Z,tmpZx,tmpZy,tmpZz);
    % curl (1/rho * div(tau))
    [tmpX,tmpY,tmpZ] = curl(X,Y,Z, rho_inv.*tmpX, rho_inv.*tmpY, rho_inv.*tmpZ);
    % omega dot curl (1/rho * div(tau))
    TMP2 = omegaX.*tmpX + omegaY.*tmpY + omegaZ.*tmpZ;
    viscous_effects = mean(mean(TMP2,3),2);

    clear tmpXx tmpXy tmpXz tmpYx tmpYy tmpYz tmpZx tmpZy tmpZz viscosity TMP2 tmpX tmpY tmpZ
    % save TMP1 (div u) for next step
    
    % Dilitation
    % $ -\omega^2 (\nabla \cdot u) $
    %TMP1 = divergence(X,Y,Z,U,V,W); % <- this is done in previous step
    tmpX = omegaX .* omegaX;
    tmpY = omegaY .* omegaY;
    tmpZ = omegaZ .* omegaZ;
    TMP2 = tmpX.*TMP1 + tmpY.*TMP1 + tmpZ.*TMP1;
    dilatation = -1.0 * mean(mean(TMP2,3),2);

    clear TMP2 TMP1 tmpX tmpY tmpZ
    
    % Baroclinic torque
    % $ \frac{\omega}{\rho^2} \cdot ( \nabla \rho \times \nabla P ) $

    % old way
    % [tmpx,tmpy,tmpz] = gradient(RHO, dx, dy, dz);
    % [tmpX,tmpY,tmpZ] = gradient(P, dx, dy, dz);
    % for i = 1:nx
    %     for j = 1:ny
    %         for k = 1:nz
    %             a = [tmpx(i,j,k) tmpy(i,j,k) tmpz(i,j,k)];
    %             b = [tmpX(i,j,k) tmpY(i,j,k) tmpZ(i,j,k)];
    %             TMP2(i,j,k) = sum(RHO_inv(i,j,k)^2.*cross(a,b));
    %         end
    %     end
    % end
    % TMP1 = omegaX.*TMP2 + omegaY.*TMP2 + omegaZ.*TMP2;

    pressure = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 4);

    [tmpx,tmpy,tmpz] = gradient(rho, dx, dy, dz);
    [tmpX,tmpY,tmpZ] = gradient(pressure, dx, dy, dz);
    TMP1 = tmpy.*tmpZ - tmpz.*tmpY; % x component of $(\nabla \rho \times \nabla P)$
    TMP2 = tmpz.*tmpX - tmpx.*tmpZ; % y component of $(\nabla \rho \times \nabla P)$
    TMP3 = tmpx.*tmpY - tmpy.*tmpX; % z component of $(\nabla \rho \times \nabla P)$
    TMP1 = (omegaX.*TMP1 + omegaY.*TMP2 + omegaZ.*TMP3) .* (rho_inv.^2);
    baroclinic = mean(mean(TMP1,3),2);

    clear rho pressure rho_inv TMP1 TMP2 TMP3 tmpX tmpY tmpZ tmpx tmpy tmpz
    
    % Forcing

    % old (incorrect) forcing term
    %tmpx = RHO_inv.*srcU;
    %tmpy = RHO_inv.*srcV;
    %tmpz = RHO_inv.*srcW;
    %[tmpX,tmpY,tmpZ] = curl(X,Y,Z,tmpx,tmpy,tmpz);
    %TMP1 = omegaX.*tmpX+omegaY.*tmpY+omegaZ.*tmpZ;
    %Forcing = mean(mean(TMP1,3),2);

    % forcing term calculation based on Bobbit 2016
    % $A \omega^2 = \omega^2 / (2 \tau_0)$
    TMP2 = (omegaX.*omegaX + omegaY.*omegaY + omegaZ.*omegaZ).*(973.05);
    forcing = mean(mean(TMP2,3),2);
    
    savename = join(['enstrophy_', num2str(n), '_', transport_model, '.mat']);
    save(savename, 'LHS', 'stretch', 'dilatation', 'baroclinic', 'viscous_effects', 'forcing')

    clear omegaX omegaY omegaZ TMP2
    
end

toc
