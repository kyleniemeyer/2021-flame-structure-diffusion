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
[X,Y,Z] = meshgrid(xm,ym,zm);

%% Read data

for n = 1:numfiles
    U = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 1);
    V = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 2);
    W = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 3);

    % calculate gradient of velocity field
    [tmpXy,tmpXx,tmpXz] = gradient(U,dx,dy,dz);
    [tmpYy,tmpYx,tmpYz] = gradient(V,dx,dy,dz);
    [tmpZy,tmpZx,tmpZz] = gradient(W,dx,dy,dz);
    
    % vorticity
    %[omegaY,omegaX,omegaZ] = curl(X,Y,Z,U,V,W);
    % calculate vorticity directly
    omegaX = tmpZy - tmpYz;
    omegaY = tmpXz - tmpZx;
    omegaZ = tmpYx - tmpXy;
    
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
    density = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 5);
    density_inv = 1.0 ./ density;
    
    % this is div(u)
    %TMP1 = divergence(X,Y,Z,U,V,W);
    div_velocity = tmpXx + tmpYy + tmpZz;
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                % this is local grad(u)
                A = [tmpXx(i,j,k), tmpXy(i,j,k), tmpXz(i,j,k);
                     tmpYx(i,j,k), tmpYy(i,j,k), tmpYz(i,j,k);
                     tmpZx(i,j,k), tmpZy(i,j,k), tmpZz(i,j,k)];
                % this is local tau, the viscous stress tensor
                B = viscosity(i,j,k).*((A+A') - (2./3.).*div_velocity(i,j,k)*eye(3));
                
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
    % tmpXx, tmpXy, ... now has the viscous stress tensor field (tau)

    % divergence of tau
    % avoid using the divergence function
    %tmpX = divergence(Y,X,Z,tmpXy,tmpXx,tmpXz);
    %tmpY = divergence(Y,X,Z,tmpYy,tmpYx,tmpYz);
    %tmpZ = divergence(Y,X,Z,tmpZy,tmpZx,tmpZz);
    [~,tmp_x,~] = gradient(tmpXx,dx,dy,dz);
    [tmp_y,~,~] = gradient(tmpXy,dx,dy,dz);
    [~,~,tmp_z] = gradient(tmpXz,dx,dy,dz);
    tmpX = tmp_x + tmp_y + tmp_z;

    [~,tmp_x,~] = gradient(tmpYx,dx,dy,dz);
    [tmp_y,~,~] = gradient(tmpYy,dx,dy,dz);
    [~,~,tmp_z] = gradient(tmpYz,dx,dy,dz);
    tmpY = tmp_x + tmp_y + tmp_z;

    [~,tmp_x,~] = gradient(tmpZx,dx,dy,dz);
    [tmp_y,~,~] = gradient(tmpZy,dx,dy,dz);
    [~,~,tmp_z] = gradient(tmpZz,dx,dy,dz);
    tmpZ = tmp_x + tmp_y + tmp_z;
    clear tmp_x tmp_y tmp_z
    clear tmpXx tmpXy tmpXz tmpYx tmpYy tmpYz tmpZx tmpZy tmpZz

    % curl (1/rho * div(tau))
    %[tmpY,tmpX,tmpZ] = curl(X,Y,Z, density_inv.*tmpX, density_inv.*tmpY, density_inv.*tmpZ);
    [tmpXy,tmpXx,tmpXz] = gradient(density_inv.*tmpX, dx,dy,dz);
    [tmpYy,tmpYx,tmpYz] = gradient(density_inv.*tmpY, dx,dy,dz);
    [tmpZy,tmpZx,tmpZz] = gradient(density_inv.*tmpZ, dx,dy,dz);
    tmpX = tmpZy - tmpYz;
    tmpY = tmpXz - tmpZx;
    tmpZ = tmpYx - tmpXy;

    % omega dot curl (1/rho * div(tau))
    TMP2 = omegaX.*tmpX + omegaY.*tmpY + omegaZ.*tmpZ;
    viscous_effects = mean(mean(TMP2,3),2);

    clear viscosity TMP2 tmpX tmpY tmpZ
    clear tmpXx tmpXy tmpXz tmpYx tmpYy tmpYz tmpZx tmpZy tmpZz
    % save div u for next step
    
    % enstrophy
    TMP1 = omegaX.*omegaX + omegaY.*omegaY + omegaZ.*omegaZ;
    enstrophy = mean(mean(TMP1,3),2);

    % Dilitation
    % $ -\omega^2 (\nabla \cdot u) $
    TMP2 = TMP1 .* div_velocity;
    dilatation = -1.0 * mean(mean(TMP2,3),2);

    clear TMP1 TMP2 div_velocity
    
    % Baroclinic torque
    % $ \frac{\omega}{\rho^2} \cdot ( \nabla \rho \times \nabla P ) $

    pressure = NGAdatareader_large(fullfile(data_path, filenames(n,:)), 4);

    [tmpy,tmpx,tmpz] = gradient(density, dx, dy, dz);
    [tmpY,tmpX,tmpZ] = gradient(pressure, dx, dy, dz);
    TMP1 = tmpy.*tmpZ - tmpz.*tmpY; % x component of $(\nabla \rho \times \nabla P)$
    TMP2 = tmpz.*tmpX - tmpx.*tmpZ; % y component of $(\nabla \rho \times \nabla P)$
    TMP3 = tmpx.*tmpY - tmpy.*tmpX; % z component of $(\nabla \rho \times \nabla P)$
    TMP1 = (omegaX.*TMP1 + omegaY.*TMP2 + omegaZ.*TMP3) .* (density_inv.^2);
    baroclinic = mean(mean(TMP1,3),2);

    clear density pressure density_inv TMP1 TMP2 TMP3 tmpX tmpY tmpZ tmpx tmpy tmpz
    
    % Forcing
    % forcing term calculation based on Bobbit 2016
    % $A \omega^2 = \omega^2 / (2 \tau_0)$
    TMP2 = (omegaX.*omegaX + omegaY.*omegaY + omegaZ.*omegaZ) .* 973.05;
    forcing = mean(mean(TMP2,3),2);

    clear omegaX omegaY omegaZ TMP2
    
    savename = join(['enstrophy_', num2str(n), '_', transport_model, '.mat']);
    save(savename, 'enstrophy', 'stretch', 'dilatation', 'baroclinic', 'viscous_effects', 'forcing')

    disp(['Done with ', num2str(n)])
    
end

toc
