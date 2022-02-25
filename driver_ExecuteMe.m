% Rafael Villamor Lora
% February 15, 2022
% SIMULATE LINEAR FLOW FROM APERTURE MAPS
%
% This script shows, step-by-step, how to simulate 2D steady laminar flow 
% of a Newtonian fluid with constant density using the Local Cubic Law
% (LCL) approximation.

clear all, close all, clc

% LOAD DATA: APERTURE FIELDS
%{
apertureField : 2D matrix of the local fracture apertures
dx            : grid spacing (pixel size)
units         : length units of apertureField and dx
%}
load("apertureField.mat")

%% FLOW SIMULATIONS
% Fluid properties
rho          = 997;                   % Density (water @ 25ºC)                [kg/m3]
mu           = 0.000889;              % Dynamic viscosity (water @ 25ºC)      [Pa·s]
% Boundary conditions
Q_in         = 0.01 * 1e-6;           % Inlet flow                            [m3/s]
p_out        = 50   * 1e3;            % Outlet pressure                       [Pa]
% Domain geometry
dx           = dx * 1e-3;             % Grid spacing in x-direction           [m]
dy           = dx;                    % Grid spacing in y-direction           [m]
dz           = apertureField * 1e-3;  % Aperture field                        [m]
bmin         = 2 * 1e-6;              % Effective aperture at contact regions [m]
dz(dz< bmin) = bmin;
% Convert the local apertures to local permeabilities using the LCL
% approximation
k            = (dz .^ 2) / 12;        % From Local Cubic Law model
% Determine the pressure and flow fields using the Local Cubic Law
% approximation
[Pmat, ux, uy] = LCL(dx,dy,dz,k,mu,Q_in,p_out);

%% PLOT RESULTS
figure; set(gcf,'color','w')
% Aperture field
ax(1) = subplot(3,1,1); 
Z = dz;
h = pcolor(Z'*1e6); set(h, 'EdgeColor', 'none'); axis equal, axis tight, axis off
hcb = colorbar; title(hcb, '(\mum)'), caxis([0 100]), title('Aperture map'), colormap(ax(1), parula);
% Pressure drop
ax(2) = subplot(3,1,2); 
Pfield = (Pmat - mean(Pmat(:,end))); Pfield = flipud(Pfield ./ mean(Pfield(:,1)));
h = pcolor(Pfield); set(h, 'EdgeColor', 'none'); axis equal, axis tight, axis off
hcb = colorbar; title(hcb, '(-)'), title('Normalized pressure, P / P_{in}'), colormap(ax(2), summer);
% Fluxes
ax(3) = subplot(3,1,3); 
u    = (ux.^2 + uy.^2).^0.5;   % Flux/pixel
urel = u ./ (Q_in/size(dz,2)); % Flux relative to uniform form
h = pcolor(urel); set(h, 'EdgeColor', 'none'); axis equal, axis tight, axis off
hcb = colorbar; title(hcb, '(-)'), caxis([0 2.5]), title('Flux / Flux_{uniform}'), colormap(ax(3), hot);
startx = 1:round(size(dz,2)/20):size(dz,2);
starty = ones(length(startx),1);
a = streamline(uy,ux,starty,startx); set(a, 'Color', 'white'); % streamlines