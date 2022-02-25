function [Pmat,ux, uy, Ux,Uy,Tmat,b] = LCL(dx,dy,dz,k,mu,Q_in,p_out)
% LOCAL CUBIC LAW: Reynolds approximation
% [Pmat,ux, uy, Ux,Uy,Tmat,b] = LCL(dx,dy,dz,k,mu,Q_in,p_out)
%
% This funtion computes the Pressure Field and evaluates the Fluxes in the
% 2D-domain specified by the user. The numerical solution is computed using
% a FVM with piece-wise constant mobilities and the two-point flux
% approximation.
% 
% Boundary conditions: Bottom = Constant pressure (p_in)
%                      Top    = Constant pressure (p_out)
%
% Note that this function solvese the boundary problem assuming constant 
% pressure Pin, and Pout, and then rescales the solution depending on the 
% user BCs.
%
%
% Inputs:
% -------------------------------------------------------------------------
% dx    = Grid spacing in x-direction                [Scalar]  [L]
% dy    = Grid spacing in y-direction                [Scalar]  [L]
% dz    = Aperture field (along z-direction)         [Matrix]  [L]
% k     = 2D Permeability field                      [Matrix]  [L^2]
% mu    = Dynamic viscosity pure liquid              [Scalar]  [M/(L·T)]
% Q_in  = Prescribed flow rate @ bottom              [Scalar]  [L^3/T)]
% p_out = Prescribed pressure @ top                  [Scalar]  [M/(L·T^2)]
%
% Outputs:
% -------------------------------------------------------------------------
% Pmat  = Pressure field                             [Matrix]  [M/(L·T^2)]
% ux    = Fluxes @ center of element in x-direction  [Matrix]  [Units]
% uy    = Fluxes @ center of element in y-direction  [Matrix]  [Units]
% Ux    = Fluxes across element in x-direction       [Matrix]  [Units]
% Uy    = Fluxes across element in y-direction       [Matrix]  [Units]
% Tmax  = System transmissivity                      [Matrix]
% b     = Boundary conditions                        [Vector]
%
% Author:
% -------------------------------------------------------------------------
% Rafael Villamor Lora
% April 8, 2020 [Last modified (04/09/2020)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER DEFINITIONS                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = size(k,2);           % Number of elements in x-direction
Ny = size(k,1);           % Number of elements in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSAMBLE TRANSMISSIVITY MATRIX                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix that contains the mobilities
lambda = k./mu;
% Matrices that contain the local (vertical) apertures. Use harmonic avg.
dz_x = 2 ./ [1./dz(:,1:end-1) + 1./dz(:,2:end)]';
dz_y = 2 ./ [1./dz(1:end-1,:) + 1./dz(2:end,:)]';

% TRANSMISIBILITY MATRICES (Tx_l , Tx_r, Ty_b, Ty_t,T)
% Matrix that contains x-transmisibilities (Assume no flow boundaries)
Tx   = [zeros(1,Ny); 2*dy.*dz_x/dx./ [1./lambda(:,1:end-1) + 1./lambda(:,2:end)]'; zeros(1,Ny)];
% Matrix that contains y-transmisibilities (Assume constant pressure boundaries)
Ty   = [[2*dx/dy.*dz(1,:).*lambda(1,:)]'  2*dx/dy.*dz_y./[1./lambda(1:end-1,:)+...
      1./lambda(2:end,:)]' [2*dx/dy.*dz(end,:).*lambda(end,:)]'];
% x-Transmibilibities left (Tx_l)
Tx_l = reshape(Tx(1:end-1,:)',Nx*Ny,1);
T3   = [Tx_l(Ny+1:end);zeros(Ny,1)];
% x-Transmibilibities right (Tx_r)
Tx_r = reshape(Tx(2:end,:)',Nx*Ny,1);
T4   = [zeros(Ny,1);Tx_r(1:end-Ny)];
% y-Transmibilibities bottom (Ty_b)
Ty_b = reshape(Ty(:,1:end-1)',Nx*Ny,1);
T1   = reshape([Ty(:,2:end-1) zeros(Nx,1)]',Nx*Ny,1); % Disconnect top and bottom elements
% y-Transmibilibities top (Ty_t)
Ty_t = reshape(Ty(:,2:end)',Nx*Ny,1);
T2   = reshape([zeros(Nx,1) Ty(:,2:end-1)]',Nx*Ny,1); % Disconnect top and bottom elements
% System Transmisibility Matrix (T)
T0   = Tx_l + Tx_r + Ty_b + Ty_t;
Tmat = spdiags([-T3,-T1, T0,-T2,-T4],...
    [-Ny -1 0 1 Ny],Nx*Ny,Nx*Ny);
disp('Transmissivity matrix assembled')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITIONS                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITIONS
p_in = 1.5*p_out;

% RIGHT HAND SIDE VECTOR (b)
b = zeros(Nx*Ny,1);
% Find the inlet cell numbers
inlet     = Ny * [0:1:Nx-1] + 1;
% Find the outlet cell numbers
outlet    = Ny * [1:1:Nx];
% Correct for Boundary Conditions
% Inlet: Prescribed pressure
b(inlet)  = Ty_b(inlet)  *p_in;
% Outlet: Prescribed pressure
b(outlet) = Ty_t(outlet) *p_out;
disp('Boundary Condition vector assembled')
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE THE PRESSURE EQUATION AND EVALUATE FLUXES                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% SOLVE STEADY STATE PROBLEM (P)
p = Tmat\b;
Pmat = [reshape(p,Ny,Nx)]';
disp('Pressure field solved')
  
% EVALUATE INTEGRATED FLUXES 
% x-fluxes
Ux = Tx.*([zeros(1,Ny);Pmat]-[Pmat;zeros(1,Ny)]);
% y-fluxes
Uy = Ty.*([p_in*ones(Nx,1) Pmat]-[Pmat p_out*ones(Nx,1)]);
disp('Velocity field solved')

% EVALUATE MEAN FLUXES AT THE CENTER OF THE CELL
ux = (Ux(1:end-1,:) + Ux(2:end,:)) / 2;
uy = (Uy(:,1:end-1) + Uy(:,2:end)) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           RE-SCALE THE PROBLEM                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So far the linear problem has been solved assuming constant pressure at
% the inlet/outlet boundaries. However, the original problem had other BCs:
% Inlet:  constant flux
% Outlet: constant pressure
%
% Since it is a linear system, one can rescale the problem just by
% multiplying both Pmat and ux, uy by the same factor

% CROSS-SECTIONAL FLUX ALONG FLOW DIRECTION
Qy      = sum(Uy,1);
Q_error = (max(Qy) - min(Qy)) ./ min(Qy) *100;
disp(['Relative error in mass conservation: ',num2str(Q_error,'%0.5f'), '%'])

% RESCALE
ScalingFactor = Q_in ./ mean(Qy); % Scaling factor
Pmat = Pmat * ScalingFactor + (1 - ScalingFactor) * p_out;
Ux   = Ux * ScalingFactor;
Uy   = Uy * ScalingFactor;
ux   = ux * ScalingFactor;
uy   = uy * ScalingFactor;
end

