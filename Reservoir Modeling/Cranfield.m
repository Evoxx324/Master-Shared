clear all;
close all;
mymrst = '/Users/moortgat/Box Sync/Teaching/5751_2018/mrst-2017a';
addpath(genpath(mymrst));

load Cranfield.mat;

format long g;
Grid.Nx=64; Grid.Lx=195.07031; Grid.hx = Grid.Lx/Grid.Nx;            % Dimension in x-direction
Grid.Ny=51; Grid.Ly=155.45312; Grid.hy = Grid.Ly/Grid.Ny;            % Dimension in y-direction
Grid.Nz=79;  Grid.Lz=25;  Grid.hz = Grid.Lz/Grid.Nz;            % Dimension in z-direction
N=Grid.Nx*Grid.Ny*Grid.Nz;                                    % Total number of grid blocks

Grid.V   = Grid.hx*Grid.hy*Grid.hz;                           % Cell volumes
Grid.K(1,1:Grid.Nx,1:Grid.Ny,1:Grid.Nz) = reshape(perm(1:N),Grid.Nx,Grid.Ny,Grid.Nz); 
Grid.K(2,1:Grid.Nx,1:Grid.Ny,1:Grid.Nz) = reshape(perm(1:N),Grid.Nx,Grid.Ny,Grid.Nz); 
Grid.K(3,1:Grid.Nx,1:Grid.Ny,1:Grid.Nz) = reshape(perm(1:N),Grid.Nx,Grid.Ny,Grid.Nz); 
Grid.por = reshape(poro(1:N),Grid.Nx,Grid.Ny,Grid.Nz);                 % Unit porosity
Grid.pv  = Grid.V(:).*Grid.por(:);                            % pore volume=cell volume*porosity
PVtot = sum(Grid.pv);
%
Q=zeros(N,1);
Q([1 N])=[0.47 -0.47];                     % Production/injection
%
Fluid.vw=0.022;  Fluid.vo=0.22;                         % Viscosities
Fluid.swc=0.; Fluid.sor=0.4;                      % Irreducible saturations
Fluid.kr0w = 0.8; 
Fluid.kr0o = 1;
Fluid.nw = 2.6;
Fluid.no = 4.2;
BHP = 300;
%
S=Fluid.swc*ones(N,1);                             % Initial saturation

nt = 30;                                           % Time steps
tmax = 0.47*PVtot;                                  % nr of pressure time-steps 
cflx = 1;                                          % how many implicit saturation time-steps per pressure time-step
timemethod = 1;                                    % explicit=1, implicit=2; explicit Runge-Kutta 4
presmethod = 1;                                    % TPFA=1, MFE = 2
plotstyle = 4;
tic
Simulator( Grid, Fluid, Q, tmax, nt,timemethod, cflx , presmethod,plotstyle,BHP)
toc