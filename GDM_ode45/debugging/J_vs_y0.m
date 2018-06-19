% This code calculates layser averaged particle density based on a 1D drift-diffusion model. 
% For a HOLE-ONLY device
% length units in m, energy units in eV
% Mobility dependence is set with a GDM framework. Parameters are adapted from van Mensfoort 2010 paper.
% Incorporate both type I and II solutions.

clear all;
close all;
clc

tic

%% ********************** Parameter settings ********************
phi1 = 0.4;     % Barrier at left side
phi2 = 0.6;     % Barrier at right side: In this code, set them equal!

epir=3.0;
Nt = 2e26;   % unit in m^(-3)
Nx=100;          % Unit in a (lattice constant)
a=(1/Nt)^(1/3);       % lattice constant, unit: m
L = Nx*a;
epi0=8.854187817e-12;
e=1.602176565e-19;     % unit in Coulomb
T = 189;                % temperature in K
kbT=1.38e-23*T/e;          % unit: eV
esig = 0.14;        % in eV
esig_interface=esig;       % Disorder strength, unit: eV

%*********************** End Parameter settings **************************

%% =========================== GDM parameters ===========================
fun_h_1=@(x) Nt/(sqrt(2*pi)*esig_interface) * exp(-x.^2/(2*esig_interface^2)) .* (1./(1+exp((x+phi1)/kbT)));
n1 = integral(fun_h_1,-inf,inf);
% n1 = 2.13e26;   %match Fig6
fun_h_2=@(x) Nt/(sqrt(2*pi)*esig_interface) * exp(-x.^2/(2*esig_interface^2)) .* (1./(1+exp((x+phi2)/kbT)));
n2 = integral(fun_h_2,-inf,inf);
% n2 = 2.13e26;   %match Fig6 of Mensfoort 2008
% n2 = n1*1e-4;
mu0_star=22e-6;        % Units in m^2/(Vs)
C = 0.42;
mu0 = mu0_star * exp(-C*(esig/kbT)^2);
% mu0 = 1e-10;
delta = 2*(log((esig/kbT)^2-esig/kbT)-log(log(4)))/(esig/kbT)^2;

% ================================= End GDM parameters ==============================

warning_flag = 0;
%% ========================================= Enters main loop =============================================

% J = logspace(-1,4,20)';
J = 1;
tolerence = 1e-3;
f_upper = 12800;    % tspan parameter in ode45 (see Matlab manual)
f_lower = -12800;   % has to be large so that the cutoff f limit does not exceed the range [f_lower,f_upper]

%% ================================== ENTER TYPE II SOLUTION SOLVER: ==================================
ymin = linspace(5e-5,1e-6,10);
interp_para = zeros(1,4);
for p=1:length(ymin)
    if mod(p,4)==0
        disp(p);
    end
    % dimensionless parameters defined in Appendix B
    gamma1 = e*L^2/(epi0*epir*kbT)*n1;
    gamma2 = e*L^2/(epi0*epir*kbT)*n2; 
    i = 1/kbT^2*(Nx*a)^3/(epi0*epir*mu0)*J;
    y1 = gamma1/i^(2/3);
    y2 = gamma2/i^(2/3);

    coefficients = [delta,kbT,mu0,L,Nt,epi0,epir,e,a,esig,J,Inf,ymin(p),y1,y2,f_lower,f_upper];
    [i1,fmin_1,warning_flag,vpa_flag] = calc_i_type2_ode45_GDM_original_df_dy(coefficients,1);
    if warning_flag~=0
        fprintf('Serious Warning: f range [f_lower, f_upper] not large enough, or y2 too big!\n');
        fprintf('ymin_initial = %f\n',ymin);
    end
    J1 = epi0*epir*mu0/L^3*(kbT)^2 * i1;
    interp_para = [interp_para;ymin(p), fmin_1, i1, J1];

end

% ================================== END TYPE II SOLUTION SOLVER ==================================
interp_para(1,:)=[];

toc;