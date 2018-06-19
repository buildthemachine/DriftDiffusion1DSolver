% This code tries to calculate numerical solutions to equation (B14), and
% compare it with vpasolve results.

% CONCLUSION: VPASOLVE results are reliable

clear all;
close all;

tic

%% ********************** Parameter settings ********************
phi1 = 1.2;     % Barrier at left side
phi2 = 1.2;     % Barrier at right side: In this code, set them equal!

epir=3.0;
esig_interface=0.3;       % Disorder strength, unit: eV
Nt = 2e26;   % unit in m^(-3)
Nx=100;          % Unit in a (lattice constant)
a=1.0e-9;       % lattice constant, unit: m
L = Nx*a;
epi0=8.854187817e-12;
e=1.602176565e-19;     % unit in Coulomb
T = 298;                % temperature in K
kbT=1.38e-23*T/e;          % unit: eV
esig = 0.14;        % in eV
%*********************** End Parameter settings **************************

%% =========================== GDM parameters ===========================
fun_h_1=@(x) Nt/(sqrt(2*pi)*esig_interface) * exp(-x.^2/(2*esig_interface^2)) .* (1./(1+exp((x+phi1)/kbT)));
% n1 = integral(fun_h_1,-inf,inf);
n1 = 2.13e26;   %match Fig6
fun_h_2=@(x) Nt/(sqrt(2*pi)*esig_interface) * exp(-x.^2/(2*esig_interface^2)) .* (1./(1+exp((x+phi2)/kbT)));
% n2 = integral(fun_h_2,-inf,inf);
% n2 = 2.13e26;   %match Fig6
n2 = n1*1e-4;
mu0_star=22e-6;        % Units in m^2/(Vs)
C = 0.42;
mu0 = mu0_star * exp(-C*(esig/kbT)^2);
% mu0 = 1e-10;
delta = 2*(log((esig/kbT)^2-esig/kbT)-log(log(4)))/(esig/kbT)^2;

% ================================= End GDM parameters ==============================

warning_flag = 0;
%% ========================================= Enters main loop =============================================

% J = logspace(-2,7,20)';
J = 1e0;
tolerence = 1e-6;
f_upper = 12800;    % tspan parameter in ode45 (see Matlab manual)
f_lower = -12800;   % has to be large so that the cutoff f limit does not exceed the range [f_lower,f_upper]

%% ================================== ENTER TYPE II SOLUTION SOLVER: ==================================
% dimensionless parameters defined in Appendix B
gamma1 = e*L^2/(epi0*epir*kbT)*n1;
gamma2 = e*L^2/(epi0*epir*kbT)*n2; 
i = 1/kbT^2*L^3/(epi0*epir*mu0)*J;
y1 = gamma1/i^(2/3);
y2 = gamma2/i^(2/3);
ymin_1 = 0.1*y2;
ymin_scaling_factor = 10;

coefficients = [delta,kbT,mu0,L,Nt,epi0,epir,e,a,esig,J,Inf,ymin_1,y1,y2,f_lower,f_upper];
[i1,fmin_1,warning_flag,vpa_flag] = calc_i_type2_ode45_GDM(coefficients,1);
if warning_flag~=0
    fprintf('Serious Warning: f range [f_lower, f_upper] not large enough, or y2 too big!\n');
end
fmin_1

coeff_1 = epi0*epir*kbT/(e*L^2)*i^(2/3);         % n = coeff_1 * y;
coeff_2 = i^(1/3)/L*kbT;                        % F = coeff_2 * f;
g1_pre = 0.5*((esig/kbT)^2-esig/kbT);           % prefactor for g1 inside exponential term
g2_pre = 0.44*((esig/kbT)^1.5-2.2);             % prefactor for g2
dg1_pre = epi0*epir*kbT/(Nt*e*L^2) * i^(2/3) * ((esig/kbT)^2-esig/kbT) * delta;     % prefactor for dg1/dy
dg2_pre = kbT/L*i^(1/3) * g2_pre * 0.8*(a/esig)^2;

f_try = (linspace(0.2,10,100))';
func = @(f) exp(g1_pre*(2*coeff_1*ymin_1/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*ymin_1^3 ...,
           + ( 1 +ymin_1/exp(g1_pre*(2*coeff_1*ymin_1/Nt)^delta)*dg1_pre*(2*coeff_1*ymin_1/Nt)^(delta-1)*exp(g1_pre*(2*coeff_1*ymin_1/Nt)^delta)) ...,
           * ( ymin_1*f - 1/(exp(g1_pre*(2*coeff_1*ymin_1/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))) ) ...,
           + ymin_1^2/exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*f)^2) * coeff_2*f *exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1));
for q=1:length(f_try)
    func_value(q,1) = func(f_try(q,1));
end

zhi = [f_try, func_value];


% ================================== END TYPE II SOLUTION SOLVER ==================================
% result

toc;