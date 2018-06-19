function [density] = calc_density(coef_n,f_end,y_end)

%% Get constants from input array
delta = coef_n(1);
kbT = coef_n(2);
mu0 = coef_n(3);
L = coef_n(4);
Nt = coef_n(5);
epi0 = coef_n(6);
epir = coef_n(7);
e = coef_n(8);
a = coef_n(9);
esig = coef_n(10);
nL = coef_n(11);
J = coef_n(12);

i = 1/kbT^2*L^3/(epi0*epir*mu0)*J;
coeff_1 = epi0*epir*kbT/(e*L^2)*i^(2/3);         % n = coeff_1 * y;
coeff_2 = i^(1/3)/L*kbT;                        % F = coeff_2 * f;
g1_pre = 0.5*((esig/kbT)^2-esig/kbT);           % prefactor for g1 inside exponential term
g2_pre = 0.44*((esig/kbT)^1.5-2.2);             % prefactor for g2

%% faster algorithm
num_points = 2*L/a;       % Total number of points along the n(x)~x curve
index_space = floor(length(f_end)/num_points);

offset=2;
for k=0:num_points
    y_k = k*index_space+offset;
    n(k+1) = coeff_1*y_end(y_k);
    s(k+1) = -trapz(y_end(1:y_k), exp(g1_pre*(2*coeff_1*y_end(1:y_k)/Nt).^delta)...,
            .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_end(1:y_k)).^2)-1)) ...,
            ./ (i^(1/3) * (1-exp(g1_pre*(2*coeff_1*y_end(1:y_k)/Nt).^delta) ...,
            .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_end(1:y_k)).^2)-1)) ...,
            .* f_end(1:y_k) .* y_end(1:y_k) ) ) );
    x(k+1) = s(k+1);
end
s_end = s(end)
density = [x' n'];

%% much slower algorithm
% s = zeros(length(y_end),1);
% n = zeros(length(y_end),1);
% n(1) = nL;
% for p=2:length(s)
%     s(p) = trapz(y_end(1:p), ...,
%            exp(g1_pre*(2*coeff_1*y_end(1:p)/Nt).^delta).*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_end(1:p)).^2)-1)) ...,
%            ./ (i^(1/3)*(exp(g1_pre*(2*coeff_1*y_end(1:p)/Nt).^delta).*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_end(1:p)).^2)-1)) ...,
%            .*y_end(1:p).*f_end(1:p) - 1)) );
%     n(p) = coeff_1*y_end(p);
% end
% density = [s n];
    