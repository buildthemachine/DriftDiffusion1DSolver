function [f_result,y_result,f_min,integ] = calc_i_type2_ode45_returnFY_GDM_original_df_dy(coefficients,f_or_y)

%% Get constants from input array
delta = coefficients(1);
kbT = coefficients(2);
mu0 = coefficients(3);
L = coefficients(4);
Nt = coefficients(5);
epi0 = coefficients(6);
epir = coefficients(7);
e = coefficients(8);
a = coefficients(9);
esig = coefficients(10);
J = coefficients(11);
fmin = coefficients(12);
ymin = coefficients(13);
y1 = coefficients(14);
y2 = coefficients(15);

ode_solver = 'ode45';
ode_rel_tol = 1e-12;
ode_abs_tol = 1e-14;

i = 1/kbT^2*L^3/(epi0*epir*mu0)*J;
coeff_1 = epi0*epir*kbT/(e*L^2)*i^(2/3);         % n = coeff_1 * y;
coeff_2 = i^(1/3)/L*kbT;                        % F = coeff_2 * f;
g1_pre = 0.5*((esig/kbT)^2-esig/kbT);           % prefactor for g1 inside exponential term
g2_pre = 0.44*((esig/kbT)^1.5-2.2);             % prefactor for g2
dg1_pre = epi0*epir*kbT/(Nt*e*L^2) * i^(2/3) * ((esig/kbT)^2-esig/kbT) * delta;     % prefactor for dg1/dy
dg2_pre = kbT/L*i^(1/3) * g2_pre * 0.8*(a/esig)^2;

% dy_df = @(f,y) f - 1/exp(g1_pre*(2*coeff_1*y/Nt)^delta) ...,
%                 / exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) / y;          %anomynous function handel for eqn (8) in Mensfoor 2008
df_dy = @(y,f) exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y ...,
            /(exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y*f -1);
% df_dy = @(y,f) 1/(f - 1/(exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) *y) ); 
        
%% ================================================== solving for inflection point ======================================================
if f_or_y == 0
    syms y
    ymin = vpasolve( exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1))*y^3 ...,
               + ( 1 +y/exp(g1_pre*(2*coeff_1*y/Nt)^delta)*dg1_pre*(2*coeff_1*y/Nt)^(delta-1)*exp(g1_pre*(2*coeff_1*y/Nt)^delta)) ...,
               * ( y*fmin - 1/(exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1))) ) ...,
               + y^2/exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1)) * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*fmin)^2) * coeff_2*fmin *exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1)) == 0, y);
           
    ymin = double(ymin);
    f_min = ymin;
elseif f_or_y == 1
    syms f
    fmin = vpasolve( exp(g1_pre*(2*coeff_1*ymin/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*ymin^3 ...,
               + ( 1 +ymin/exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)*dg1_pre*(2*coeff_1*ymin/Nt)^(delta-1)*exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)) ...,
               * ( ymin*f - 1/(exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))) ) ...,
               + ymin^2/exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*f)^2) * coeff_2*f *exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) == 0, f);
    fmin = double(fmin);
    f_min = fmin;
else
    fpintf('Fatal Error: Input condition not recognized!!!\n');
end

%% ================================================ ode solver for differential equation ==================================================
if ymin<y2
    [y_plus, f_plus] = ode45(df_dy, [ymin y1], fmin);
    y_result = y_plus;
    f_result = f_plus;
    
    options = odeset('RelTol',ode_rel_tol,'AbsTol',ode_abs_tol);
    if strcmp(ode_solver, 'ode45')
        [~, f_temp] = ode45(df_dy, [ymin y2], fmin, options);
        [y_minus,f_minus] = ode45(df_dy, [y2 y1],f_temp(end),options);
    elseif strcmp(ode_solver, 'ode15s')
        [~, f_temp] = ode15s(df_dy, [ymin y2], fmin, options);
        [y_minus,f_minus] = ode15s(df_dy, [y2 y1],f_temp(end),options);
    else
        fprintf('Please enter a valid ode solver!\n');
        return
    end
    y_result = y_minus;
    f_result = f_minus;
else
    [y_minus, f_minus] = ode45(df_dy, [ymin y2], fmin);
    [y_plus, f_plus] = ode45(df_dy, [ymin y1], fmin);
    f_result = [flipud(f_minus);f_plus];
    y_result = [flipud(y_minus);y_plus];
    
    options = odeset('RelTol',ode_rel_tol,'AbsTol',ode_abs_tol);
    if strcmp(ode_solver, 'ode45')
        [y_plus, f_plus] = ode45(df_dy, [ymin y2], fmin,options);
        [y_minus, f_minus] = ode45(df_dy, [ymin y1], fmin,options);
    elseif strcmp(ode_solver, 'ode15s')
        [y_plus, f_plus] = ode15s(df_dy, [ymin y2], fmin,options);
        [y_minus, f_minus] = ode15s(df_dy, [ymin y1], fmin,options);
    else
        fprintf('Please enter a valid ode solver!\n');
        return
    end
    f_result = [flipud(f_minus);f_plus];
    y_result = [flipud(y_minus);y_plus];
end

integ = -trapz(y_result,exp(g1_pre*(2*coeff_1*y_result/Nt).^delta) .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_result).^2)-1)) ...,
        .* f_result ./ (1 - exp(g1_pre*(2*coeff_1*y_result/Nt).^delta) .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_result).^2)-1)) ...,
        .* y_result .* f_result));
    
    
    
    