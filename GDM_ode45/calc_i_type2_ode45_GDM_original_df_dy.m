function [i0,value,warning_flag,vpa_flag] = calc_i_type2_ode45_GDM_original_df_dy(coefficients,f_or_y)

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

warning_flag = 0;
vpa_flag = 0;
i0 = NaN;
value = NaN;

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
    value = ymin;
elseif f_or_y == 1
    syms f
    fmin = vpasolve( exp(g1_pre*(2*coeff_1*ymin/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*ymin^3 ...,
               + ( 1 +ymin/exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)*dg1_pre*(2*coeff_1*ymin/Nt)^(delta-1)*exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)) ...,
               * ( ymin*f - 1/(exp(g1_pre*(2*coeff_1*ymin/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))) ) ...,
               + ymin^2/exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*f)^2) * coeff_2*f *exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) == 0, f);
           
    if size(fmin,1)==0
        fprintf('vpasolve failed!\n')
        vpa_flag=1;
        warning_flag=1;
        return  
    end
    fmin = double(fmin);
    value = fmin;
else
    fpintf('Fatal Error: Input condition not recognized!!!\n');
end
% ================================================== end solving for inflection point ========================================================
% dy2_df2 = @(f,y) exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y^3 ...,
%            + ( 1 +y/exp(g1_pre*(2*coeff_1*y/Nt)^delta)*dg1_pre*(2*coeff_1*y/Nt)^(delta-1)*exp(g1_pre*(2*coeff_1*y/Nt)^delta)) ...,
%            * ( y*f - 1/(exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))) ) ...,
%            + y^2/exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) * dg2_pre/sqrt(1+0.8*(a/esig*coeff_2*f)^2) * coeff_2*f *exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1));
% second_order = dy2_df2(fmin,ymin)

%% ================================================ ode solver for differential equation ==================================================
if ymin<y2
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
%     fprintf('Case ymin < y2:\n')
%     fprintf('ymin is: %f, and y_result(1) is %f\n',ymin,y_result(1));
%     fprintf('y1 is: %f, and y_result(end) is %f\n',y1,y_result(end));           % since always have y2 < y1 in type II solution
    if abs((y_result(end)-y1)/y1)>0.05
        warning_flag = 1;
        fprintf('y1 is: %f, and y_result(1) is %f\n',y1,y_result(1));
        fprintf('y2 is: %f, and y_result(end) is %f\n\n',y2,y_result(end));           % since always have y2 < y1 in type II solution
    end
else
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
%     first_order = f_minus(end)-1/(exp(g1_pre*(2*coeff_1*y_minus(end)/Nt)^delta) * exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f_minus(end))^2)-1)) *y_minus(end))
%     first_order_min = fmin-1/(exp(g1_pre*(2*coeff_1*ymin/Nt)^delta) * exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1)) *ymin)
    f_result = [flipud(f_minus);f_plus];
    y_result = [flipud(y_minus);y_plus];
    %first_order = df_dy(ymin,fmin)
%     semilogy(f_result,1./y_result)
    
%     fprintf('Case y2 < ymin < y1:\n');
%     fprintf('y1 is: %f, and y_result(1) is %f\n',y2,y_result(1));
%     fprintf('y2 is: %f, and y_result(end) is %f\n',y1,y_result(end));           % since always have y2 < y1 in type II solution
    if abs((y_result(1)-y1)/y1)>0.05 || abs((y_result(end)-y2)/y2)>0.05
        warning_flag = 1;
        fprintf('y1 is: %f, and y_result(1) is %f\n',y1,y_result(1));
        fprintf('y2 is: %f, and y_result(end) is %f\n\n',y2,y_result(end));           % since always have y2 < y1 in type II solution
    end
end
% ====================================================== end ode45 for differential equation ========================================================

% lalala = [f_result y_result]
i0 = (trapz(y_result,exp(g1_pre*(2*coeff_1*y_result/Nt).^delta) .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_result).^2)-1)) ...,
   ./ (exp(g1_pre*(2*coeff_1*y_result/Nt).^delta) .* exp(g2_pre*(sqrt(1+0.8*(a/esig*coeff_2*f_result).^2)-1)) .* y_result .* f_result -1) ))^3;
% i0 = (trapz(f_result, 1./y_result))^3

% semilogy(f_plus,y_plus,'r');
% hold on
% semilogy(f_minus,y_minus,'c');
% pause

