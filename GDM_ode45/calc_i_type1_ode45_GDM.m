function [i0,value,warning_flag] = calc_i_type1_ode45_GDM(coefficients,f_or_y)

% Get constants from input array
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
f_lower = coefficients(16);
f_upper = coefficients(17);

warning_flag = 0;
ode_rel_tol = 1e-12;
ode_abs_tol = 1e-14;


i = 1/kbT^2*L^3/(epi0*epir*mu0)*J;
coeff_1 = epi0*epir*kbT/(e*L^2)*i^(2/3);         % n = coeff_1 * y;
coeff_2 = i^(1/3)/L*kbT;                % F = coeff_2 * f;


dy_df = @(f,y) f - 1/exp(0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*y/Nt)^delta) ...,
                / exp(0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) / y;
if f_or_y == 0
    syms y
    if fmin<1e-1
        ymin = vpasolve( 0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*y/Nt)^delta ...,
                    + 0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1) + log(y) == -log(fmin), y);
    else
        ymin = vpasolve( exp(0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*y/Nt)^delta ...,
                    + 0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1))*y == 1/fmin,y);
    end
%     ymin = vpasolve( 0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*y/Nt)^delta ...,
%                 + 0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1) + log(y) == -log(fmin), y);
    ymin = double(ymin);
    value = ymin;
elseif f_or_y == 1
    syms f
    fmin = vpasolve(f - 1/exp(0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*ymin/Nt)^delta) ...,
                / exp(0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) / ymin == 0, f);
    fmin = double(fmin);
    value = fmin;
else
    fpintf('Fatal Error: Input condition not recognized!!!\n');
end


overFcn2 = @(f,y) MyEventFcn(f,y,y2);
options2=odeset('Events',overFcn2,'RelTol',ode_rel_tol,'AbsTol',ode_abs_tol);
% (2*coeff_1*ymin/Nt)^delta
%exp(0.5*((esig/kbT)^2-esig/kbT)*(2*coeff_1*ymin/Nt)^delta)
%coeff_1*ymin/Nt
%exp(0.44*((esig/kbT)^1.5-2.2) * (sqrt(1+0.8*(a/esig*coeff_2*fmin)^2)-1))
%a/esig*coeff_2*fmin

[f_plus, y_plus] = ode45(dy_df, [fmin f_upper], ymin, options2);

overFcn1 = @(f,y) MyEventFcn(f,y,y1);
options1=odeset('Events',overFcn1,'RelTol',ode_rel_tol,'AbsTol',ode_abs_tol);
[f_minus, y_minus] = ode45(dy_df, [fmin f_lower], ymin, options1);

f_result = [flipud(f_minus);f_plus];
y_result = [flipud(y_minus);y_plus];

if abs((y_result(end)-y2)/y2)>0.05 || abs((y_result(1)-y1)/y1)>0.05
    warning_flag = 1;
    fprintf('y1 is: %f, and y_result(1) is %f\n',y1,y_result(1));
    fprintf('y2 is: %f, and y_result(end) is %f\n',y2,y_result(end));
end
i0 = (trapz(f_result,1./y_result))^3;


function [value,isterminal,direction] = MyEventFcn(~,y,y_boundary)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1)-y_boundary; % detect y-y_boundary = 0
isterminal = 1; % stop the integration
direction = 0; % ALL direction