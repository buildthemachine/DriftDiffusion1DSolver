% In this version, I still try to calculate y=y(f), rather than f=f(y) as suggested 
% in van Mensfoort paper (2008). The divergence of dy/df does not seem to pose a problem here.

function [f_result,y_result,value] = calc_i_type2_ode45_returnFY_GDM(coefficients,f_or_y)

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

i = 1/kbT^2*L^3/(epi0*epir*mu0)*J;
coeff_1 = epi0*epir*kbT/(e*L^2)*i^(2/3);         % n = coeff_1 * y;
coeff_2 = i^(1/3)/L*kbT;                        % F = coeff_2 * f;
g1_pre = 0.5*((esig/kbT)^2-esig/kbT);           % prefactor for g1 inside exponential term
g2_pre = 0.44*((esig/kbT)^1.5-2.2);             % prefactor for g2
dg1_pre = epi0*epir*kbT/(Nt*e*L^2) * i^(2/3) * ((esig/kbT)^2-esig/kbT) * delta;     % prefactor for dg1/dy
dg2_pre = kbT/L*i^(1/3) * g2_pre * 0.8*(a/esig)^2;

dy_df = @(f,y) f - 1/exp(g1_pre*(2*coeff_1*y/Nt)^delta) ...,
                / exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) / y;          %anomynous function handel for eqn (8) in Mensfoor 2008
% df_dy = @(y,f) exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y ...,
%             /(exp(g1_pre*(2*coeff_1*y/Nt)^delta)*exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1))*y*f -1);
% df_dy = @(y,f) 1/(f - 1/(exp(g1_pre*(2*coeff_1*y/Nt)^delta) * exp(g2_pre * (sqrt(1+0.8*(a/esig*coeff_2*f)^2)-1)) *y) ); 
        
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
    fmin = double(fmin);
    value = fmin;
else
    fpintf('Fatal Error: Input condition not recognized!!!\n');
end

if ymin<y2
    overFcn1 = @(f,y) MyEventFcn(f,y,y1);
    options = odeset('Events',overFcn1);
    [f_minus, y_minus] = ode45(dy_df, [fmin f_lower], ymin, options);
    y_result = flipud(y_minus);
    f_result = flipud(f_minus);
else
    overFcn1 = @(f,y) MyEventFcn(f,y,y1);
    options1 = odeset('Events',overFcn1);
    [f_minus, y_minus] = ode45(dy_df, [fmin f_lower], ymin, options1);
    
    overFcn2 = @(f,y) MyEventFcn(f,y,y2);
    options2 = odeset('Events',overFcn2);
    [f_plus, y_plus] = ode45(dy_df, [fmin f_upper], ymin, options2);
    f_result = [flipud(f_minus);f_plus];
    y_result = [flipud(y_minus);y_plus];
end

function [value,isterminal,direction] = MyEventFcn(~,y,y_boundary)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1)-y_boundary; % detect y-y_boundary = 0
isterminal = 1; % stop the integration
direction = 0; % ALL direction