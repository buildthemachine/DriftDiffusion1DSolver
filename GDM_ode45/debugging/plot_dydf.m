% Ploting dy/df == 0 function, trouble shoot why getting complex numbers
tic;

clear all;
close all;

%fmin = transpose(linspace(1e-6,1e-4,20));
% fmin = 1e-8;
% ymin = [];
% for count=1:length(fmin)
%     if mod(count,4)==0
%         count
%     end
%     syms y
%     % yval = vpasolve(fmin(count) - 1/(y*exp(15*(2.97e-5*y)^0.1708 + 94.072*(sqrt(1+0.8*(0.0064*fmin(count))^2)-1))) ==0,y);
%     % yval = vpasolve( (y*exp(15*(2.97e-5*y)^0.5)) == 1e5,y);
%     yval = vpasolve( log(y)+15*(2.97e-5*y)^0.1708 == log(1/fmin),y);
%     yval = double(yval);
%     ymin = [ymin;yval];
% end
% ymin
%semilogy(fmin,ymin);

fmin = 0.05963;
syms y
yval = vpasolve(15*y^0.1708 + log(y) == -log(fmin)-94*(sqrt(1+0.8*(1.29*fmin)^2)-1),y);
yval = double(yval)

toc;