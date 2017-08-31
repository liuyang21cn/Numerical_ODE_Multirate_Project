close all
clear all
clc

% y0 = [3.9, 1.1, 2.8]';
% y0 = [1.2, 3.1, 3]';
y0 = [3, 3, 3.5]';


tspan = [0 10];
options = odeset('RelTol',1e-12);
[t, Y] = ode113('Brusselator', tspan, y0, options);
plot(t,Y);
title('Brusselator Equation');
xlabel('t');ylabel('Y');
legend('Fast','Fast','Slow');