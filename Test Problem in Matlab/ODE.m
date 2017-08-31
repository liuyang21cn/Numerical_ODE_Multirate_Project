close all
clear all
clc

t = linspace(0,1,129);
x = cos(50*t);
y = sin(50*t);
z = 5051/2501*exp(-t) - 49/2501*cos(50*t)+51/2501*sin(50*t);

hold on
plot(t,x, 'r');
plot(t,y, 'b');
plot(t,z, 'g');
hold off
title('System of ode problem');
xlabel('time');ylabel('y');
legend('Fast', 'Fast', 'Slow');