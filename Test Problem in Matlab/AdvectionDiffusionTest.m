% Advection-Diffusion equation

close all
clear all

L = 1;
N = 401;
dx = 5e-3;
x = -1:dx:1;

delta = 1e-12;
a = 5;
d = 0.01;
c = 100;

y0 = zeros(N,1);

tspan = [0 0.8];
options = odeset('RelTol',1e-12);
[t, Y] = ode113('AdvectionDiffusion', tspan, y0, options);
n = size(Y);
y = Y(n(1),:);
plot(x,y);
y(210)
y(220)