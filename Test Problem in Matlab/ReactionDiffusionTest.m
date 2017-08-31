% Reaction Diffusion equation
close all
clear all

L = 3.5;
N = 401;
dx = 0.00875;
x = 0:dx:L;

delta = 1e-12;
epsilon = 0.01;
gamma = 100;
lambda = 0.5*sqrt(2*gamma/epsilon);

y0 = zeros(N,1);
for i = 1:N
    y0(i) = 1/(1+exp(lambda*(x(i)-1)));
end

tspan = [0 0.5];

options = odeset('RelTol',1e-12);
[t, Y] = ode23('ReactionDiffusion', tspan, y0, options);
n = size(Y);
y = Y(n(1),:);
plot(x,y);
y(120)
y(130)


