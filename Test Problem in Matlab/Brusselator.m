function [f] = Brusselator(t,y)

% a = 1.2;
% b = 2.5;
% epi = 0.01;

% a = 1.0;
% b = 3.5;
% epi = 0.01;

a = 0.5;
b = 3;
epi = 0.01;


u = y(1);
v = y(2);
w = y(3);

f = zeros(3,1);
f(1) = a - (w+1.0)*u + v*u*u;
f(2) = w*u - v*u*u;
f(3) = (b-w)/epi - w*u;


end