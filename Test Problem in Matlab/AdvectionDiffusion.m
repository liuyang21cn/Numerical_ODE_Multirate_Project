function [yt] = AdvectionDiffusion(t, y)

% Advection-Diffusion equation

a = 5;
d = 0.01;
c = 100;

N = length(y);
dx = 2/(N-1);
x = -1:dx:1;

yt = zeros(N,1);
g = zeros(N,1);

for i = 3:N-2
    g(i) = 1000*(cos(pi*x(i)/2))^200*sin(pi*t);
    yt(i) = ( -a*(-y(i+2)+8*y(i+1)-8*y(i-1)+y(i-2))/12/dx + ...
        d*(-y(i+2)+16*y(i+1)-30*y(i)+16*y(i-1)-y(i-2))/12/dx/dx - ...
        c*y(i) + g(i) );
end

yt(1) = 0;
yt(2) = 0;
yt(N-1) = 0;
yt(N) = 0;


end