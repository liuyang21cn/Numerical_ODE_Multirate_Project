% Burger's equation
% Lax-Friedrichs Method
% Yang Liu
%  dy/dt + d(y^2/2)/dx = 0
%
%     0 < t < T = 0.6,  ?1 < x < 1
%     u(x, 0) =   1  if -0.6 < x < 0
%                 0   otherwise
%     u(-1,t) = 0
%     u(1, t) = 0;

clear all
close all

h = 1/640;
x = -1:h:1;
nx = length(x);
u = zeros(nx ,1);

t0 = 0;
tf = 0.6;
nu = 1;
dt = nu*h;
nt = (tf-t0)/dt+1;

for i = 1:nx
    xx = -1+(i-1)*h;
    if xx > -0.6 && xx < 0
        u(i) = 1;
    else
        u(i) = 0;
    end
end

plot(x, u, 'r');
pause
legend('Initial condition', 'Solution at t = 0.6');
xlabel('x');ylabel('u');title('Burger''s equation');

for i = 1:nt
    for j = 2:nx-1
        
        unew(j) = (u(j+1)+u(j-1))/2 - dt/h/4*(u(j+1)*u(j+1)-u(j-1)*u(j-1));
    end
    
    u = unew;
    u(1) = 0;
    u(nx) = 0;
    plot(x, u);
    t = (i-1)*dt;
    fprintf('at t = %d\n', t);
    pause
end
