function [yt] = ReactionDiffusion(t, y)

N = length(y);
epsilon = 0.01;
gamma = 100;
dx = 3.5/(N-1);

yt = zeros(N,1);
yt(1) = 0;
yt(N) = 0;
for i = 2:N-1
    yt(i) = epsilon/dx/dx*(y(i+1)-2*y(i)+y(i-1)) + gamma*y(i)*y(i)*(1-y(i));
end


end