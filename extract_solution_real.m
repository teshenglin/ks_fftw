%close all
clear all

NDIM = 127; % NB!!! This value is given in line 4 of b.all

LAB= [1:2:11]; % This the lablel of the solution that we wish to plot

NUM = length(LAB);
uu = zeros(NDIM, NUM);
leg = zeros(NUM,1);

for ii=1:NUM
    sol=load(['sol_LAB_',int2str(LAB(ii)),'.dat']);
    U(2:NDIM) = sol(3:NDIM+1);
    U(1) = -sum(U(2:NDIM));
    uu(:,ii) = U;
    leg(ii) = LAB(ii);
end

x = linspace(0, 2*pi, NDIM+1); x = x(1:NDIM)'; x = x*ones(1,NUM);

plot(x, uu); hold on

legend(num2str(leg))
