%close all
clear all

NDIM = 64; % NB!!! This value is given in line 4 of b.all

LAB= [1]; % This the lablel of the solution that we wish to plot

NUM = length(LAB);
uu = zeros(NDIM, NUM);
leg = zeros(NUM,1);

for ii=1:NUM
    sol=load(['sol_LAB_',int2str(LAB(ii)),'.dat']);
    U(2:NDIM) = sol(3:NDIM+1);
    U2(1)=0;
    if(mod(NDIM, 2)==0)
        U2(2:NDIM/2) = U(2:2:NDIM-1)+1i*U(3:2:NDIM-1);
        U2(NDIM/2+1) = real(U(NDIM));
        U2(NDIM/2+2:NDIM) = conj(U2(NDIM/2:-1:2));
    else
        U2(2:(NDIM+1)/2) = U(2:2:NDIM)+1i*U(3:2:NDIM);
        U2((NDIM+1)/2+1:NDIM) = conj(U2((NDIM+1)/2:-1:2));
    end
    uu(:,ii) = real(ifft(U2));
    leg(ii) = LAB(ii);
end

x = linspace(0, 2*pi, NDIM+1); x = x(1:NDIM)'; x = x*ones(1,NUM);

figure(2)
plot(x, uu); hold on

legend(num2str(leg))
