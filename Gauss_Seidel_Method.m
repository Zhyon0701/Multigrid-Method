close all
clear all

N = 64; L = 1;
dx = L/N;

phi = zeros(1,N+1);
new = zeros(1,N+1);
tmp = (sin(pi*[1:N-1]*dx)+sin(16*pi*[1:N-1]*dx))/2;
r0 = zeros(1,N+1); r10 = zeros(1,N+1); r100 = zeros(1,N+1);
r0(2:N) = tmp;
resi(1) = max(abs(tmp));

for t = 1:100
    for j = 2:N
        new(j) = (phi(j+1)+new(j-1)-dx^2*tmp(j-1))/2;
    end
    new(1) = 0; new(N+1) = 0;
    r = tmp-(new(1:N-1)-2*new(2:N)+new(3:N+1))/dx^2;
    resi(t+1) = max(abs(r));
    phi = new;
    if t == 10
        r10(2:N) = r;
    elseif t == 100
        r100(2:N) = r;
    end
end

figure
plot([0:length(resi)-1],resi,'+-');
xlabel('Number of Iterations')
ylabel('max(|r_j|)')
title('Convergence Curve')

x = linspace(0,1,N+1);
figure
plot(x,r0,'-',x,r10,'+-',x,r100,'x-')
legend('0 iterations','10 iterations','100 iterations')
xlabel('x_j')
ylabel('r_j')
title('r_j against x_j')