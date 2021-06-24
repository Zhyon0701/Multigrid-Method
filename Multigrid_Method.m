close all
clear all

N = 64; L = 1;
h = L/N;
phi = zeros(1,N+1);
f = (sin(pi*[0:N]*h)+sin(16*pi*[0:N]*h))/2;

for cnt = 1:1000
    phi = V_Cycle(phi,f,h);
    r = residual(phi,f,h);
    if max(abs(r)) < 0.001
        break
    end
end

function phi = V_Cycle(phi,f,h)
 % Recursive V-Cycle Multigrid for solving the Poisson equation (\nabla^2 phi = f) on a uniform grid of spacing h

 % Pre-Smoothing
 phi = smoothing(phi,f,h);
 
 % Compute Residual Errors
 r = residual(phi,f,h);
 
 % Restriction
 rhs = restriction(r);

 eps = zeros(size(rhs));

 % stop recursion at smallest grid size, otherwise continue recursion
 if length(eps)-1 == 2
         eps = smoothing(eps,rhs,2*h);
 else        
         eps = V_Cycle(eps,rhs,2*h);        
 end
 
 % Prolongation and Correction
 phi = phi + prolongation(eps);
 
 % Post-Smoothing
 phi = smoothing(phi,f,h);    
end

function res = smoothing(phi,f,h)
    N = length(phi)-1;
    res = zeros(1,N+1);
    for j = 2:N
        res(j) = (phi(j+1)+res(j-1)-h^2*f(j))/2;
    end
end

function res = residual(phi,f,h)
    N = length(phi)-1;
    res = zeros(1,N+1);
    res(2:N) = f(2:N)-(phi(1:N-1)-2*phi(2:N)+phi(3:N+1))/h^2;
end

function res = restriction(r)
    N = (length(r)-1)/2;
    res = zeros(1,N+1);
    for j = 2:N
        res(j) = (r(2*j-2)+2*r(2*j-1)+r(2*j))/4;
    end
end

function res = prolongation(eps)
    N = (length(eps)-1)*2;
    res = zeros(1,N+1);
    for j = 2:2:N
        res(j) = (eps(j/2)+eps(j/2+1))/2;
    end
    for j = 1:2:N+1
        res(j) = eps((j+1)/2);
    end
end

