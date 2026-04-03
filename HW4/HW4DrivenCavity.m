% CFD HW4 Driven Cavity
clear
l = 1;
Lx = 1;
Ly = 1;
N = 128;
tolerance = 1e-8;
dt = 0.001;
Re1 = 100;
[u,v,uold,vold,ustar,vstar,Ustar,Vstar,Uold,Vold,U,V,X,Y,p,dx,dy] = discertize(l,N);
error = 1;
iter = 0;
%% Calculation 
tic
while error>tolerance
    % Step 1 Finding u* & v* with u,v,U,v,uold, vold, Uold, and Vold
    ustar = uGauss(ustar,u,uold,U,V,Uold,Vold,dx,Re1,dt,tolerance);
    vstar = vGauss(vstar,v,vold,U,V,Uold,Vold,dx,Re1,dt,tolerance);
    Ustar = U_star(ustar,Ustar);
    Vstar = V_star(vstar,Vstar);
    % Step 2 Calculate p n+1 by u* v* U* and V*
    pnew = pressure(p,Ustar,Vstar,dx,dt,tolerance);
    % Step 3 Update u,v,U,V at n+1 by u* v* U* and V* and P n+1
    uplus1 = uplusone(ustar,pnew,dx,dt);
    vplus1 = vplusone(vstar,pnew,dy,dt);
    Uplus1 = Uplusone(Ustar,pnew,dx,dt);
    Vplus1 = Vplusone(Vstar,pnew,dy,dt);
    % Calculating Error 
    error = residual(uplus1,u)+residual(vplus1,v)+residual(pnew,p);
    iter = iter +1;
    % u,v,U,V used in this iteration will be u n-1 and v n-1 in next
    % iteration
    uold = u;
    vold = v;
    Uold = U;
    Vold = V;
    % Updated u,v,U,V will be u,v,U,V at n in next iteration
    u = uplus1;
    v = vplus1;
    U = Uplus1;
    V = Vplus1;
    % undate pressure
    p = pnew;   
end
toc
%% Graphing
hold on
streamslice(X,Y,u,v)
hold off
xlim([0 1])
ylim([0 1])
xlabel('X')
ylabel('Y')
title(['Lid Driven Cavity with N = ' num2str(N) ' Re =' num2str(Re1)])

%% Function Part
function [u,v,uold,vold,ustar,vstar,Ustar,Vstar,Uold,Vold,U,V,X,Y,p,dx,dy] = discertize(l,N)
    dx = l/N;
    dy = dx;
    x = linspace(-dx/2,l+dx/2,N+2);
    y = linspace(-dy/2,l+dy/2,N+2);
    [X,Y]=meshgrid(x,y);
    u = 0*X;
    v = 0*Y;
    uold = 0*X;
    vold = 0*Y;
    ustar = 0*X;
    vstar = 0*Y;
    p = 0*X;
    u = ghostu(u);
    v = ghostv(v);
    ustar = ghostu(ustar);
    vstar = ghostv(vstar);
    Ustar=(ustar(2:end-1,1:end-1)+ustar(2:end-1,2:end))/2;
    Vstar=(vstar(1:end-1,2:end-1)+vstar(2:end,2:end-1))/2;
    Uold = Ustar;
    Vold = Vstar;
    U = Ustar;
    V = Vstar;
end

function res = residual(u_new,u_old)
    udiff = u_new(2:end-1,2:end-1) - u_old(2:end-1,2:end-1);
    res = sqrt(sum(sum(udiff.^2)))/(size(udiff,1)*size(udiff,2));
end

function u = ghostu(u)
    u(end,2:end-1) = 2 - u(end-1,2:end-1);
    u(1,2:end-1) = 0 - u(2,2:end-1);
    u(2:end-1,1) = 0 - u(2:end-1,2);
    u(2:end-1,end) = 0 - u(2:end-1,end-1);
end

function v = ghostv(v)
    v(end,2:end-1) = 0 - v(end-1,2:end-1);
    v(1,2:end-1) = 0 - v(2,2:end-1);
    v(2:end-1,1) = 0 - v(2:end-1,2);
    v(2:end-1,end) = 0 - v(2:end-1,end-1);
end

function [fnew] = uGauss(fstar,f,fold,U,V,Uold,Vold,dx,Re,dt,tolerance)
    res = 1;
    fnew = fstar;
    sz = size(fnew);
    NL = nonlinear(f,fold,U,V,Uold,Vold,dx);
    c = constant(f,Re,dt,dx);
    while res > tolerance
        fold = fnew;
        for i=2:sz(1)-1
            for j=2:sz(2)-1
                fnew(i,j) = ((dt/(2*Re*dx^2))*(fnew(i+1,j)+fnew(i-1,j)+fnew(i,j+1)+fnew(i,j-1))...
                    +c(i,j)+dt*NL(i,j))/(1+(4*dt)/(2*Re*dx^2));
            end
        end
        fnew = ghostu(fnew);
        res = residual(fold,fnew);
    end   
end

function [fnew] = vGauss(fstar,f,fold,U,V,Uold,Vold,dx,Re,dt,tolerance)
    res = 1;
    fnew = fstar;
    sz = size(fnew);
    NL = nonlinear(f,fold,U,V,Uold,Vold,dx);
    c = constant(f,Re,dt,dx);
    while res > tolerance
        fold = fnew;
        for i=2:sz(1)-1
            for j=2:sz(2)-1
                fnew(i,j) = ((dt/(2*Re*dx^2))*(fnew(i+1,j)+fnew(i-1,j)+fnew(i,j+1)+fnew(i,j-1))...
                    +c(i,j)+dt*NL(i,j))/(1+(4*dt)/(2*Re*dx^2));
            end
        end
        fnew = ghostv(fnew);
        res = residual(fold,fnew);
   end   
end

function constant = constant(f,Re,dt,dx)
    constant = 0*f;
    sz = size(constant);
    for i = 2:sz(1)-1
        for j = 2:sz(2)-1
            constant(i,j) = (1-4*dt/(2*Re*dx^2))*f(i,j)...
                +(dt/(2*Re*dx^2))*(f(i,j+1)+f(i,j-1)+f(i+1,j)+f(i-1,j));
        end
    end
end

function NL = nonlinear(f,fl,U,V,Uold,Vold,dx)
    NL = 0*f;
    sz = size(NL);
    for i = 2:sz(1)-1
        for j = 2:sz(2)-1
            NL(i,j) = (-1.5*(U(i-1,j)*(f(i,j+1)+f(i,j))/2-U(i-1,j-1)*(f(i,j-1)+f(i,j))/2)/dx) ...
                +(-1.5*(V(i,j-1)*(f(i+1,j)+f(i,j))/2-V(i-1,j-1)*(f(i-1,j)+f(i,j))/2)/dx)...
                +(0.5*(Uold(i-1,j)*(fl(i,j+1)+fl(i,j))/2-Uold(i-1,j-1)*(fl(i,j-1)+fl(i,j))/2)/dx)...
                +(0.5*(Vold(i,j-1)*(fl(i+1,j)+fl(i,j))/2-Vold(i-1,j-1)*(fl(i-1,j)+fl(i,j))/2)/dx);
        end
    end
end

function Fstar = U_star(fstar,Fstar)
    Fsz = size(Fstar);
    for i = 1:Fsz(1)
        for j = 1:Fsz(2)
            Fstar(i,j) = (fstar(i+1,j)+fstar(i+1,j+1))/2;
        end
    end
end

function Fstar = V_star(fstar,Fstar)
    Fsz = size(Fstar);
    for i = 1:Fsz(1)
        for j = 1:Fsz(2)
            Fstar(i,j) = (fstar(i+1,j+1)+fstar(i,j+1))/2;
        end
    end
end

function pnew = pressure(p,Ustar,Vstar,dx,dt,tolerance)
    err = 1;
    pnew = p;
    sz = size(p);
    while err > tolerance
        pold = pnew;
        for i = 2:sz(1)-1
            for j = 2:sz(2)-1
                pnew(i,j) = (pnew(i,j+1)+pnew(i,j-1)+pnew(i+1,j)+pnew(i-1,j))/4 ...
                -(dx/dt)*(Ustar(i-1,j)-Ustar(i-1,j-1)+Vstar(i,j-1)-Vstar(i-1,j-1))/4;
            end
        end
        pnew = ghostp(pnew);
        err = residual(pold,pnew);
    end
end

function p = ghostp(p)
    p(end,2:end-1) = p(end-1,2:end-1);
    p(1,2:end-1) = p(2,2:end-1);
    p(2:end-1,1) = p(2:end-1,2);
    p(2:end-1,end) = p(2:end-1,end-1);
end

function uplus1 = uplusone(ustar,p,dx,dt)
    uplus1 = ustar;
    sz = size(ustar);
    for i = 2:sz(1)-1
        for j = 2:sz(2)-1
            uplus1(i,j) = ustar(i,j) -dt*((p(i,j+1)-p(i,j-1))/(2*dx));
        end
    end
    uplus1 = ghostu(uplus1);
end

function vplus1 = vplusone(vstar,p,dy,dt)
    vplus1 = vstar;
    sz = size(vstar);
    for i = 2:sz(1)-1
        for j = 2:sz(2)-1
            vplus1(i,j) = vstar(i,j) -dt*((p(i+1,j)-p(i-1,j))/(2*dy));
        end
    end
    vplus1 = ghostv(vplus1);
end

function Uplus1 = Uplusone(Ustar,p,dx,dt)
    Uplus1 = Ustar;
    sz = size(Ustar);
    for i = 1:sz(1)
        for j = 1:sz(2)
            Uplus1(i,j) = Ustar(i,j) -dt*((p(i+1,j+1)-p(i+1,j))/dx);
        end
    end
end

function Vplus1 = Vplusone(Vstar,p,dy,dt)
    Vplus1 = Vstar;
    sz = size(Vstar);
    for i = 1:sz(1)
        for j = 1:sz(2)
            Vplus1(i,j) = Vstar(i,j) -dt*((p(i+1,j+1)-p(i,j+1))/dy);
        end
    end
end

