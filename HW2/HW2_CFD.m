clear; clc;
Pe = 50;
dt = 0.0001;
T = 10;
scheme = 1;

%% Convergence test
dx = 1/20;
figure
[xn, un] = FTCS(dt, dx, Pe, 0.01, scheme);
plot(xn, un, 'linewidth', 1)
hold on
[xn, un] = FTCS(dt, dx, Pe, 0.05, scheme);
plot(xn, un, 'linewidth', 1)
hold on
[xn, un] = FTCS(dt, dx, Pe, 0.1, scheme);
plot(xn, un, 'linewidth', 1)
hold on
[xn, un] = FTCS(dt, dx, Pe, 1, scheme);
plot(xn, un, 'linewidth', 1)
hold on
[xn, un] = FTCS(dt, dx, Pe, 2, scheme);
plot(xn, un, '--', 'linewidth', 1)

title('Convergence test of solutions for integration time at \Deltax = 1/20')
xlabel('x')
ylabel('u')
legend('T = 0.01', 'T = 0.05', 'T = 0.1', 'T = 1', 'T = 2')

%% Plot solutions 1
x = linspace(0, 1, 1000);
uexa = (exp(Pe.*x)-1)./(exp(Pe)-1);
figure;
plot(x, uexa, 'linewidth', 1)

dx = 1/10;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

dx = 1/20;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

dx = 1/50;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

dx = 1/100;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

title('Solutions of viscous Burger''s equation using QUICK scheme')
xlabel('x')
ylabel('u')
legend('Exact solution', '\Deltax = 1/10', '\Deltax = 1/20', '\Deltax = 1/50', '\Deltax = 1/100')

%% Plot solutions 2
x = linspace(0, 1, 1000);
uexa = (exp(Pe.*x)-1)./(exp(Pe)-1);
dx = 1/100;
figure;
plot(x, uexa, 'linewidth', 1)

scheme = 1;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

scheme = 2;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

scheme = 3;
[xn, un] = FTCS(dt, dx, Pe, T, scheme);
hold on
plot(xn, un, 'linewidth', 1)

title('Solutions of viscous Burger''s equation using different schemes')
xlabel('x')
ylabel('u')
legend('Exact solution', 'FTCS', '1^st-order upwind', 'QUICK')

%% Functions
function [xn, un] = FTCS(dt, dx, Pe, T, scheme)
N=floor(1/dx);
xn = linspace(0, 1, N+1);
n = T/dt;

un=zeros(N+1,1);
un(N+1) = 1;
unew = un;
C = dt/dx;
r = dt/Pe/dx^2;

for i = 1:n
    for j = 2:N
        if scheme == 1
            unew(j) = (r-C/2)*un(j+1) + (1-2*r)*un(j) + (r+C/2)*un(j-1);
        end
        if scheme == 2
            unew(j) = r*un(j+1) + (1-C-2*r)*un(j) + (r+C)*un(j-1);
        end

    end
    un = unew;
end

if scheme == 3
    ug = 2*un(1)-un(2);
    for i = 1:n
        j=2;
        unew(j) = (r-3*C)*un(j+1) + (1-3*C-2*r)*un(j) + (r-5*C)*un(j-1) + C*ug;
        for j = 3:N
            unew(j) = (r-3*C/8)*un(j+1) + (1-3*C/8-2*r)*un(j) + (r+7*C/8)*un(j-1) + (-C/8)*un(j-2);
        end
        un = unew;
        ug = 2*un(1)-un(2);
    end
end
un = un';

end