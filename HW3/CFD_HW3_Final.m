%% Vorticity contours
% Plot individually
clear; clc;
mu = 0.01;
dx = 0.01;
dt = 0.001;
domain_size = 1;
k=2;
T = 1.5;
bc = 2;
[u,v] = gauss(mu, dt, dx, T, bc, k);
w = vorticity(u,v,dx,k);
w(w<1e-5)=0;
[X, Y] = meshgrid(linspace(0,k*domain_size,size(w, 2)), linspace(0,1,size(w, 1)));
maxValue = max(max(w));
minValue = min(min(w));
levels = linspace(minValue, maxValue, 8);

figure;
contour(X, Y, w, levels);
cd = colorbar;
axis equal
xticks(0:0.5:6);
yticks(0:0.5:3);
grid on;
caxis([minValue maxValue]);

title([sprintf('t = %.1f, ',T), '\Deltax = ', sprintf('%.3f, ',dx), '\mu = ', sprintf('%.3f, ',mu)]);
xlabel('x');
ylabel('y');

%% Plot together
clear; clc;
mu = 0.001;
dx = 0.02;
dt = 0.001;
k=2;
domain_size = 1;
for i = 1 : 5
    for j = 1 : 2
        bc = j;
        T = 0.5 * (i - 1);
        subplot(5,2, (i-1)*2 + j)
        [u,v] = gauss(mu, dt, dx, T, bc, k);
        w = vorticity(u,v,dx,k);
        [X, Y] = meshgrid(linspace(0,2*domain_size,size(w, 2)), linspace(0,domain_size,size(w, 1)));
        maxValue = max(max(w));
        minValue = min(min(w));
        levels = linspace(minValue, maxValue, 10);
        caxis([minValue maxValue]);
        contour(X, Y, w, levels);
        colorbar;
        axis equal
        xticks(0:0.5:6);
        yticks(0:0.5:3);
        grid on;
        title([sprintf('t = %.1f, ',T), '\Deltax = ', sprintf('%.3f, ',dx), '\mu = ', sprintf('%.3f, ',mu)]);
        xlabel('x');
        ylabel('y');
    end
end

%% Modified wave numbers
beta = linspace(0, pi, 100);
k_exact = beta;
ki = sin(beta);
kr = -(2-2*cos(beta));
figure;
plot(beta,k_exact,'linewidth',1)
hold on
plot(beta,ki,'linewidth',1)
hold on
plot(beta,kr,'linewidth',1)

x = linspace(0, pi, 100);
y = 0.*x;
hold on
plot(x,y,'k-','linewidth',1)

xlim([0 pi])
xlabel('\beta (k\Deltax)')
ylabel('k''\Deltax')
title('Modified wave numbers of central difference scheme',' for convection and diffusion terms');
legend('k_{exact}','k_{real}','k_{imaginary}')

%% Functions
function [u,v] =gauss(mu, dt, dx, T, bc, k)
N = 1/dx;
dy = dx;
x = (0:k*N).*dx;
y = (0:N).*dy;
Vt = 0.25;
x0 = 0.5;
y0 = 0.5;
r0 = 0.1;
R = mu*dt/dx^2;
c = dt/(dx*2);

for i = 1:N+1
    for j = 1:k*N+1
        r(i,j) = sqrt((x(j)-x0)^2+(y(i)-y0)^2);
        u_init(i,j) = 1-Vt*(y(i)-y0)*exp((1-(r(i,j)/r0)^2)/2);
        v_init(i,j) = Vt*(x(j)-x0)*exp((1-(r(i,j)/r0)^2)/2);
    end
end

if bc == 1
    u_init(1,1:k*N+1) = 1;
    u_init(N+1,1:k*N+1) = 1;
    u_init(1:N+1,1) = 1;
    u_init(1:N+1,k*N+1) = 1;

    v_init(1,1:k*N+1) = 0;
    v_init(N+1,1:k*N+1) = 0;
    v_init(1:N+1,1) = 0;
    v_init(1:N+1,k*N+1) = 0;
end

if bc == 2
    u_init(1,1:k*N+1) = 1;
    u_init(N+1,1:k*N+1) = 1;
    u_init(1:N+1,1) = 1;
    u_init(1:N+1,k*N+1) = u_init(1:N+1,k*N);

    v_init(1,1:k*N+1) = 0;
    v_init(N+1,1:k*N+1) = 0;
    v_init(1:N+1,1) = 0;
    v_init(1:N+1,k*N+1) = v_init(1:N+1,k*N);
end

u=u_init;
unew=u;
v=v_init;
vnew=v;
n = T/dt;

for t = 1:n
    rtot=1;
    rtota=[];

    vrtot=1;
    vrtota=[];
    
    % u velocity
    uk = u;
    while rtot >= 1e-5

        rtot = 0;
        temp=uk;
        for i = 2:N
            for j = 2:k*N
                temp(i,j)=(u(i,j) ...
                    -(u(i,j)*(u(i,j+1)-u(i,j-1)))*c ...
                    -(v(i,j)*(u(i+1,j)-u(i-1,j)))*c ...
                    + R*(temp(i-1,j)+temp(i+1,j)+temp(i,j+1)+temp(i,j-1)))/(1+4*R);
                if bc == 2
                    temp(i,end) = temp(i,end-1);
                end
            end
        end
        uk_1 = temp;
        for i = 2:N
            for j = 2:k*N
                res=abs(uk_1(i,j)-uk(i,j));
                rtot = rtot + res;
            end
        end
        rtota=[rtota,rtot];
        uk = uk_1;
    end
    u = uk_1;
    
    % v velocity
    vk = v;
    while vrtot >= 1e-5

        vrtot = 0;
        vtemp = vk;
        for i = 2:N
            for j = 2:k*N
                vtemp(i,j)=(v(i,j) ...
                    -(u(i,j).*(v(i,j+1)-v(i,j-1)))*c ...
                    -(v(i,j).*(v(i+1,j)-v(i-1,j)))*c...
                    + R*(vtemp(i-1,j)+vtemp(i+1,j)+vtemp(i,j+1)+vtemp(i,j-1)))/(1+4*R); 
                if bc == 2
                    vtemp(i,end) = vtemp(i,end-1);
                end
            end
        end
        vk_1 = vtemp;
        for i = 2:N
            for j = 2:k*N
                res=abs(vk_1(i,j)-vk(i,j));
                vrtot = vrtot + res;
            end
        end
        vrtota=[vrtota,vrtot];
        vk = vk_1;
    end
    v = vk_1;

end

end

function w = vorticity(u,v,dx,k)
    N = 1/dx;
    dy = dx;
    w = zeros(N+1,k*N+1);
    for i = 2:N
        for j = 2:k*N
            dudy(i,j) = (u(i+1,j)-u(i-1,j))/2/dy;
            dvdx(i,j) = (v(i,j+1)-v(i,j-1))/2/dx;
            w(i,j) = dvdx(i,j)-dudy(i,j);
        end
    end

end