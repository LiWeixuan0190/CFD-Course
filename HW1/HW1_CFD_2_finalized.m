clear; clc;
%% Plot the energy spectra
dx = 2*pi/20;
dt = 0.1;
L = floor(2*pi/dx);                         
figure;

for i = 1:10
    t = i;
    un1_ali = ForwardEuler(dt, dx, 1, t);                      
    
    U = un1_ali;               
    Y = fft(U);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:(L/2))/L * L;
    
    loglog(f, P1)
    hold on

end

title(['Fully aliased energy spectra for equation 1 with each time step ', num2str(dt)]);
xlabel('Wave Number k')
ylabel('Amplitude')
legend('1','2','3','4','5','6','7','8','9','10')

%% Plot the fully aliased and unaliased energy spectra 
dt = 0.01;
L = floor(2*pi/dx); 
x = linspace(2, 2*pi, L+1);
colors = lines(10);

for i = 1:10
    t = i;
    dx = 2*pi/20;
    L = floor(2*pi/dx); 
    un1_ali = ForwardEuler(dt, dx, 1, t); 

    U = un1_ali;               
    Y = fft(U);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_ali = zeros(11,1);
    P1_ali = P1;
    f = (0:(L/2))/L * L;

    figure;
    plot(f, P1_ali, 'ro', 'MarkerSize', 8);
    hold on

    dx = 2*pi/20/100;
    L = floor(2*pi/dx); 
    un1_unali = ForwardEuler(dt, dx, 1, t);
    U = un1_unali;               
    Y = fft(U);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1_unali = zeros(11,1);
    P1_unali(:) = P1(1:11);

    plot(f, P1_unali, 'bx', 'MarkerSize', 8); 
    title(['Energy spectra of fully aliased and unaliased scheme ' sprintf('\n') ' for equation 1 at time step ', num2str(t)], 'FontWeight', 'normal');
    xlabel('Wave number k')
    ylabel('Amplitude')
    legend('Fully aliased', 'Unaliased')
    
end

%% Functions
function [un] = ForwardEuler(dt, dx, eqnum, t)
N=floor(2*pi/dx);
n=t;

if eqnum == 1
    un=zeros(N+1,1);
    for i = 1:N+1
        x = (i-1)*dx;
        un(i) = sin(x) + 0.5*sin(4*x);
    end
    
    for i = 1:n
        for j = 1:N+1
            unew(j)=-un(j)^2*dt+un(j);
        end
        un = unew;
    end
    un = un';
end

if eqnum == 2
    un=zeros(N+1,1);
    for i = 1:N+1
        x = (i-1)*dx;
        un(i) = sin(x) + 0.5*sin(4*x);
    end
    
    for i = 1:n
        for j = 1:N+1
            x = (j-1)*dx;
            unew(j)=-sin(3*x)*un(j)*dt + un(j);
        end
        un = unew;
    end
    un = un';
end

end