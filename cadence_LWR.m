clear all
close all
% Parameters
T = 3; % LWR starts blow up at T ~ 1.65
M = 512;
N = M;
dx = 2*pi/N;
dt = T/M;
L = pi/32; % l=2*pi means LWR
a = round(L/dx);
t = [0:M]/M*T;
x = [1:N]'/N*2*pi;
% Solutions
rho = zeros(N,M+1);
% Initial condition
%rho(:,1)=sin(x);
% Initialize an empty matrix to store the functions
functions = zeros(M, 10);


%%%%%%%%%%%%%%%% Create 10 smooth functions and store them in the 'functions' matrix
initials = GenerateICs(10,N,0.3,0.8);
Data = NaN(8, M, N+1);  %%%%%%%%%%%%%%%%%%%%%DATA(I,J,K) IS I-TH INITIAL CONDITION, J-TH X VALUE, K-TH TIME STEP)
%{ [IF YOU WANT TO PLOT INITIAL CONDITIONS]
% Plot the 10 smooth functions
figure;
hold on;
for i = 1:10
    plot(x, initials(i,:));
end
hold off;

% Set axis limits
xlim([0, 2*pi]);
ylim([0.3, 0.8]);

% Add labels and title
xlabel('x');
ylabel('f(x)');
title('10 Smooth Functions in [0, 2\pi]');

% Add a legend
legend('Function 1', 'Function 2', 'Function 3', 'Function 4', 'Function 5', ...
    'Function 6', 'Function 7', 'Function 8', 'Function 9', 'Function 10');
%}

%%%%%%%%%%%%%%%%%%%%%%%%% BELOW: NUMERICAL SCHEME
for init = 1:8
    
    %DEFINE INITIAL CONDITION
    rho(:,1)=initials(init,:);       %0.8+ 0.1*cos(x);
    %plot(x,rho(:,1), 'DisplayName','Initial')
    %hold on
    
    % Weight
    w = zeros(1,a);
    w(:) = 1/L;
    %w = 2/l-2/l^2*([0:a-1]/a*l);
    % Kernel
    wr = zeros(1,N);
    wr(1:a) = w;
    wc = zeros(1,N);
    wc(1) = w(1);
    wc(N:-1:N-a+2)=w(2:a);
    W=dx*toeplitz(wc,wr);
    
    
    % Evolution
    for k=1:M %k is time step
    % Finite difference scheme
    % rho(:,k+1) = rho(:,k) - dt/dx*(rho([2:N,1],k).*(1-rho([2:N,1],k)) - rho(:,k).*(1-rho(:,k)));
    % Upwind scheme
    %rho(:,k+1) = rho(:,k) - dt/dx*(rho([2:N,1],k).*(1-rho([2:N,1],k)) - rho(:,k).*(1-rho(:,k))).*(1-2*rho(:,k)<0) - dt/dx*(rho(:,k).*(1-rho(:,k)) - rho([N,1:N-1],k).*(1-rho([N,1:N-1],k))).*(1-2*rho(:,k)>=0);
    % Finite volume scheme
    % Upwind flux
    Fdiff = (rho([2:N,1],k).*(1-rho([2:N,1],k)) - rho(:,k).*(1-rho(:,k)))./(rho([2:N,1],k)-rho(:,k));
    F = rho(:,k).*(1-rho(:,k)).*(Fdiff>=0)+rho([2:N,1],k).*(1-rho([2:N,1],k)).*(Fdiff<0);
    % Nonlocal slowdown
    %F = F.*exp(-W*rho(:,k));
    %LWR
    %F = F.*exp(-1/2);
    % FV scheme
    rho(:,k+1) = rho(:,k) - dt/dx*(F-F([N,1:N-1]));
    % if mod(k,100) == 0
    %    plot(x,rho(:,k+1));
    %    hold on;
    %end
    
    
    end
    % Plot the solution at T
    plot(x,rho(:,M+1));
    hold on
    Data(init,:,:) = rho;
end
save('LWRData-init-xval-tval.mat', 'Data')