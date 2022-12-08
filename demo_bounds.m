% Dirichlet kernel bounds for "Bounded Magnitude DFT" - IEEE Signal Processing
% Magazine Tipps and Tricks
% 
% Sebastian J. Schlecht, Thursday, 08 December 2022
clear; clc; close all;

%% Parameters
N2 = 15;
N = N2*2+1;
piN = pi/N;

%% DTFT frequency points
wk = linspace(-N2,N2+1-0.001,1021).'; % w - k
k = ceil(-wk); % -wk = (k - w)
w = wk + k;

% set indices for odd and even
even = mod(k,2) == 0;
odd = mod(k,2) == 1;
nyOdd = (mod(N2,2)==1 & (k == -N2));
nyEven = (mod(N2,2)==0 & (k == -N2));
zero = (k == 0 | k == 1);

% Sinusoidal interpolation
sinF = sin(pi*(wk - floor(wk)));

% Exact Dirichlet kernel weights
Dexact = (-1).^k ./ sin(pi/N.*(wk)) / N;

% Upper and lower bound to the kernel
DUpper = 0*Dexact;
DUpper(even) = 1 ./ sin(pi/N.*(0 - k(even))) / N ;
DUpper(odd) = -1 ./ sin(pi/N.*(1 - k(odd))) / N ;
DUpper(nyOdd & odd) = -1 ./ sin(pi/N.*(0.5 - k(nyOdd & odd))) / N ;
DUpper(zero) = 0;

DLower = 0*Dexact;
DLower(even) = 1 ./ sin(pi/N.*(1 - k(even))) / N ;
DLower(odd) = -1 ./ sin(pi/N.*(0 - k(odd))) / N ;
DLower(nyEven & even) = 1 ./ sin(pi/N.*(0.5 - k(nyEven & even))) / N ;

% remove neighbors for upper bound
Dexact0 = Dexact;
Dexact0(zero) = 0;

deleteCloseNeighbor = 0*wk;
deleteCloseNeighbor(floor(wk) == 0) = nan;

% check that the upper and lower bounds are actually respected
min(DUpper - Dexact0,[],1) % should be 0
max(DLower - Dexact,[],1) % should be 0


%% plot
plotOptions = {'LineWidth',2};
dotOptions = {'.', 'MarkerSize',40};

figure(1); hold on; grid on;
colorOrder = get(gca, 'ColorOrder');
set(gca, 'ColorOrder', colorOrder([1 3 2],:))

plot(wk, Dexact, '-', 'LineWidth',4)
plot(wk, DLower, '-', plotOptions{:})
plot(wk, DUpper + deleteCloseNeighbor, '-', plotOptions{:})

legend({'$D_N(f - f_k) / \sin(\pi f)$','$L_N(f_k)$','$U_N(f_k)$'},'interpreter','latex')
xlabel('Frequency (bins)','interpreter','latex')
ylabel('Amplitude (linear)','interpreter','latex')
ylim([-0.4 1])
xlim([0 max(wk)])
set(gcf,'Position',1*[100 100 400 200])
% exportgraphics(gcf,'./plots/boundedComponent.png','Resolution',600)

figure(2); hold on; grid on;
set(gca, 'ColorOrder', colorOrder([1 3 2],:))

plot(wk, Dexact .* sinF, 'LineWidth',4)
plot(wk, DLower .* sinF, '-', plotOptions{:})
plot(wk, DUpper .* sinF + deleteCloseNeighbor, '-', plotOptions{:})

legend({'Dirichlet Kernel','Lower Bound Kernel','Upper Bound Kernel'},'interpreter','latex')
xlabel('Frequency (bins)','interpreter','latex')
ylabel('Amplitude (linear)','interpreter','latex')
ylim([-0.4 1])
xlim([0 max(wk)])
set(gcf,'Position',1*[100 100 400 200])
% exportgraphics(gcf,'./plots/DirichletUpperLower.png','Resolution',600)
