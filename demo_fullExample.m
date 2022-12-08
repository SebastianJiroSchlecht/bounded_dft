% Full example for "Bounded Magnitude DFT" - IEEE Signal Processing
% Magazine Tipps and Tricks
% 
% Sebastian J. Schlecht, Thursday, 08 December 2022
clear; clc; close all;

%% Parameters
rng(11); % fix random seed
N2 = 16; % lengths of FIR filter
N = N2*2+1; % lenghts of the symmetric FIR filter

firCase = 2; % choose FIR sequence (1 or 2)

%% Define scaled Dirichlet
scaledDiric = @(x) sin(pi .* x) ./ (N .* sin(pi/N.* x));

%% create FIR sequence h
switch firCase
    case 1 % random
        h = randn(N2+1,1); 
        h = h / norm(h) / sqrt(sqrt(N2)) * 1.1;
    case 2 % sinc
        h = sinc(4.1*linspace(-1,1,N2+1).'); 
        h = h / norm(h) / sqrt(sqrt(N2)) * 1.3;
end
hh = conv(h, flipud(h)); % symmetric FIR
h0 = ifftshift(hh); % make it zero phase 
H0 = real(fftshift(fft(h0))); % H = H0
f0 = circspace(N).'/(2*pi)*N; % DFT frequency bins

%% Compute DTFT
num_zeroPad = 301; % choose arbitary high number
h1_zero = [zeros(num_zeroPad,1);hh;zeros(num_zeroPad,1)];
H1_zero = real(fftshift(fft(ifftshift(h1_zero))));

DTFTlen = length(h1_zero);
f1_zero = circspace(DTFTlen).'/(2*pi)*N; 

%% Method
[~,lx,ux] = diricBounds(N2);

% compute interpolation weights
wfrac = linspace(0.001,0.999,100); % number of points (100) between each DFT bin 
w = (wfrac + f0).';
sinw = sin(pi*wfrac);

% compute direct neighbors with actual scalted Dirichlet kernel
neighbor1 = scaledDiric(wfrac);
neighbor2 = scaledDiric(1 - wfrac);

% apply bounds to FIR filter
H2U = real(circshift(fft(ux .* h0),-1));
H2L = real(circshift(fft(lx .* h0),-1));

% construct entire bounds
neighbors = (H0 .* neighbor1 + circshift(H0,-1) .* neighbor2).';
upperSin = (H2U .* sinw).';
lowerSin = (H2L .* sinw).';
upperBound = neighbors + upperSin;
lowerBound = neighbors + lowerSin;

%% plot
plotOptions = {'LineWidth',2};
dotOptions = {'s', 'MarkerSize',9,'LineWidth',2};
figure(1); hold on; grid on;
plot(f0,H0,dotOptions{:});
plot(f1_zero,H1_zero,'-',plotOptions{:});
plot(f0,H0,'-k',plotOptions{:});
legend({'DFT points $X(f_k)$','DTFT $X(f)$','DFT linear interpolation'},'Location','northoutside','NumColumns',3,'interpreter','latex')
xlabel('Frequency (bins)','interpreter','latex')
ylabel('Magnitude (linear)','interpreter','latex')
xlim([0,N2])
set(gcf,'Position',[100 100 400 200])
% exportgraphics(gcf,'./plots/DTFT.png','Resolution',300)

figure(11); hold on; grid on;
plot(f0,H0,dotOptions{:});
plot(f1_zero,H1_zero, plotOptions{:})
h = shade(w(:),lowerBound(:),w(:),upperBound(:),'FillType',[1 2;2 1],'FillColor',[1 1 1]*0.2); % offset required for zero crossings
delete(h(1:2)) % remove upper and lower bound lines

legend({'DFT points $X(f_k)$','DTFT $X(f)$','DFT bounds'},'Location','northoutside','NumColumns',3,'interpreter','latex')
xlabel('Frequency (bins)','interpreter','latex')
ylabel('Magnitude (linear)','interpreter','latex')
xlim([0,N2])
yticks([0:.2:1, 1.1])
set(gcf,'Position',[100 100 400 200])
% exportgraphics(gcf,sprintf('./plots/boundedDTFall%d.png',firCase),'Resolution',300)
