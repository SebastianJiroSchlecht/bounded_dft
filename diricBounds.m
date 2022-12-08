function [f0c,lx,ux] = diricBounds(N2)
% Compute the lower and upper bounds for a given DFT lengths

N = N2*2+1;
piN = pi/N; 

f0 = (-N2:N2).';
f0c = circshift(f0,N2+1); % centered f0, start at 0
Ueven = mod(f0c,2)==0;
Uodd = mod(f0c,2)==1;
Uzero = (f0c == -1 | f0c == 0);
UX = 0*f0c;
UX(Ueven) = -1 ./ ( sin(0)* cos(piN*f0c(Ueven)) - cos(0)*sin(piN*f0c(Ueven)) ); % evenneg / evenpos
UX(Uodd) = 1 ./ ( sin(-piN)* cos(piN*f0c(Uodd)) - cos(piN)*sin(piN*f0c(Uodd)) ); % oddneg / oddpos
UX(Uzero) = 0; % not a good bound for neighbors

LX = 0*f0c;
LX(Uodd) = 1 ./ ( sin(0)* cos(piN*f0c(Uodd)) - cos(0)*sin(piN*f0c(Uodd)) ); % evenneg / evenpos
LX(Ueven) = -1 ./ ( sin(-piN)* cos(piN*f0c(Ueven)) - cos(piN)*sin(piN*f0c(Ueven)) ); % oddneg / oddpos
LX(Uzero) = 0; % not a good bound for neighbors

ux = ifft(ifftshift(UX));
lx = ifft(ifftshift(LX));

