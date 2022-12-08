function c = circspace(N)
% function description
% N points as rad degree without repetition of pi
% 
% (c) Sebastian Jiro Schlecht:  16. October 2018

c = linspace(0,2*pi,N+1);
c = c(1:end-1);

c = mod(c+pi,2*pi)-pi;
c = fftshift(c);