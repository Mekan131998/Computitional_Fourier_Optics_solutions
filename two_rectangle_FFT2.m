%% Exercise 5.2
clear;
close all;

%% Aperture a, rectange with sides wx and wy
lambda=0.5*10^-6;           % illumination wavelength
L1=2e-3;          %in mm, side length
M=500;          % number of samples
dx1=L1/M;     

wy=0.1e-3;     %m
wx=0.05e-3;    %m

x1=-L1/2:dx1:L1/2-dx1; % x coor
y1=x1;              % y coor

x01=-0.3e-3;
y01=0.3e-3;

[X1, Y1]=meshgrid(x1, y1);
[X2, Y2]=meshgrid(x1+x01, y1+y01);

u1=rect(X1/(2*wx)).*rect(Y1/(2*wy))+rect(X2/(2*wx)).*rect(Y2/(2*wy));
I1=abs(u1).^2;

figure(1); clf;
imagesc(x1, y1, I1);   %image display
colormap("gray");
axis square;
axis xy;
xlabel('x (m)'); ylabel('y (m)');

%% FFT of this thing 
%Taking fourier transform: 
g0=fftshift(u1);
G0=fft2(g0)*dx1^2;
G=fftshift(G0);

%freq. coordinates: 
fx=-1/(2*dx1):1/L1:1/(2*dx1)-1/L1;  
fy=fx;

figure(2); clf; 
imagesc(fx, fy, abs(G))
colormap('hot')
shading interp
ylabel('fy (cyc/m) '); xlabel('fx (cyc/m)');




