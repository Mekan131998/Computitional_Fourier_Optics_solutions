% Exercise 4.4 c
%%
clear;
close all;

L = 12e-3;         % in m, side length
M = 1000;          % number of samples
dx = L / M;     

w = 1e-3;     % m
s = 4e-3;     % m 
x = -(L)/2:dx:(L)/2-dx; % x coor
y = x; % y coor

[X, Y] = meshgrid(x, y);

R = sqrt((X).^2 + Y.^2);
R1 = sqrt((X-s).^2 + Y.^2);
c = circ(R/w)+circ(R1/w);

figure(1); clf;
imagesc(x, y, c);   %image display
colormap('parula');
axis square;
axis xy;

%Taking fourier transform: 
g0=fftshift(c);
G0=fft2(g0)*dx^2;
G=fftshift(G0);

%freq. coordinates: 
fx=-1/(2*dx):1/L:1/(2*dx)-1/L;  
fy=fx;

figure(2); clf; 
imagesc(fx, fy, abs(G))
camlight left; lighting phong
colormap('parula')
shading interp
ylabel('fy (cyc/m) '); xlabel('fx (cyc/m)');

lambda=0.633e-6;   %in m 
k=2*pi/lambda;     %in 1/m
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fresnel number N_E=w^2/(lmbda*z) 
% where 
%       lambda -wavelength
%       w- size of the aperture 
%       z- propagation distance
%  condition for Fresnel zone: 
%          N_E<0.1 
%   so 
%           z>10*w^2/lambda
% % % % % % % % % % % % % % % % % % % % % % % 

z=1.8*10*w^2/lambda;        

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  D_lobe=1.22*lambda/w
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
D_lobe=1.22*lambda*z/w;

L_side=5*D_lobe;
M=1000; %# samples 
dx=L_side/M; %sample interval 
x=-L_side/2:dx:L_side/2-dx; y=x; %coords 
[X,Y]=meshgrid(x,y); 


U2=exp(1j*k*z)/(1j*lambda*z)*exp(1j*k/(2*z)*(X^2+Y^2)).*G;

% Irradiance:

I2=abs(U2)^2;

figure(3);clf;
imagesc(x, y, nthroot(I2, 3));
xlabel('x (m)');
ylabel('y (m)');
colormap("parula");
axis square;
axis xy;

Iy=sum(I2, 1)/max(sum(I2, 1));
Ix=sum(I2, 2)/max(sum(I2, 2));

figure();clf;
subplot(1, 2, 1)
plot(x, Iy);
% xlim(-0.01, 0.01);

subplot(1, 2, 2)
plot(y, Ix);
xlim(-0.01, 0.01);


%% Analog solution:
r=sqrt(X.^2+Y.^2);
I_anal=(2/(lambda*z))^2*(w^2*jinc(w*r/(lambda*z)))^2.*(cos(pi*s*X/(lambda*z)))^2;
I_anal=abs(I_anal);


figure(4); clf;
imagesc(x, y, nthroot(I_anal, 3));
xlabel('x (m)');
ylabel('y (m)');
title("Analog solution")
colormap("parula");
axis xy;


%profile:
Ix=sum(I_anal, 1);
Iy=sum(I_anal, 2);

figure(5); clf;
subplot(2, 1, 1)
plot(x, Iy);
title("Intensity profile y axis")

subplot(2, 1, 2)
plot(y, Ix)
title("Intensity profile along x axis")

