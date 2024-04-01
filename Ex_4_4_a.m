clear;
close all;

L=0.5e-3;          %in mm, side length
M=500;          % number of samples
dx=L/M;     

wx=0.1e-3;     %m
wy=0.05e-3;    %m

x=-L/2:dx:L/2-dx; % x coor
y=x;              % y coor

[X, Y]=meshgrid(x, y);

c=rect(X/(2*wx)).*rect(Y/(2*wy));

figure(1); clf;
imagesc(x, y, c);   %image display
colormap('parula');
axis square;
axis xy
xlabel('x (m)'); ylabel('y (m)');

%Taking fourier transform: 
g0=fftshift(c);
G0=fft2(g0)*dx^2;
G=fftshift(G0);

%freq. coordinates: 
fx=-1/(2*dx):1/L:1/(2*dx)-1/L;  
fy=fx;

figure(1); clf; 
imagesc(fx, fy, abs(G))
camlight left; lighting phong
colormap('sky')
shading interp
ylabel('fy (cyc/m) '); xlabel('fx (cyc/m)');

lambda=0.633e-6;   %in m 
k=2*pi/lambda;     %in 1/m

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

z=1.8*10*wy^2/lambda;        

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  D_lobe=1.22*lambda/w
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
D_lobe=1.22*lambda*z/wy;

L_side=5*D_lobe;

dx1=L_side/M; %sample interval 
x1=-L_side/2:dx1:L_side/2-dx1; y1=x1; %coords 
[X1,Y1]=meshgrid(x1,y1); 


U2=exp(1j*k*z)/(1j*lambda*z)*exp(1j*k/(2*z)*(X1.^2+Y1.^2)).*abs(G);

% Irradiance:

I2=abs(U2).^2;

figure(2); clf;
imagesc(x1, y1, nthroot(I2, 3));
xlabel('x (m)');
ylabel('y (m)');
title("Numeric")
colormap("sky");
axis square;
axis xy;

%% Analytic solution:

I_anal=((4*wx*wy/(lambda*z))^2)*(sinc((2*wx/(lambda*z))*X1).*sinc((2*wy/(lambda*z))*Y1))^2;
I_anal=abs(I_anal);

figure(3); clf;
imagesc(x1, y1, nthroot(I_anal, 3));
xlabel('x (m)');
ylabel('y (m)');
title("Analog solution")
colormap("sky");
axis xy;


%profile:
Ix=sum(I_anal, 1);
Iy=sum(I_anal, 2);

figure(4); clf;
subplot(2, 1, 1)
plot(x1, Iy);
title("Intensity profile y axis")

subplot(2, 1, 2)
plot(y1, Ix)
title("Intensity profile along x axis")





