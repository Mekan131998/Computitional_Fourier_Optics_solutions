% grating_cos diffraction grating example 

lambda=0.5e-6;          % wavelength
f=0.5;                  % propagation distance 
P=1e-4;                 % grating period 
D1=1.02e-3;             % grating side length 

L1=1e-2;                % array side length 
M=500;                  % # samples 
dx1=L1/M; 
x1=-L1/2:dx1:L1/2-dx1; y1=x1;  % source coords 
[X1, Y1]=meshgrid(x1, y1);

% graing field and irradiance: 
u1=1/2*(1-cos(2*pi*X1/P)).*rect(X1/D1).*rect(Y1/D1);
I1=abs(u1).^2;

figure(1); clf; 
imagesc(x1, y1, I1);
colormap("parula"); 

% Fraunhofer pattern 
[u2, L2]=propFF(u1, L1, lambda, f);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2;  y2=x2; 

I2=abs(u2).^2;

figure(2); clf; 
imagesc(x2, y2, I2);
colormap("parula"); 

