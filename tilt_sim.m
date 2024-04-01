% Square beam example
close all; 
%% Source simulation: 

% Consider a source plane with dimensions 0.5 x 0.5 m (L1=0.5 m) 

L1=0.5;                     % side length
M=250;                      % number of samples
dx1=L1/M;                   % src sample interval 
x1=-L1/2:dx1:L1/2-dx1;      % src coords
y1=x1;

% assume a square aperture with a half width of 0.051 m (51 mm) illuminated
% by a unit-amplitude plane wave from the backside 

lambda=0.5*10^-6;           % illumination wavelength
k=2*pi/lambda;              % wavenumber
w=0.051;                    % source half width (m) 
z=2000;                     % propagation dist (m)

[X1, Y1]=meshgrid(x1, y1);
u1=rect(X1/(2*w)).*rect(Y1/w);      % src field 
I1=abs(u1.^2);                      % src irradiance 

% Source simulation plot
figure(1); clf; 
imagesc(x1, y1, I1);
axis square; axis xy;
colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
title("z=0 m");

%% simulate the tilt: 

alpha=5e-5;
theta=deg2rad(45);

u2=tilt(u1, L1, lambda, alpha, theta);

I2=abs(u2).^2;

% NF tilted beam simulation plot
figure(2); clf; 
imagesc(x1, y1, I2);
axis square; axis xy;
colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
title("tilted beam z=0 m");

%% propagate tilted beam:
z=2000;
u3=propIR(u2, L1, lambda, z);    % propagation 

x2=x1;                      % observation plane coordinates
y2=y1;
I3=abs(u3).^2;              % observation plane irradiation

% plotting irradiance
figure(3); clf; 
imagesc(x2, y2, I3);
axis square; axis xy; 
colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
title(["z=", num2str(z), " m"]);
