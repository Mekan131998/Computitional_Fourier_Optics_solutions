%% Fraunhofer propagation simulation 
clear;
close all;
%% Define the source 
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
w=0.011;                    % source half width (m) 
z=2000;                     % propagation dist (m)

[X1, Y1]=meshgrid(x1, y1);
u1=rect(X1/(2*w)).*rect(Y1/w);      % src field 
I1=abs(u1.^2);     
w=0.011;            % source half width (m)


%% Fraunhofer propagator 
[u2, L2]=propFF(u1, L1, lambda, z);

dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2;      % observation plane coordinates
y2=x2;

I2=abs(u2)^2;

%% Plottting Irradiance and magnitude using Fraunhofer propagator
figure(1); clf;
imagesc(x2, y2, nthroot(I2, 3));
axis square; axis xy; 
colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
title(sprintf("Irradiance using Fraunhofer propagator, z=%d", z));

figure(2); clf; 
plot(x2, abs(u2(M/2+1, :)));
xlabel("x (m)"); ylabel(" Magnitude ");
title(sprintf('Magnitude Fraunhofer propagator z=%d',z));
axis square; 