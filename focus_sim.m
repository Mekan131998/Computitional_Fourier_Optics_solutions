% Square beam example
close all; 
clear all;
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
%% Focus the beam:

zf=2000;
u2=focus(u1, L1, lambda, zf);
I2=abs(u2).^2;

% Source simulation plot
figure(2); clf; 
imagesc(x1, y1, I2);
axis square; axis xy;
colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
title("focused z=0 m");

%% progagate focusing beam:

u3_before=propTF(u2, L1, lambda, zf-800);
u3=propTF(u2, L1, lambda, zf);
u3_after=propTF(u2, L1, lambda, zf+800);

I3=abs(u3).^2;

figure(3); clf;

subplot(1, 3, 1)
imagesc(x1, y1, abs(u3_before).^2);
axis square; axis xy;
colormap("parula"); xlabel("x (m)"); ylabel("y (m)");
title(sprintf("before focux beam z=%d", zf-800))

subplot(1, 3, 2)
imagesc(x1, y1, I3);
axis square; axis xy;
colormap("parula"); xlabel("x (m)"); ylabel("y (m)");
title(sprintf("Focused beam zf=%d", zf))

subplot(1, 3, 3)
imagesc(x1, y1, abs(u3_after).^2);
axis square; axis xy;
colormap("parula"); xlabel("x (m)"); ylabel("y (m)");
title(sprintf("after focus z=%d", zf+800))

%% Trying negative focus:
zf=-2000;
u2_neg=focus(u1, L1, lambda, zf);

% progagate: 
z=3000;
u3_prop=propTF(u2_neg, L1, lambda, z);

I3_pr=abs(u3_prop).^2;

figure(4);
imagesc(x1, y1, I3_pr);
axis square; axis xy;
colormap("parula"); xlabel("x (m)"); ylabel("y (m)");
title(sprintf("Negative focus beam zf=%d", zf))
