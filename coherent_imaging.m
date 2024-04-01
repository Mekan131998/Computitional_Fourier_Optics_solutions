%% Coherence imaging example 
clear;
close all;

%% simulating the ideal image:

A=imread('USAF-1951.gif');

[M, N ]=size(A);                    % get image sample rate 
A=flipud(A);                        % reverse row order 
Ig=single(A);                       % integer to floating 
Ig=Ig(1:362, 100:461, :);           % dirty way of cropping an image
Ig=Ig./max(max(Ig));                 % normalize ideal image 

ug=sqrt(Ig);                        % ideal image field 
L=0.3e-3;                          % image plane side length (m)
du=L/(M-1);                             % sample interval 
u=-L/2:du:L/2-du; v=u;

figure(1); clf;                     % check ideal image
imagesc(u, v, Ig);
colormap("gray"); xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy 

%% Generating coherent transfer function: 
lambda=0.5*10^-6;               % wavelength
wxp=6.25e-3;                    % exit pupil radius
zxp=125e-3;                     % exit pupil distance 
f0=wxp/(lambda*zxp);            % cutoff frequency 

fu=-1/(2*du):1/L:1/(2*du)-(1/L); fv=fu;

[Fu, Fv]=meshgrid(fu, fv);
H=circ(sqrt(Fu.^2+Fv.^2)./f0);

figure(2); clf; 
surf(fu, fv, H.*.99);
camlight left; lighting phong; 
colormap("parula")
shading interp 
ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');

%% Simulating the image:
H=fftshift(H);
Gg=fft2(fftshift(ug));
Gi=Gg.*H;
ui=ifftshift(ifft2(Gi));
Ii=(abs(ui)).^2;

figure(3); clf;
imagesc(u, v, nthroot(Ii, 2)); 
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy 

% figure(4); clf;
% vvalue=-0.8e-4; 
% vindex=

%% Cross section comparison of the image:
figure(4); clf;
vvalue=-1.3e-4;         % select row (y value)
vindex=round(vvalue/du+(M/2+1));      % convert row index 
plot(u, Ii(vindex, :), u, Ig(vindex, :), ':');
xlabel('u (m)'); ylabel('Irradiance');


%% Illuminating with real light 

% applying a random complex exponential phase term to the object (ideal
% image) field:

ug=sqrt(Ig).*exp(1j*2*pi*rand(M-1));          % ideal image 

Gg=fft2(fftshift(ug));
Gi=Gg.*H;
ui=ifftshift(ifft2(Gi));
Ii=(abs(ui)).^2;

figure(5); clf;                         % plotting the ideal image illumination
imagesc(u, v, nthroot(Ii, 2)); 
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square
axis xy 



