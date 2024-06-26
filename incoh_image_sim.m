%% Incoherent Imaging code:
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

%% Incoherent imaging simulation:

lambda=0.5*10^-6;           % wavelength 
wxp=6.25e-3;                % exit pupil radius 
zxp=125e-3;                 % exit pupil distance 
f0=wxp/(lambda*zxp);        % coherent cutoff 

fu=-1/(2*du):1/L:1/(2*du)-(1/L);        % freq coords 
fv=fu; 

[Fu, Fv]=meshgrid(fu, fv);
H=circ(sqrt(Fu.^2+Fv.^2)/f0);
OTF=ifft2(abs(fft2(fftshift(H))).^2);
OTF=abs(OTF/OTF(1, 1));

figure(2); clf;         % checking OTF
imagesc(fu, fv, fftshift(OTF))
%camlight left; lighting phong; 
colormap('gray')
%shading interp
ylabel('fu (cyc/m)'); xlabel('fv (cyc/m)');

%%
Gg=fft2(fftshift(Ig));              % convolution 
Gi=Gg.*OTF; 
Ii=ifftshift(ifft2(Gi));

% remove residual imag parts, values < 0 
Ii=real(Ii); mask=Ii>=0; Ii=mask.*Ii;

figure(3); clf;                     % image result 
imagesc(u, v, nthroot(Ii,2));
colormap('gray'); xlabel('u (m)'); ylabel('v (m)');
axis square; 
axis xy; 

figure(4);                                  % horizontal image slice 
vvalue=0.2e-4;                              % select row (y value)
vindex=round(vvalue/du+(M/2+1));            % convert row index 
plot(u, Ii(vindex, :), u, Ig(vindex, :), ':');
xlabel('u (m)'); ylabel('Irradiance');





