%% This script is the same as evaluate_nicer.m but with more comments and some optional staff were removed:
% last update 05.04.2024
% 1. First we read the data:
% Int_xyt - intensity, 
% t- time, 
% x- physical x dimension of the camera:

clear;
clc;
close all;

filename=['C:\UserProgram\MATLAB\Mekan\20240320_ATTOlab\Test1_lowpower_short.h5'];

Int_xyt=h5read(filename,'/Int');
t=h5read(filename,'/t');
x=h5read(filename,'/x');

%% 2. introducing the time mask 
% define left t_left and t_right edges for the timemask 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter the left and right edge of time (fs):

t_left=-100;                        
t_right=100;                        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=size(Int_xyt);

%timemask is better when it is exponential smooth 
time_mask=single(permute(1./(1+exp((t-t_right)/3))./(1+exp((t_left-t)/3)),[1,3,2])); %time mask 

pix_x_cor=s(1)/2;                      % central pixel coordinates 
pix_y_cor=s(2)/2;

% plotting the 
figure(1);clf;
set(gcf,'color','w');
subplot(1,2,2)
test = (squeeze(Int_xyt(pix_x_cor,pix_y_cor,:))-mean(squeeze(Int_xyt(pix_x_cor,pix_y_cor,:)),[1,2,3])).*squeeze(time_mask);
plot(t,test)
xlabel("Time (fs)")
title("Windowed data")
subplot(1,2,1)
test = squeeze(Int_xyt(pix_x_cor,pix_y_cor,:)-mean(squeeze(Int_xyt(pix_x_cor,pix_y_cor,:)),[1,2,3]));
plot(t,test)
xlabel("Time (fs)")
title("Raw data")

%% (optional) see the camera image at the given time delay

% Enter the index of time delay (between 1-512):
idx_t=250;                                  % <----------------------------

figure(1);clf;
imagesc(Int_xyt(:, :, idx_t));
xlabel('x, pixel')
ylabel('y, pixel')
title(["Delay time: ", num2str(round(t(idx_t))), " fs"]);
axis square;
%% 3. Taking FFT 
% define x and y axis for the FFT
c=0.299792458;
siz=size(Int_xyt);          % size of Intxy_t
dx=mean(diff(x));           
x=4*[0:(siz(1)-1)]*dx;      % Here it is fixed 
y=x';                       % since image square y and x are the same
dt=mean(diff(t));           % time step

% Enter the resolution:
N=2^13;                     % <------- Here you can  change the resolution

T=dt*(1:N);
F=(0:N-1)/N/dt;

% Enter the upper and lower limit for the wavelength in um 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength_lower = 0.72;                          
wavelength_upper =  0.87;         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converting these limits into frequency domain 
f_upper = c/wavelength_lower;
f_lower = c/wavelength_upper;

%creating id to filter the data in a given range
id=permute((F>f_lower)&(F<f_upper),[1,3,2]);            % indexes of the given frequency range

K=sum(id);                                              % How many spectral points correspond that area
M=siz(1);                                               % number of points spatial 

Int0=mean(Int_xyt,3);                % Avg. intensity for each pixel, since the intensity in the third dimension 

Int_xyf=complex(zeros(M,M,K,'single')); %same dimension, x,y,f creating empty intensity data for frequency domain


% start FFT for each pixel: 
%%
for i=1:M
    disp(i)
    sp=fft(((Int_xyt(i,:,:))-(Int0(i,:))).*time_mask,N,3);      % slice by slice fft
    Int_xyf(i,:,:)=sp(1,:,id);                                  % save it with frequency filtering.
end

%% 4. Plotting the spectrum of the central point 

c=0.299792458;
f0=F(id);                   % selecting only filtered frequency values 
k0=2*pi*f0/(c);             % wavenumber of the filtered frequency values

% Here we find the spectral amplitude at the center: amp_c      (1x1x388)
amp_c=sqrt(abs(Int_xyf(M/2,M/2,:))).*exp(1i*angle(Int_xyf(M/2,M/2,:))); % spectral amplitude at center

% finding indexes that correspond to maximum and minimum values:
[maxval,maxind] = max(abs(squeeze(amp_c)/max(amp_c)));
w_central = c./f0(maxind);                              % central wavelength 

spectrum_intensity=abs(squeeze(amp_c))/max(abs(squeeze(amp_c)));

% Plots spectrum in center
figure(2); clf;
plot(f0,spectrum_intensity); 
xline(f0(maxind), '--r')
xlabel('frequency, PHz')
ylabel('spectral amplitude norm. unit')
title("Central wavelength: " + round(w_central*1e3) + " nm"+', ('+f0(maxind)+"PHz)")
axis square;

% the same thing but wavelength in the horizontal axis:

w0=c./f0;
figure(); clf; 
plot(w0, spectrum_intensity);
xlabel("wavelength, nm");
ylabel("spectrum intensity, norm. unit")
title("Central wavelength: " + round(w_central*1e3) + " nm")
axis square;

% Saving the found spectrum in txt file:
T=table(f0', spectrum_intensity);
writetable(T,'Spectrum_Test1_lowpower_short.txt');

% figure(); clf;
% set(gcf,'color','w');
% 
% subplot(1,2,1)
% plot(f0,abs(squeeze(amp_c)/max(amp_c)).^2);  % 3rd dimension stuff no work in matlab, equiv. to .ravel()
% title("Main wavelength: " + round(w_central*1e3) + " nm")
% xlabel("Frequency (PHz)")
% 
% subplot(1,2,2)
% sig= abs(squeeze(amp_c)).^2 .*f0'.^2;
% 
% plot(c./f0,sig/max(sig));  % 3rd dimension stuff no work in matlab, equiv. to .ravel()
% 
% name = 'MPC1_gdd_opt_1ps_camera_close_spectrum_100p_15_HR4D09681_17-53-20-462.txt';
% data = table2array(readtable(name,'Delimiter','\t', 'HeaderLines',14, 'DecimalSeparator', ','));
% wl_meas = data(:,1);
% spe = data(:,2);
% hold on
% plot(wl_meas*1e-3,spe/max(spe))
% legend("Reconstructed", "Measured")
% 
% xlim([0.980,1.080])
% ylim([0,1.1])
% 
% title("Main wavelength: " + round(w_central*1e3) + " nm")
% xlabel("Wavelength (um)")
%% 5. Subtraction of the 'known' reference beam
% Dimensions um (micrometer), PHz (petaherz), fs-(femtosecond)

% Normalize the spectral intensity relative to the value at the centre of the image: 

amp_xyf=Int_xyf./amp_c;

x=x-mean(x);                            % centering
y=x';   

% we need to find to parameters in order to subtract correct reference
% spherical beam from the interferometric image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the distance between camera and the focus and offset:

D=-3.816e5;                    

x0=100;
y0=-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(x,y);            % create a grid for the spherical wavefront, here x and y in um (real world values)

Sph=single(exp(1i*permute(k0,[1,3,2]).*((X-x0).^2+(Y-y0).^2)/2/D)); % Formula for the spherical wavefront

amp_xyf_R=amp_xyf.*Sph;         % HERE WE SUBTRACT THE SPHERICAL BEAM FROM THE INTERFEROMETRIC IMAGE: 

% plot the subtracted wavefront phase 
% in that image there should not be any fringes
% if there is a finges try to change values of D until you see image that
% contains no fringes. By changing x0 and y0 one can center the image: 

figure(5);clf;
set(gcf,'color','w');
imagesc(x/1e4, y/1e4, angle(amp_xyf_R(:,:,maxind)));                % divived 1e4 for converting into sm
xlabel('x size, sm');
ylabel('y size, sm');
axis square;
title("Here is the subtracted beam wavefront")
clear 'Sph';

% This is this lines of code for radial filtering of the image: 

R=single(sqrt(X.^2+Y.^2)); 

% Here you enter the mask radius: 
% mask radius and smooth parameter in um:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_mask=3500;                
smooth_param=500;            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is equation for mask: 
mask=1./(1+exp((R-r_mask)/smooth_param)); 

% one can also use the frequency mask
% equation for the frequency mask

mask_f=single(permute(1./(1+exp((f0-f_upper)/0.005))./(1+exp((f_lower-f0)/0.005)),[1,3,2])); 

% Here we apply both radial and frequency mask: 

amp_xyf_M=amp_xyf_R.*mask.*mask_f;
amp_xyf_M_before_filter = amp_xyf_M;  

kx=fftshift(2*pi/dx.*((-M/2+1):M/2)/M); 
imax=fft2(amp_xyf_M(:,:,maxind));
[Kx,Ky]=meshgrid(kx,kx);

% Here we plot for comparison:
figure(3);clf;
set(gcf,'color','w');

subplot(1,3,1)
imagesc(x/1e4, y/1e4, angle(amp_xyf(:,:,maxind)))         % phase of the interferogram, 
xlabel('x size, sm');   
ylabel('y size, sm');
axis xy
axis square
title("Farfield Phase")

subplot(1,3,2)
imagesc(x/1e4, y/1e4, angle(amp_xyf_R(:,:,maxind)))       % phase after subtraction of ref beam
xlabel('x size, sm');
ylabel('y size, sm');
axis xy
axis square
title("Corrected Farfield Phase")

subplot(1,3,3)
imagesc(x/1e4, y/1e4, angle(amp_xyf_R(:,:,maxind)).*mask)    % and phase when we apply frequency filtering
xlabel('x size, sm');
ylabel('y size, sm');
axis xy
axis square
title("With Spatial Filter")

%% 
figure(); clf;
imagesc(x/1e4, y/1e4, abs(amp_xyf_M(:,:,maxind)).*mask)    % wavefront after the subtraction
xlabel('x size, sm');
ylabel('y size, sm');
axis xy
axis square
title("Wavefront after the subtraction")

%% Plotting the wavefront after the reference beam subtraction:

figure(4); clf; 
imagesc(x/1e4, y/1e4, abs(amp_xyf_R(:,:,maxind))) 
xlabel('x size, sm')
ylabel('y size, sm')
axis xy
axis square
title("Wavefront after ref. beam subtraction")

%% 6. Filtering in the k space and back again 
Kr=sqrt(Kx.^2+Ky.^2);            % radial kl

% Enter the value for the divergency window for spatial filtering in rad 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

divg=0.0045;                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:K
    im0=fft2(amp_xyf_M(:,:,i))./(1+exp(single((Kr/k0(i))-divg)/0.0005)); % going to k space with divergency mask
    amp_xyf_M(:,:,i)=ifft2(im0);                                         % transform back with divergency applied
    disp(i);
end

sp=squeeze(amp_xyf_M(s(1)/2,s(2)/2,:));
sp0=sp/max(abs(sp));

%% Here we plot the comparison spectral image with maximum intensity before and after lowpass filter
% Here is important you select the maxind
after_filter=abs(amp_xyf_M(:,:,maxind)).^2;

figure(8);clf;
set(gcf,'color','w');
subplot(1,2,1)
imagesc(x/1e4, y/1e4, abs(amp_xyf_M_before_filter(:,:,maxind)).^2)
xlabel('x size, sm')
ylabel('y size, sm')
axis square
title("Before Lowpass filter")

subplot(1,2,2)
imagesc(x/1e4, y/1e4, after_filter )
xlabel('x size, sm')
ylabel('y size, sm')
axis square
title("After Lowpass filter")



% (optional) Here one can plot the spectral image for other points:

% figure(9);clf;
% set(gcf,'color','w');
% subplot(1,2,2)
% imagesc(abs(amp_xyf_M(:,:,50)).^2)
% title("After Lowpass filter loc 50 ")
% subplot(1,2,1)
% imagesc(abs(amp_xyf_M_before_filter(:,:,50)).^2)
% title("Before Lowpass filter loc 50")



%% Zernike Polynomials! 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radius in um, of the area where we are doing the zernike decomposition

rho=25000;                               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ang=atan2(Y,X);
siz=size(X);
Z=zeros(siz(1),siz(2),14);
r=R/rho;
id=r<1;
id=id/sqrt(sum(id(:)));

Z(:,:,1)=id;                     % piston
Z(:,:,2)=id.*(2.*(r).*sin(ang)); % tilt y
Z(:,:,3)=id.*(2.*(r).*cos(ang)); % tilt x
Z(:,:,4)=id.*(sqrt(6)*r.^2.*sin(2*ang)); % oblique asigmatism
Z(:,:,5)=id.*(sqrt(3)*(2*r.^2-1));       % defocus
Z(:,:,6)=id.*(sqrt(6)*r.^2.*cos(2*ang)); % vertical astimatism
Z(:,:,7)=id.*(sqrt(8)*r.^3.*sin(3*ang)); % Vertical trefoil
Z(:,:,8)=id.*(sqrt(8)*(3*r.^3-2*r).*sin(ang)); % Vertical coma
Z(:,:,9)=id.*(sqrt(8)*(3*r.^3-2*r).*cos(ang)); % Horizontal coma
Z(:,:,10)=id.*(sqrt(8)*r.^3.*cos(3*ang));      % Oblique trefoil
Z(:,:,11)=id.*(sqrt(10)*r.^4.*sin(4*ang));     % Oblique quadrafoil
Z(:,:,12)=id.*(sqrt(10)*(4*r.^4-3*r.^2).*sin(2*ang)); % Oblique secondary astigmatism
Z(:,:,13)=id.*(sqrt(5)*(6*r.^4-6*r.^2+1));             % Primary spherical
Z(:,:,14)=id.*(sqrt(10)*(4*r.^4-3*r.^2).*cos(2*ang));  % 	Vertical secondary astigmatism
Z(:,:,15)=id.*(sqrt(10)*r.^4.*cos(4*ang));             % 	Vertical quadrafoil

%%
co=zeros(15,K);
for i=1:K
    wf=angle(amp_xyf_M(:,:,i))./k0(i);
    Co=squeeze(sum(sum(wf.*Z,1),2));
    co(:,i)=Co;
    disp(i)
end

%% 
figure(185); clf; 

for i=1:15
    subplot(5, 3, i)
    hold on
    plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
    yyaxis right
    plot(f0,co(i,:), 'DisplayName',num2str(i))
    title(num2str(i))
    xlabel("f (Phz)")
    box on
    legend 
end 

%% Plotting phases for the different wavelength:

figure(186); clf; 
dC=round(K/20);

for i=dC:dC:K+dC
    wf=angle(amp_xyf_M(:, :, i))./k0(i);
    subplot(4, 5, round(i)/dC)
    imagesc(wf);
    axis square;
    title("Wavelength: "+w0(i)*1e3+'nm')
end

%% Plotting beam profile for different wavelength: 

dC=round(K/10);

figure(187); clf; 
for i=dC:dC:K+dC
    wf=abs(amp_xyf_M(:, :, i))./k0(i);
    subplot(2, 5, round(i)/dC)
    imagesc(wf);
    axis square;
    title("Beam profile: "+w0(i)*1e3+'nm')
end

%% Plotting the astigmatism: 
figure(4);clf;
set(gcf,'color','w');

subplot(2,2,1)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,co(4,:), 'DisplayName','oblique astigmatism')
title("Oblique Astigmatism")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,2)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,co(6,:), 'DisplayName','vertical astigmatism')
title("Vertical Astigmatism")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,3)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,co(5,:), 'DisplayName','defocus')
title("Defocus")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,4)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,co(14,:), 'DisplayName','Vertical secondary astigmatism')
title("Vertical secondary astigmatism")
xlabel("f (Phz)")
box on
legend 
%% Using zernike_coeffs.m function: 
addpath('C:\UserProgram\MATLAB\Mekan\STC\');

a_i=zeros(15,K);
for i=1:K
    wf=angle(amp_xyf_M(:,:,i))./k0(i);
    a_i(:,i)=zernike_coeffs(angle(amp_xyf_M(:,:,i))./k0(i), 15);
    disp(i)
end

%% Plotting coefficients for different frequencies:
figure(199); clf; 

for i=1:15
    subplot(5, 3, i)
    hold on
    plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
    yyaxis right
    plot(f0,a_i(i,:), 'DisplayName',num2str(i))
    title(num2str(i))
    xlabel("f (Phz)")
    box on
    legend 
end 


%% Doing correction test 
ver_astig_cor=-0.83*id.*(sqrt(6)*r.^2.*cos(2*ang));     % vertical atigmatism correction
test_ver=amp_xyf_M;
% applying correction:

for j=1:K
    disp(j)
    test_ver(:,:,j) = test_ver(:,:,j).*exp(1i*ver_astig_cor);
end

% finding zernike coefficients after the correction:
cor_coef=zeros(15,K);

for i=1:K
    disp(i)
    wf=angle(test_ver(:,:,i))./k0(i);
    Co=squeeze(sum(sum(wf.*Z,1),2));
    cor_coef(:,i)=Co;
end

% Plotting Zernike coefficients after the correction:
figure(5);clf;
set(gcf,'color','w');

subplot(2,2,1)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,cor_coef(4,:), 'DisplayName','oblique astigmatism')
title("Oblique Astigmatism")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,2)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,cor_coef(6,:), 'DisplayName','vertical astigmatism')
title("Vertical Astigmatism")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,3)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,cor_coef(5,:), 'DisplayName','defocus')
title("Defocus")
xlabel("f (Phz)")
box on
legend 

subplot(2,2,4)
hold on
plot(f0,abs(sp0).^2, '--k', 'DisplayName','spectrum')
yyaxis right
plot(f0,cor_coef(14,:), 'DisplayName','Vertical secondary astigmatism')
title("Vertical secondary astigmatism")
xlabel("f (Phz)")
box on
legend 
 
%% Plotting farfield spectral amplitude in k space:
% We define resolution:
P=2^14;             % <------ enter the resolution

% Here we define the grid with better resolution:
x_f=(single(linspace(min(x),max(x),P)));
y_f=(single(linspace(min(y),max(y),P)));
[X,Y]=(meshgrid(x_f,y_f));                      % new meshgrid


% Here is our old meshgrid: 
[X0,Y0]=meshgrid(x,y); 

id=abs(x_f)<100;            % select spatial part less than 100 micron
x2=x_f(id);                 % x, filtered #will_be_saved

dx=mean(diff(x_f));

% Here we define grid for the k space with the same resolution as new X and
% Y grid: 

k=single([0:(P/2),(-P/2+1):-1]/P/dx*2*pi);
[Kx,Ky]=meshgrid(k,k);
K2=(Kx.^2+Ky.^2); 
clear 'Kx' 'Ky'

% k space
kxf = fft2(amp_xyf_M(:,M/2,:));
kyf = fft2(amp_xyf_M(M/2,:,:));

figure(6);clf;
set(gcf,'color','w');

subplot(1,2,1)
imagesc(f0,fftshift(kx),(abs(fftshift(squeeze(kxf),1)).^2))
title("farfield, kx, f")
xlabel("f (PHz)")
ylabel("k (1/um)")

ylim([-0.02,0.02])
subplot(1,2,2)
imagesc(f0,fftshift(kx),(abs(fftshift(squeeze(kyf),1)).^2))
title("farfield, ky, f")
xlabel("f (PHz)")
ylabel("k (1/um)")
ylim([-0.02,0.02])


%% THE PROPAGATION!!! 
Dz=(single(100e3)); %100e3 micro meter focal length

amp_foc=complex(zeros(sum(id),sum(id),K,'single')); %focus profile spatial and spectrally
amp_foc_ideal=complex(zeros(sum(id),sum(id),K,'single')); %same "thing <3 "
diff_offset=250;

% i=325;
for i=1:K
    amp=interp2(X0,Y0,(amp_xyf_M(:,:,i)),X,Y,'linear',0).*exp(1i*k0(i)*(sqrt(Dz^2-X.^2-Y.^2)-Dz)); %we curve the wavefront with specified focal length
    amp_r=ifft2(fft2(amp).*exp(1i*sqrt(k0(i)^2-K2)*(Dz-diff_offset))); % diffraction gives offset of real focus
    amp_foc(:,:,i)=amp_r(id,id);
    amp_ideal=interp2(X0,Y0,abs(amp_xyf_M(:,:,i)),X,Y,'linear',0).*exp(1i*k0(i)*(sqrt(Dz^2-X.^2-Y.^2)-Dz)); %removes wavefront error from each frequency (abs)
    amp_r=ifft2(fft2(amp_ideal).*exp(1i*sqrt(k0(i)^2-K2)*(Dz-250)));
    amp_foc_ideal(:,:,i)=amp_r(id,id);
    disp(i);
end

%% Reconstructed image at the focus!

% Define the desired figure width and height
figureWidth = 1000; % in pixels
figureHeight = 400; % in pixels

% Create a new figure with the specified width and height
figure('Position', [100, 100, figureWidth, figureHeight]);

selected_freq_index=maxind;
set(gcf,'color','w');

subplot(1, 3, 1)
imagesc(abs(amp_foc(:,:,selected_freq_index)).^2)
title("Reconstructed image at the focus!")

siz=size(amp_foc);
q=exp(-1i*angle(amp_foc(round(siz(1)/2),round(siz(2)/2),:)));
amp_x=amp_foc.*q; %remove the phase in the middle, mulitply: -i, divide: +i
axis square 

subplot(1, 3,2)
imagesc(abs(amp_x(:,:,selected_freq_index)).^2)
title("after removing phase in the middle")

q=exp(-1i*angle(amp_foc_ideal(round(siz(1)/2),round(siz(1)/2),:))); % same "thing <3 " for ideal one
amp_x_ideal=amp_foc_ideal.*q;
axis square

subplot(1, 3, 3)
imagesc(abs(amp_x_ideal(:,:,selected_freq_index)).^2)
title("the same but for the ideal case")
axis square
%% Converting to the time domain: 

P=2^13;
df=mean(diff(f0));
t0=linspace(-0.5/df,0.5/df,P);
id=abs(t0)<150; %everthing plus minus 200 fs
N=sum(id); %how many points

amp_ft=(zeros(siz(1),siz(2),N,'single')); %spatial grid
amp_ft_ideal=(zeros(siz(1),siz(2),N,'single'));

% for i=1:siz(1)
%     for j=1:siz(1)
%         at=fftshift(abs(fft(amp_x(i,j,:),P,3)).^2);
%         amp_ft(i,j,:)=at(:,:,id);
%     end
% end
for i=1:siz(1)
    at=fftshift(abs(fft(amp_x(i,:,:),P,3)).^2,3);
    amp_ft(i,:,:)=at(:,:,id);
    at=fftshift(abs(fft(amp_x_ideal(i,:,:),P,3)).^2,3);
    amp_ft_ideal(i,:,:)=at(:,:,id);
    i
end
%%
int_xyf=sum(abs(amp_x).^2,3);
int_yf=squeeze(sum(abs(amp_x).^2,1));
int_xf=squeeze(sum(abs(amp_x).^2,2));

int_xyf=int_xyf/max(max(int_xyf));
int_x=int_xf/max(max(int_xf));
int_y=int_yf/max(max(int_yf));


%% Plot in space time 
figure(88);clf;
subplot(2,2,1)
imagesc(f0,x2,int_xt)
xlabel("f (PHz)")
ylabel("x (um)")

subplot(2,2,2)
imagesc(f0,x2,int_yt)
xlabel("f (PHz)")
ylabel("y (um)")

%% plot in space frequency 
figure(89);clf;
subplot(1,2,1)
imagesc(f0,x2,int_xf)
xlabel("f PHz")
ylabel("x (um)")

subplot(1,2,2)
imagesc(f0,x2,int_yf)
xlabel("f, PHz")
ylabel("y (um)")

%% Plotting in 3D 


figure(90);
isosurface(x2,x2,t0, amp_ft, 1e-4);
xlabel('x (\mum)');
ylabel('y (\mum)');
zlabel('t (fs)');


%% saving the stuff
% meta info should contain: divergence filter param, radial filter param,
% x0,y0, distance from focus to detector.
% x2: spatial axis in um, here: +- 200 um
% t0: time window selected in fs, here plus minus 200 fs
% f0: frequency window, pHz
% % amp_ft and amp_ft_ideal: 3D in xyt, x,y in terms of focus coordinates
% % amp_x and amp_x_ideal: 3D in xyf, x,y in terms of focus coordinates
t0 = t;
% folder = "C:\UserProgram\MATLAB\STC DESY\Results\Results with CORRECTED diagnostic code\";
folder = "C:\UserProgram\MATLAB\Mekan\STC\results\attolab_STC_results\";
name = "Attolab_lowPower_shortDuration";

h5create(folder + name + "_zernike.h5",'/f0',[size(f0)]);
h5create(folder + name + "_zernike.h5",'/co',[size(co)]);

h5write(folder + name + "_zernike.h5",'/f0',f0);
h5write(folder + name + "_zernike.h5",'/co',co);

h5create(folder + name + "_xyf.h5",'/amplitude_xyf_real',[size(amp_x)]);
h5create(folder + name + "_xyf.h5",'/amplitude_xyf_imag',[size(amp_x)]);
h5create(folder + name + "_xyf.h5",'/f',[size(f0)]);
h5create(folder + name + "_xyf.h5",'/x',[size(x2)]);

h5create(folder + name + "_xyt.h5",'/amplitude_xyt_real',[size(amp_ft)]);
h5create(folder + name + "_xyt.h5",'/amplitude_xyt_imag',[size(amp_ft)]);
h5create(folder + name + "_xyt.h5",'/t',[size(t0)]);
h5create(folder + name + "_xyt.h5",'/x',[size(x2)]);

h5create(folder + name + "_xyf_ideal.h5",'/amplitude_xyf_ideal_real',[size(amp_x_ideal)]);
h5create(folder + name + "_xyf_ideal.h5",'/amplitude_xyf_ideal_imag',[size(amp_x_ideal)]);
h5create(folder + name + "_xyf_ideal.h5",'/f',[size(f0)]);
h5create(folder + name + "_xyf_ideal.h5",'/x',[size(x2)]);

h5create(folder + name + "_xyt_ideal.h5",'/amplitude_xyt_ideal_real',[size(amp_ft_ideal)]);
h5create(folder + name + "_xyt_ideal.h5",'/amplitude_xyt_ideal_imag',[size(amp_ft_ideal)]);
h5create(folder + name + "_xyt_ideal.h5",'/t',[size(t0)]);
h5create(folder + name + "_xyt_ideal.h5",'/x',[size(x2)]);

h5write(folder + name + "_xyf.h5",'/amplitude_xyf_real',real(amp_x));
h5write(folder + name + "_xyf.h5",'/amplitude_xyf_imag',imag(amp_x));
h5write(folder + name + "_xyf.h5",'/f',f0);
h5write(folder + name + "_xyf.h5",'/x',x2);

h5write(folder + name + "_xyt.h5",'/amplitude_xyt_real',real(amp_ft));
h5write(folder + name + "_xyt.h5",'/amplitude_xyt_imag',imag(amp_ft));
h5write(folder + name + "_xyt.h5",'/t',t0);
h5write(folder + name + "_xyt.h5",'/x',x2);

h5write(folder + name + "_xyf_ideal.h5",'/amplitude_xyf_ideal_real',real(amp_x_ideal));
h5write(folder + name + "_xyf_ideal.h5",'/amplitude_xyf_ideal_imag',imag(amp_x_ideal));
h5write(folder + name + "_xyf_ideal.h5",'/f',f0);
h5write(folder + name + "_xyf_ideal.h5",'/x',x2);

h5write(folder + name + "_xyt_ideal.h5",'/amplitude_xyt_ideal_real',real(amp_ft_ideal));
h5write(folder + name + "_xyt_ideal.h5",'/amplitude_xyt_ideal_imag',imag(amp_ft_ideal));
h5write(folder + name + "_xyt_ideal.h5",'/t',t0);
h5write(folder + name + "_xyt_ideal.h5",'/x',x2);


%Optional staff:

%% Movie
load("colormap_whitejet.dat")
cm = dlmread('colormap_whitejet.dat');
colormap(cm)

videoFileName = 'testmovie.mp4';
writerObj = VideoWriter(videoFileName, 'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj);


figure(11);clf;
set(gcf,'color','w');
maximum = max(Int_xyt(:,:,:),[],'all') ;
%m = mean(Int_xyt(:,:,:),[1,2,3]);
movie_int = Int_xyt./maximum;
% Create an axes handle for the left subplot
% Specify the height ratio for the left and right subplots

for i = 1:512
    imagesc(movie_int(:,:,i))
    title(round(t(i))+" fs")
    axis equal
    axis xy
    xlim([1,1024]);
    ylim([1,1024]);
    colormap(cm)
    colorbar()
    caxis([0 1])
    xlabel("Pixel")
    ylabel("Pixel")

    drawnow
    frame = getframe(gcf);
    writeVideo(writerObj, frame)
end

close(writerObj);
videoPlayer = VideoReader(videoFileName);
%implay(videoPlayer);

