%% Exercise 5.2
clear;
close all;

%% Aperture a, rectange with sides wx and wy
lambda=0.5*10^-6;           % illumination wavelength
L1=2e-3;          %in mm, side length
M=500;          % number of samples
dx1=L1/M;     

wx=0.1e-3;     %m
wy=0.05e-3;    %m

x1=-L1/2:dx1:L1/2-dx1; % x coor
y1=x1;              % y coor

[X1, Y1]=meshgrid(x1, y1);

u1=rect(X1/(2*wx)).*rect(Y1/(2*wy));
I1=abs(u1).^2;

figure(1); clf;
imagesc(x1, y1, I1);   %image display
colormap("gray");
axis square;
axis xy;
xlabel('x (m)'); ylabel('y (m)');

%% b. Propagation simulation 
% Illustration of limits of the TF and IR propagation algorithms: 

x2=x1;                      % observation plane coordinates
y2=y1;

z=[0.5e-2, 0.5e-2, 1e-2, 1e-2, 5e-2, 5e-2];

% illustrating observation plane irradiances:
figure(9); clf;
tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:6
    if mod(i, 2)
        nexttile
        imagesc(x2, y2, abs(propTF(u1, L1, lambda, z(i)).^2)); 
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('TF z=%d',z(i)));
    else 
        nexttile
        imagesc(x2, y2, abs(propIR(u1, L1, lambda, z(i)).^2));
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('IR z=%d',z(i)));
    end
end

%% Aperture b) ring 
% Source simulation


lambda=0.5*10^-6;           % illumination wavelength
L1=2.5e-2;          %in m, side length
M=500;          % number of samples
dx1=L1/M;     

w1=1e-3;     %m
w2=0.2e-3;    %m

x1=-L1/2:dx1:L1/2-dx1; % x coor
y1=x1;              % y coor

[X1, Y1]=meshgrid(x1, y1);
R1=sqrt(X1.^2+Y1.^2);

u1=circ(R1/(w1))-circ(R1/w2);

I1=abs(u1).^2;

figure(1); clf;
imagesc(x1, y1, I1);   %image display
colormap('gray');
axis square;
axis xy
xlabel('x (m)'); ylabel('y (m)');

%% propagation
x2=x1;                      % observation plane coordinates
y2=y1;

z=[0.5, 0.5, 2, 2, 5, 5];

% illustrating observation plane irradiances:
figure(9); clf;
tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:6
    if mod(i, 2)
        nexttile
        imagesc(x2, y2, abs(propTF(u1, L1, lambda, z(i)).^2)); 
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('TF z=%d',z(i)));
    else 
        nexttile
        imagesc(x2, y2, abs(propIR(u1, L1, lambda, z(i)).^2));
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('IR z=%d',z(i)));
    end
end

%% Aperture c) Two circle located at s distance between each other 
% 
L1 = 2.5e-2;         % in m, side length
M = 500;          % number of samples
dx1 = L1 / M;     

w = 1e-3;     % m
s = 4e-3;     % m 
x1 = -(L1)/2:dx1:(L1)/2-dx1; % x coor
y1 = x1; % y coor

[X1, Y1] = meshgrid(x1, y1);

R = sqrt((X1).^2 + Y1.^2);
R1 = sqrt((X1-s).^2 + Y1.^2);
u1 = circ(R/w)+circ(R1/w);
I1=abs(u1).^2;

figure(1); clf;
imagesc(x1, y1, I1);   %image display
colormap('gray');
axis square;
axis xy;

%% propagation:
x2=x1;                      % observation plane coordinates
y2=y1;

z=[0.5, 0.5, 2, 2, 5, 5];

% illustrating observation plane irradiances:
figure(9); clf;
tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'compact');
for i=1:6
    if mod(i, 2)
        nexttile
        imagesc(x2, y2, abs(propTF(u1, L1, lambda, z(i)).^2)); 
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('TF z=%d',z(i)));
    else 
        nexttile
        imagesc(x2, y2, abs(propIR(u1, L1, lambda, z(i)).^2));
        colormap("gray"); xlabel("x (m)"); ylabel("y (m)");
        title(sprintf('IR z=%d',z(i)));
    end
end