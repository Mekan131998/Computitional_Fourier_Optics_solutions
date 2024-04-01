clear;
close all;

L = 0.5e-3;          % in mm, side length
M = 500;             % number of samples
dx = L / M;

wx = 0.1e-3;         % m
wy = 0.05e-3;        % m

x = -L / 2:dx:L / 2 - dx; % x coor
y = x;                    % y coor

[X, Y] = meshgrid(x, y);

c = rect(X / (2 * wx)) .* rect(Y / (2 * wy));

figure(1); clf;
imagesc(x, y, c);   % image display
colormap('parula');
axis square;
axis xy
xlabel('x (m)'); ylabel('y (m)');

% Taking fourier transform:
g0 = fftshift(c);
G0 = fft2(g0) * dx^2;
G = fftshift(G0);

% freq. coordinates:
fx = -1 / (2 * dx):1 / L:1 / (2 * dx) - 1 / L;
fy = fx;

figure(2); clf;
imagesc(fx, fy, abs(G))
camlight left; lighting phong
colormap('sky')
shading interp
ylabel('fy (cyc/m) '); xlabel('fx (cyc/m)');

lambda = 0.633e-6;   % in m
k = 2 * pi / lambda; % in 1/m

z = 1.8 * 10 * wy^2 / lambda;

D_lobe = 1.22 * lambda * z / wy;

L_side = 5 * D_lobe;

dx = L_side / M; % sample interval
x = -L_side / 2:dx:L_side / 2 - dx; y = x; % coords
[X, Y] = meshgrid(x, y);

U2 = exp(1j * k * z) / (1j * lambda * z) * exp(1j * k / (2 * z) * (X.^2 + Y.^2)) .* G;

% Irradiance:

I2 = abs(U2).^2;

figure(3); clf;
imagesc(x, y, nthroot(I2, 3));
xlabel('x (m)');
ylabel('y (m)');
title("Numeric")
colormap("sky");
axis square;
axis xy;

%% Analytic solution:

I_anal = ((4 * wx * wy / (lambda * z))^2) * (sinc((2 * wx / (lambda * z)) * X) .* sinc((2 * wy / (lambda * z)) * Y)).^2;
I_anal = abs(I_anal);

figure(4); clf;
imagesc(x, y, nthroot(I_anal, 3));
xlabel('x (m)');
ylabel('y (m)');
title("Analog solution")
colormap("sky");
axis square;
axis xy;


% profile:
Ix = sum(I_anal, 1);
Iy = sum(I_anal, 2);

figure(5); clf;
subplot(2, 1, 1)
plot(x, Iy);
title("Intensity profile y axis")

subplot(2, 1, 2)
plot(x, Ix)
title("Intensity profile along x axis")
