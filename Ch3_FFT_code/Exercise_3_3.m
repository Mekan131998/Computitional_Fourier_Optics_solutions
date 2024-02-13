%Exercise 3.3 

radius=0.015;   %in m 
L=0.2;          %in m, side length
M=200;          % number of samples
dx=L/M;         

x=-L/2:dx:L/2-dx; % x coor
y=x;              % y coor

[X, Y]=meshgrid(x, y);

c=circle_func_0015(X, Y, radius);

figure(1); clf;
imagesc(x, y, c);   %image display
colormap('parula');
axis square;
axis xy
xlabel('x (m)'); ylabel('y (m)');

g0=fftshift(c);
G0=fft2(g0)*dx^2;
G=fftshift(G0);

fx=-1/(2*dx):1/L:1/(2*dx)-1/L;  %freq coordinates
fy=fx;

figure(2); clf; 
surf(fx, fy, abs(G))
camlight left; lighting phong
colormap('parula')
shading interp
ylabel('fy (cyc/m) '); xlabel('fx (cyc/m)');





