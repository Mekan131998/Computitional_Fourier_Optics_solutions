%Fraunhofer irradiance plot
L=0.2;
M=250;

dx=L/M;
x=-L/2:dx:L/2-dx; %x coor
y=x;              %y coor
[X, Y]=meshgrid(x, y);

w=1e-3;           % x half-width
lambda=0.633e-6;  % wavelength
z=50;             % prop distance
k=2*pi/lambda; 
lz=lambda*z;

%irradiance
I2=(w^2/lz)^2.*(jinc(w/lz*sqrt(X.^2+Y.^2))).^2;

figure(1); clf; %irradiance image
imagesc(x, y, nthroot(I2,3));
xlabel('x (m)'); ylabel('y (m)');
colormap('turbo');
axis square;
axis xy; 

figure(2); clf; 
plot(x, I2(M/2+1, :));
xlabel('x (m)'); ylabel('Irradiance');




