% Seidel polynomials simulation 
L1=4;                               % in sm 
M=1000; 
dx1=L1/M;

x1=-L1/2:dx1:L1/2-dx1; y1=x1;          % source coords

[X1, Y1]=meshgrid(x1, y1);

u0=1;
v0=0;
wd=0;
w040=1;
w131=0;
w222=0;
w220=0;
w331=0;


w=seidel_5(u0, v0, X1, Y1, wd, w040, w131, w222, w220, w331);

P=circ(sqrt(X1.^2+Y1.^2));
mask=(P==0);
w(mask)=NaN;

figure(1); clf;
subplot(2, 2, 1)
surfc(x1, y1, w );
camlight left; lighting phong; 
colormap("parula"); shading interp; 
xlabel('x'); ylabel('y');

subplot(2, 2, 2)
surfc(x1, y1, seidel_5(1, 0, X1, Y1, 0, 1, 0, 1, 0, 0));
camlight left; lighting phong; 
colormap("parula"); shading interp; 
xlabel('x'); ylabel('y');

subplot(2, 2, 3)
surfc(x1, y1, seidel_5(0, 1, X1, Y1, 1, 1, 1, 0, 0, 0));
camlight left; lighting phong; 
colormap("parula"); shading interp; 
xlabel('x'); ylabel('y');

subplot(2, 2, 4)
surfc(x1, y1, seidel_5(1, 0, X1, Y1, 0, 0, 0, 0, 1, 0));
camlight left; lighting phong; 
colormap("parula"); shading interp; 
xlabel('x'); ylabel('y');
