% STC simulations

rw_width=5;
M=1000;
drw=rw_width/M;

rw=-rw_width/2:drw:rw_width/2;    % r/w

T=8;

t=-T/2:T/M:T/2;                 % t/t0

[RW, T]=meshgrid(rw, t);

omt=5:10/1000:15;

[RW, OMT]=meshgrid(rw, omt);

E=exp(-RW.^2).*exp(-T.^2).*exp(1j*OMT);

I=abs(E).^2;

figure(1); clf;
imagesc(t, rw, I);
colormap("parula");
axis square;


% AD/PFT

gamma=0.000001;     % magnitude of the AD/PFT
w0=10;

phi_w=gamma.*RW.*(OMT-w0);

figure(2); clf;
imagesc(omt, rw, phi_w);
colormap("parula");
axis square;

