% Exercise 5.3 

clear; 
close all; 

% A useful diagnostic for propagation simulations is to compute the power
% in the source and observation planes. Assuming no asorbtion or scatter of
% the light, which is true for the simulations presented in this book the
% power should be conserved. In other words, the source and observation
% planes should contain the same optical power. if not, there may be a code
% error or a sampling problem
%        ____   _____         
%        \   '  \    '
%    P=   \      \    I(x, y) dx dy 
%         /      /
%        /____' /____'
%      

%% Source simulation: Square beam

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


%% Initial power 

P1=sum(I1, 1)*sum(I1, 2)*L1^2;
disp(P1);

I2=abs(propTF(u1, L1, lambda, z)).^2;
P2=sum(I2, 1)*sum(I2, 2)*L1^2;
disp(P2);
clear P2;

%%
z=[1000, 2000, 4000, 20000];

disp("IR power")
for i=1:4
    disp(z(i))
    I2=abs(propTF(u1, L1, lambda, z(i))).^2;
    P2=sum(I2, 1)*sum(I2, 2)*L1^2;
    disp(P2)
end

disp("IR power")
for i=1:4
    disp(z(i))
    I2=abs(propTF(u1, L1, lambda, z(i))).^2;
    P2=sum(I2, 1)*sum(I2, 2)*L1^2;
    disp(P2)
end
