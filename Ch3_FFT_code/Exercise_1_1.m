%Exercises on computitional Fourier optics
% 1.7 Exercises:
% 1.1. Scetch the following functions:
% Defining the rectangle function
% plotting
t = -8:0.1:8;

[X, Y] = meshgrid(t, t);

circle_values=circle_func(X-2, Y)+circle_func(X+2, Y);
ex_e=(ucomb(t/4).*triangle(t)).*rectan(t/12);

figure(1); clf; 
subplot(3, 2, 1)
plot(t, rectan(t))
title('Rectangle Function')

subplot(3, 2, 2)
plot(t, triangle(t))
title('Triangle Function')

subplot(3, 2, 3)
surf(X, Y, double(circle_values));

subplot(3, 2, 4)
plot(t, ex_e);


