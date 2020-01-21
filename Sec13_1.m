clear all; close all; clc

L=20; % define the computational domain [-L/2,L/2] 
n=128; % define the number of Fourier modes 2^n
x2=linspace(-L/2,L/2,n+1); % define the domain discretization 
x=x2(1:n); % consider only the first n points: periodicity

u=exp(-x.*x); % function to take a derivative of 
ut=fft(u); % FFT the function 
utshift=fftshift(ut); % shift FFT

figure(1), plot(x,u) % plot initial gaussian 
figure(2), plot(abs(ut)) % plot unshifted transform 
figure(3), plot(abs(utshift)) % plot shifted transform

%% TAKING DERIVATIVES

clear all; % clear all variables and figures
L=20; % define the computational domain [-L/2,L/2] 
n=128; % define the number of Fourier modes 2^n
x2=linspace(-L/2,L/2,n+1); % define the domain discretization 
x=x2(1:n); % consider only the first n points: periodicity 
dx=x(2)-x(1); % dx value needed for finite difference 

u=sech(x); % function to take a derivative of
ut=fft(u); % FFT the function
k=(2*pi/L)*[0:(n/2-1) (-n/2):-1]; % k rescaled to 2pi domain

% FFT calculation of derivatives
ut1=i*k.*ut; % first derivative
ut2=-k.*k.*ut; % second derivative 
u1=real(ifft(ut1)); u2=real(ifft(ut2)); % inverse transform 
u1exact=-sech(x).*tanh(x); % analytic first derivative 
u2exact=sech(x)-2*sech(x).^3; % analytic second derivative

% Finite difference calculation of first derivative
% 2nd-order accurate 
ux(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx); 
for j=2:n-1
    ux(j)=(u(j+1)-u(j-1))/(2*dx); 
end
ux(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx);

% 4th-order accurate 
ux2(1)=(-3*u(1)+4*u(2)-u(3))/(2*dx); 
ux2(2)=(-3*u(2)+4*u(3)-u(4))/(2*dx); 
for j=3:n-2
    ux2(j)=(-u(j+2)+8*u(j+1)-8*u(j-1)+u(j-2))/(12*dx); 
end
ux2(n-1)=(3*u(n-1)-4*u(n-2)+u(n-3))/(2*dx); 
ux2(n)=(3*u(n)-4*u(n-1)+u(n-2))/(2*dx);

figure(4), plot(x,u,'r',x,u1,'g',x,u1exact,'go',x,u2,'b',x,u2exact,'bo') 
figure(5)
subplot(3,1,1), plot(x,u1exact,'ks-',x,u1,'k',x,ux,'ko',x,ux2,'k*'), axis([-1.15 -0.75 0.47 0.5])
subplot(3,1,2), plot(x,u1exact,'ks-',x,u1,'kv',x,ux,'ko',x,ux2,'k*'), axis([-0.9376 -0.9374 0.49848 0.49850])
subplot(3,1,3), plot(x,u1exact,'ks-',x,u1,'kv',x,ux,'ko',x,ux2,'k*'), axis([-0.9376 -0.9374 0.498487 0.498488])
