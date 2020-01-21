clear all; close all; clc

%% NOISE ON SIGNAL

L=30; % time slot to transform 
n=512; % number of Fourier modes 2^9
t2=linspace(-L,L,n+1); t=t2(1:n); % time discretization 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; % frequency components of FFT
u=sech(t); % ideal signal in the time domain 
figure(1), subplot(3,1,1), plot(t,u,'k'), hold on


noise=1;
ut=fft(u);
utn=ut+noise*(randn(1,n)+i*randn(1,n));
un=ifft(utn);
figure(1), subplot(3,1,2), plot(t,un,'k'), hold on

noise=5;
ut=fft(u);
utn=ut+noise*(randn(1,n)+i*randn(1,n));
un=ifft(utn);
figure(1), subplot(3,1,3), plot(t,un,'k'), hold on


%%  SIGNAL AND SPECTRUM

figure(2)
noise=10;
ut=fft(u); 
unt=ut+noise*(randn(1,n)+i*randn(1,n)); 
un=ifft(unt);
subplot(2,1,1), plot(t,abs(un),'k')
axis([-30 30 0 2])
xlabel('time (t)'), ylabel('|u|')
subplot(2,1,2) 
plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k') 
axis([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

%% FILTER ON BAND

figure(3)
subplot(3,1,1)
filter=exp(-0.2*(k).^2); 
unft=filter.*unt; 
unf=ifft(unft);
plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k')
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',[2])
axis([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'k')
axis([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

subplot(3,1,3)
plot(t,unf,'k','Linewidth',[2]), hold on
plot(t,u,'k')
plot([-30 30],[0.5 0.5],'k:')
axis([-30 30 0 1.2])
xlabel('time (t)'), ylabel('|u|')

%% FILTER OFF BAND

figure(4)
subplot(3,1,1)
filter=exp(-0.2*(k-15).^2); 
unft=filter.*unt; 
unf=ifft(unft);
plot(fftshift(k),abs(fftshift(unt))/max(abs(fftshift(unt))),'k')
hold on
plot(fftshift(k),fftshift(filter),'k','Linewidth',[2])
axis([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

subplot(3,1,2)
plot(fftshift(k),abs(fftshift(unft))/max(abs(fftshift(unft))),'k')
axis([-25 25 0 1])
xlabel('wavenumber (k)'), ylabel('|ut|/max(|ut|)')

subplot(3,1,3)
plot(t,unf,'k','Linewidth',[2]), hold on
plot(t,u,'k')
plot([-30 30],[0.5 0.5],'k:')
axis([-30 30 0 1.2])
xlabel('time (t)'), ylabel('|u|')
