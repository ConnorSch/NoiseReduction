clear all; close all; clc

L=30; % time slot
n=512; % Fourier modes
t2=linspace(-L,L,n+1); t=t2(1:n); 
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k);
noise=10;

labels=['(a)';'(b)';'(c)';'(d)']; 
realize=[1 2 5 100];

for jj=1:length(realize)
  u=sech(t); ave=zeros(1,n); ut=fft(u);
  for j=1:realize(jj)
    utn(j,:)=ut+noise*(randn(1,n)+i*randn(1,n)); 
    ave=ave+utn(j,:); 
    dat(j,:)=abs(fftshift(utn(j,:)))/max(abs(utn(j,:))); 
    un(j,:)=ifft(utn(j,:));
  end
  ave=abs(fftshift(ave))/realize(jj);
  subplot(4,1,jj)
  plot(ks,ave/max(ave),'k') 
  set(gca,'Fontsize',[15])
  axis([-20 20 0 1]) 
  text(-18,0.7,labels(jj,:),'Fontsize',[15]) 
  ylabel('|fft(u)|','Fontsize',[15])
end

hold on 
plot(ks,abs(fftshift(ut))/max(abs(ut)),'k:','Linewidth',[2]) 
set(gca,'Fontsize',[15])
xlabel('frequency (k)')


%%  SIGNAL

figure(2)
waterfall(ks,1:8,dat(1:8,:)), 
colormap([0 0 0]), 
view(-15,80) 
set(gca,'Fontsize',[15],'Xlim',[-28 28],'Ylim',[1 8]) 
xlabel('frequency (k)')
ylabel('realization')
zlabel('|fft(u)|')

%% MOVING SIGNAL

figure(3)
slice=[0:0.5:10];
[T,S]=meshgrid(t,slice); 
[K,S]=meshgrid(k,slice);
U=sech(T-10*sin(S)).*exp(i*0*T);

subplot(2,1,1)
waterfall(T,S,U), colormap([0 0 0]), view(-15,70) 
set(gca,'Fontsize',[15],'Xlim',[-30 30],'Zlim',[0 2]) 
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|')
for j=1:length(slice)
    Ut(j,:)=fft(U(j,:));
    Kp(j,:)=fftshift(K(j,:)); 
    Utp(j,:)=fftshift(Ut(j,:)); 
    Utn(j,:)=Ut(j,:)+noise*(randn(1,n)+i*randn(1,n)); 
    Utnp(j,:)=fftshift(Utn(j,:))/max(abs(Utn(j,:))); 
    Un(j,:)=ifft(Utn(j,:));
end
subplot(2,1,2)
waterfall(Kp,S,abs(Utp)/max(abs(Utp(1,:))))
colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',[15],'Xlim',[-28 28])
xlabel('frequency (k)'), ylabel('realizations'), zlabel('|fft(u)|')

%%

figure(4)
subplot(2,1,1)
waterfall(T,S,abs(Un)), 
colormap([0 0 0]), view(-15,70) 
set(gca,'Fontsize',[15],'Xlim',[-30 30],'Zlim',[0 2]) 
xlabel('time (t)'), ylabel('realizations'), zlabel('|u|') 
subplot(2,1,2)
waterfall(Kp,S,abs(Utnp)), colormap([0 0 0]), view(-15,70)
set(gca,'Fontsize',[15],'Xlim',[-28 28])
xlabel('frequency (k)'), ylabel('realizations'), zlabel('|fft(u)|')

%% AVERAGE SIGNAL

figure(5)
Uave=zeros(1,n); 
Utave=zeros(1,n); 
for j=1:length(slice)
    Uave=Uave+Un(j,:);
    Utave=Utave+Utn(j,:); 
end
Uave=Uave/length(slice); 
Utave=fftshift(Utave)/length(slice);
subplot(2,1,1) 
plot(t,abs(Uave),'k') 
set(gca,'Fontsize',[15]) 
xlabel('time (t)'), ylabel('|u|')

subplot(2,1,2)
plot(ks,abs(Utave)/max(abs(Utave)),'k'), hold on 
plot(ks,abs(fftshift(Ut(1,:)))/max(abs(Ut(1,:))),'k:','Linewidth',[2]) 
axis([-20 20 0 1])
set(gca,'Fontsize',[15])
xlabel('frequency (k)'), ylabel('|fft(u)|')

