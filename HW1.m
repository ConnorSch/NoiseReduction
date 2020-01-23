%Connor Schleicher AMATH 582 HW 1

%% Initialize the program 
% initial code from HW prompt for initializing data
clear all; close all; clc;
load Testdata

L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% initial figure to show the unusefulness of the noisy signal
figure(1)
Un(:,:,:)=reshape(Undata(1,:),n,n,n);
close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Raw data at first signal pass')

%% Average the noisy singal and find central frequency
% iterate through the first 20 realizations
% averaging the signal to find the central frequency
% iso surface to create visualization to see the averaging happen

Uave = zeros(n,n,n); % initialize the averaged variable
for j=1:20
    figure(2);
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Uave = Uave + fftn(Un); % averaging the shifted transform of the data
    isosurface(Kx,Ky,Kz,abs(fftshift(Uave))/max(abs(fftshift(Uave(:)))),0.6)
    axis([-6 6 -6 6 -6 6]), grid on, drawnow
    xlabel('Kx'); ylabel('Ky'); zlabel('Kz');
    title('Averaging of Frequency Data')
    pause(1)
    if j < 20 % closes figure for next iteration isosurface won't erase old data
        close(2);
    end 
end

Uave = abs(fftshift(Uave)); %shift and take abs value of averaged data

% find the maximum value of the averaged signal
[value,idx] = max(Uave(:));
[r,c,p] = ind2sub(size(Uave),idx);

% sanity check of max value calculation 
fprintf('Max value of function (abs): %d \n',abs(value));
fprintf('Uave value at position r,c,p (abs): %d \n', abs(Uave(r,c,p)));

% print locaton of central frequency
fprintf('Location of central frequency: %d, %d, %d \n', r,c,p);

%% Create Filter and Plot Marble Path
% create the filter to de-noise the data around the central frequency
tau = 0.5; % bandwidth of filter
% a simple gaussian filter across the three spatial dimensions
filter = exp(-tau*((Kx - ks(c)).^2 + (Ky - ks(r)).^2 + (Kz - ks(p)).^2));

% initialize some vectors to store values from each signal pass
value2 = zeros(1,20);
idx2 = zeros(1,20);
r2 = zeros(1,20);
c2 = zeros(1,20);
p2 = zeros(1,20);
locations = zeros(20,3);

% Apply the filter to each signal (loop of 20 iterations)
for i = 1:20
    Un(:,:,:)=reshape(Undata(i,:),n,n,n); 
    
    Unt = fftshift(fftn(Un)); % applying fft and shift to the raw data
    
    Utnf = Unt.*filter; % apply the Gaussian filter to each element of the transformed data
   
    U = ifftn(Utnf); % inverse transform to get back to the original signal
    
    % getting value and index position of the max value
    [value2(i), idx2(i)] = max(U(:)); 
    
    % translating the index to the row, column, slice's of the matrix
    [r2(i),c2(i),p2(i)] = ind2sub(size(U),idx2(i)); 
    locations(i,:) = [X(r2(i),c2(i),p2(i)), Y(r2(i),c2(i),p2(i)),...
        Z(r2(i),c2(i),p2(i))]; % sets XYZ coordinates in a matrix to be used for plotting
end

figure (3) 
% Plotting the position of the marble for each singal recording
% Marble is starting at the top (highest Z position) and working downwards
plot3(locations(:,1),locations(:,2),locations(:,3),...
    '-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
hold on

% Replotting the last point of the marble to color it red to stand out 
plot3(locations(end,1),locations(end,2),locations(end,3),...
    '-o','Color','r','MarkerSize',10,'MarkerFaceColor','#FFD9D9')
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Movement of Marble');
grid on;

fprintf('The final location of the marble is (X,Y,Z): %0.4f, %0.4f, %0.4f \n' ...
    ,locations(end,1),locations(end,2),locations(end,3))
