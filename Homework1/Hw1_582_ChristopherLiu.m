clear all; close all; clc;
load Testdata

L=15; % spatial domain
n=64; % Fourier modes

x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k); %wave numbers

%mesh grids used for plotting
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%%
%Averaging the spectrum and finding the center frequency.
    
%reshape each row into a 64x64x64 and take the sum of each row.
sum = zeros(64,64,64);
for j = 1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    sum = sum + fftn(Un);
end

avg = sum/20; %divide total sum by 20 to get average
[cfreq,idx] = max(abs(avg(:))); %finding the center freq
[x_c,y_c,z_c] = ind2sub([64 64 64],idx); %finding the index for max freq

%%
% creating the filter and plotting position

%Where the filter will be centered at in the Fourier domain
kx_i = Kx(x_c,y_c,z_c);
ky_i = Ky(x_c,y_c,z_c);
kz_i = Kz(x_c,y_c,z_c);

%Create a gaussian filter in the fourier domain centered at 
%(kx_i,ky_i,kz_i)
width = 1;
filter=exp(-width*((Kx - kx_i).^2+(Ky-ky_i).^2+(Kz-kz_i).^2));

pos = zeros(20,3); %position of the marble 
for j = 1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    ut = fftn(Un);
    ut_filter = ut.*filter; %Applying the filter to each time 
    un_filter = ifftn(ut_filter);
    
    %finding where the marble is and converting it back to x,y,z
    [max_pos, m_indx] = max(abs(un_filter(:))); 
    [x_m,y_m,z_m] = ind2sub([64 64 64], m_indx);
    pos(j,:) = [X(x_m, y_m, z_m),Y(x_m, y_m, z_m),Z(x_m, y_m, z_m)];
    
end

%final position of the marble
final_pos = pos(20,:);

%rounding the postions for report table
round_pos = round(pos,3,'significant');

%plotting the trajector of the marble
figure(1)
plot3(pos(:,1),pos(:,2),pos(:,3))
hold on
plot3(final_pos(1),final_pos(2),final_pos(3),'x')
title('Marble Path')
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')

%%
% Supplmentary plots for report

% Plot for the filter
filt_shift = fftshift(filter);
figure(2)
isosurface(Kx,Ky,Kz,filt_shift,0.4)
title('Filter')
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')
grid on

%Plot the avg in fourier space (normalized)
figure(3)
isosurface(Kx,Ky,Kz,fftshift(abs(avg))./max(abs(avg(:))),0.5)
grid on
title('Data after Averaging in Fourier Space')
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')

%%
% Plotting the initial noisy data
   
Un(:,:,:)=reshape(Undata(20,:),n,n,n);
figure(4)
isosurface(X,Y,Z,abs(Un),0.4)
axis([-10 15 -5 10 -10 15]), grid on, drawnow, title('Unfiltered Data')
xlabel('X')
ylabel('Y')
zlabel('Z')



%%
% Plotting the transformed data

Un(:,:,:)=reshape(Undata(20,:),n,n,n);
figure(5)
Ut = abs(fftn(Un));
isosurface(Kx,Ky,Kz,fftshift(Ut)./max(abs(Ut(:))),0.4)
grid on, drawnow, title('Transformed Data')
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')


%%
% Plotting the filtered data (untransformed)

for j = 1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    ut2 = fftn(Un);
    ut_filter2 = ut2.*filter; %Applying the filter to each time 
    un_filter2 = ifftn(ut_filter2);

    figure(6)
    isosurface(X,Y,Z,abs(un_filter2),0.1)
    grid on, drawnow, title('Filtered Data')
    axis([-10 15 -5 10 -10 15])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end


