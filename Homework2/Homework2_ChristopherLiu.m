clear all; clc; close all;

load handel.mat

v = y'/2;

L=round(max(abs((1:length(v))/Fs))); n=length(v)-1;
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(1/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

v = v(1:n);
%%
% Initial plots of the signal and the FFT
figure(1)
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Handel Signal, v(n)');

figure(2)
Vt = fft(v(1:n));
plot(ks,abs(fftshift(Vt))/max(abs(Vt)),'k');
set(gca,'Fontsize',[14])
title('Handel Signal after FFT');
xlabel('frequency (\omega)'), ylabel('FFT(S)')

%%
% Sliding Gabor Transform with fixed width

%figure(3)
Sgt_spec=[]; 
tslide=0:0.1:L;
width = [0.5,1,5,10]; %change as needed

for k = 1:length(width)
    for j=1:length(tslide)
        g=exp(-width(k)*(t-tslide(j)).^2); % Gabor 
        Vg=g.*v; 
        Vgt=fft(Vg); 
        Sgt_spec=[Sgt_spec; 
        abs(fftshift(Vgt))]; 
    end
    
    subplot(2,2,k)
    pcolor(tslide,ks,Sgt_spec.'), 
    shading interp 
    set(gca,'Fontsize',[14])
    ylim([0 4000])
    xlabel('Time(s)')
    ylabel('Frequency (\omega)')
    title(['Width = ' num2str(width(k))])
    colormap(hot)
    Sgt_spec=[];

end

%%
%Over and Undersampling

%figure(4)
Sgt_spec2=[];
increment = [0.05 0.1 0.5 1]; %change increment

width2 = 1; %change as needed
for k = 1:length(increment)
    tslide2=0:increment(k):L;
    
    for j=1:length(tslide2)
        g=exp(-width2*(t-tslide2(j)).^2); % Gabor 
        Vg=g.*v; 
        Vgt=fft(Vg); 
        Sgt_spec2=[Sgt_spec2; 
        abs(fftshift(Vgt))]; 
    end
    
    subplot(2,2,k)
    pcolor(tslide2,ks,Sgt_spec2.'), 
    shading interp 
    set(gca,'Fontsize',[14])
    ylim([0 4000])
    xlabel('Time(s)')
    ylabel('Frequency (\omega)')
    title(['Samples = ' num2str(length(tslide2))])
    colormap(hot)
    Sgt_spec2=[];
end

%%
% Mexican hat wavelet

Sgt_spec3=[]; 
tslide3=0:0.1:L;
widthmh = 1;  %change as needed

for j=1:length(tslide3)
    mhat = (1-(t-tslide3(j)).^2).*exp(-widthmh*(t-tslide3(j)).^2);
    Vm=mhat.*v; 
    Vmt=fft(Vm); 
    Sgt_spec3=[Sgt_spec3; 
    abs(fftshift(Vmt))]; 
end

%%
% Plotting spectrogram for mexican hat wavelet

figure(7)
pcolor(tslide3,ks,Sgt_spec3.'), 
shading interp 
set(gca,'Fontsize',[14]) 
ylim([0 4000])
xlabel('Time(s)')
ylabel('Frequency (\omega)')
title('Mexican Hat filter')
colormap(hot)

%%
% Sliding step-function

widthshn = 500;
m = n; %lenght of the filter
step = 500;
Sgt_spec4=[];
Vs=zeros(n,1);

for j = 1:step:m
    if j < widthshn+1
        Vs = zeros(m,1);
        Vs(j:1:j+widthshn) = v(j:1:j+widthshn);
        
        if j>1
        Vs(1:j) = v(1:j);
        end
        
    elseif (j+widthshn) >= m
        Vs = zeros(m,1);
        Vs(j - widthshn:1:m) = v(j - widthshn:1:m);
    else   
        Vs = zeros(m,1);
        Vs(j - widthshn:1:j+widthshn) = v(j - widthshn:1:j+widthshn);
    end
    
    Vst=fft(Vs); 
    Sgt_spec4=[Sgt_spec4;abs(fftshift(Vst.'))];
end

%%
% Plotting spectrogram for shannon window

figure(8)
tslide4=linspace(0,9,length(1:step:m));
pcolor(tslide4,ks,Sgt_spec4.'), 
shading interp 
set(gca,'Fontsize',[14]) 
ylim([0 4000])
xlabel('Time(s)')
ylabel('Frequency (\omega)')
title('Step-function (Shannon) Window')
colormap(hot)

%%
% Part 2
clc; clear all; close all;

figure(9)
tr_piano=16; % record time in seconds
y2=audioread('music1.wav'); Fs2=length(y2)/tr_piano;
plot((1:length(y2))/Fs2,y2);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)'); drawnow

figure(10)
tr_rec=14; % record time in seconds
y3=audioread('music2.wav'); Fs3=length(y3)/tr_rec;
plot((1:length(y3))/Fs3,y3);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');

%%
% Frequencies for the piano overtones.

n=length(y2);
t2=linspace(0,tr_piano,n+1); tp=t2(1:n); 
kp=(1/tr_piano)*[0:n/2-1 -n/2:-1]; 
ksp=fftshift(kp);

subplot(1,2,1)
y2t = fft(y2(1:n));
plot(ksp,abs(fftshift(y2t))/max(abs(y2t)),'k');
xlim([0 5000])
title('Frequency content for piano')
%xlim([220 350])
set(gca,'Fontsize',[14])
xlabel('frequency (\omega)'), ylabel('FFT(S)')

% Frequencies for the recorder overtones.

n=length(y3);
t2=linspace(0,tr_rec,n+1); tr=t2(1:n); 
kr=(1/tr_rec)*[0:n/2-1 -n/2:-1]; 
ksr=fftshift(kr);

subplot(1,2,2)
y3t = fft(y3(1:n));
plot(ksr,abs(fftshift(y3t))/max(abs(y3t)),'k');
xlim([0 5e3])
title('Frequency content for recorder')
set(gca,'Fontsize',[14])
xlabel('Frequency (\omega)'), ylabel('FFT(S)')
%%
%Gabor filter for piano

piano_spec=[];
incr = 0.2;
tslide_p=0:incr:tr_piano;
width = 25; %change as needed
    
for j=1:length(tslide_p)
    g=exp(-width*(tp-tslide_p(j)).^2); % Gabor 
    Pg=g.*y2.'; 
    Pgt=fft(Pg); 
    piano_spec=[piano_spec; 
    abs(fftshift(Pgt))]; 
end

%%
%Spectrogram for piano
figure(11)

pcolor(tslide_p,ksp,(piano_spec.'/max(max(abs(piano_spec))))) 
shading interp 
set(gca,'Fontsize',[14]) 
ylim([0 1000])
title('Score for Piano')
xlabel('Time(s)')
ylabel('Frequency (\omega)')
colormap(hot)
%%
%Gabor filter for recorder

rec_spec=[];
incr = 0.1;
tslide_r=0:incr:tr_rec;
width = 20; %change as needed

for j=1:length(tslide_r)
    g=exp(-width*(tr-tslide_r(j)).^2); % Gabor 
    Rg=g.*y3.'; 
    Rgt=fft(Rg); 
    rec_spec=[rec_spec; 
    abs(fftshift(Rgt))]; 
end

%%
%Spectrogram for recorder
figure(12)

pcolor(tslide_r,ksr,(rec_spec.'/max(max(abs(rec_spec))))) 
shading interp 
set(gca,'Fontsize',[14]) 
ylim([0 2e3])
xlabel('Time(s)')
ylabel('Frequency (\omega)')
title('Score for Recorder')
colormap(hot)

