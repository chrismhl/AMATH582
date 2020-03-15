clear all; clc; close all;
 
%%
% Importing music, creating spectrograms and ouputting as .mat
% This was done to make it easier to move data between laptop and desktop
D = 'C:\XXX'; %change as needed
S = dir(fullfile(D,'*.wav'));
X = [];

C = {S(~[S.isdir]).name}; % files in folder.
for jj = 1:10
    start = randi([15,25]); %change as needed
    F = fullfile(D,C{jj});
    [y, Fs] = audioread(F);
    
    spec=[];
    incr = 0.2;
    tslide_p=0:incr:5;
    width = 50; %change as needed
    yc = y(start*Fs:100:(start+5)*Fs,1);
    if mod(length(yc),10) == 1
        yc = yc(1:length(yc)-1);
    end
    n = length(yc);

    t2=linspace(0,5,n+1);
    t2 = t2(1:n);
    kp=(1/5)*[0:n/2-1 -n/2:-1]; 
    ksp=fftshift(kp);


    for j=1:length(tslide_p)
        g=exp(-width*(t2-tslide_p(j)).^2); % Gabor 
        Pg=g.*yc.'; 
        Pgt=fft(Pg); 
        spec=[spec; abs(fftshift(Pgt))]; 
    end

    dim = size(spec);
    x = reshape(spec,[dim(1)*dim(2), 1]);
    
    X = [X, x];
end
dp = X;
save('dp.mat','dp');

%%
% Load data and perform SVD for Test 1
clear all; clc; close all

load acdc.mat
load acdc2.mat
load dp.mat
load dp2.mat
load jb.mat
load jb2.mat

A = [acdc, acdc2, dp, dp2, jb, jb2];

[U, S, V] = svd(A,0);

test = [V(1:5,1:10); V(11:15,1:10); V(21:25,1:10); V(31:35,1:10);V(41:45,1:10); V(51:55,1:10)];
train = [V(6:10,1:10); V(16:20,1:10); V(26:30,1:10); V(36:40,1:10);V(46:50,1:10); V(56:60,1:10)];

group = ones(30,1);
group(1:10) = 1.*group(1:10);
group(11:20) = 2.*group(11:20);
group(21:30) = 3.*group(21:30);

[class, err, logp, coeff] = classify(test, train, group);

correct =0;
for i = 1:length(class)
    if class(i) == group(i)
        correct = correct +1;
    end
end
error_1 = correct/length(group);

figure(1)
plot(group,'b.','MarkerSize',20);
hold on
plot(class,'r.','MarkerSize',10);
hold off
yticks([1,2,3])
legend('True Solution','Classified Data','Location','NorthWest')
title('Test 1')
xlabel('Songs')
ylabel('Labels')
%%
% Load data and perform SVD for Test 2
close all; clear all; clc;

load gd1.mat
load gd2.mat
load yc1.mat
load yc2.mat
load nufan1.mat
load nufan2.mat

A = [gd1,gd2,yc1,yc2,nufan1,nufan2];

[U, S, V] = svd(A,0);

test = [V(1:5,1:10); V(11:15,1:10); V(21:25,1:10); V(31:35,1:10);V(41:45,1:10); V(51:55,1:10)];
train = [V(6:10,1:10); V(16:20,1:10); V(26:30,1:10); V(36:40,1:10);V(46:50,1:10); V(56:60,1:10)];

group = ones(30,1);
group(1:10) = 1.*group(1:10);
group(11:20) = 2.*group(11:20);
group(21:30) = 3.*group(21:30);

[class, err, logp, coeff] = classify(test, train, group);

correct =0;
for i = 1:length(class)
    if class(i) == group(i)
        correct = correct +1;
    end
end
error_2 = correct/length(group);

figure(2)
plot(group,'b.','MarkerSize',20);
hold on
plot(class,'r.','MarkerSize',10);
hold off
yticks([1,2,3])
legend('True Solution','Classified Data','Location','NorthWest')
title('Test 2')
xlabel('Songs')
ylabel('Labels')
%%
% Load data and create data and perform SVD for Test 3
close all; clear all; clc;

load gd1.mat
load gd2.mat
load yc1.mat
load yc2.mat
load nufan1.mat
load nufan2.mat
load bob.mat
load kpop.mat
load b182.mat
load djdeck.mat
load ok.mat

punk = [gd1,yc1,nufan1,b182];
rap = [bob, djdeck, ok];

pperm = randperm(size(punk,2),20);
rperm = randperm(size(rap,2),20);
kperm = randperm(size(kpop,2),20);

A = [punk(:,pperm),rap(:,rperm),kpop(:,kperm)];
[U, S, V3] = svd(A,0);

save('V3.mat','V3');
save('S3.mat','S');
%%
% Classify test 3
close all; clear all; clc;
load V3.mat
V = V3;

train = [V(1:5,1:10); V(11:15,1:10); V(21:25,1:10); V(31:35,1:10);V(41:45,1:10); V(51:55,1:10)];
test = [V(6:10,1:10); V(16:20,1:10); V(26:30,1:10); V(36:40,1:10);V(46:50,1:10); V(56:60,1:10)];

group = ones(30,1);
group(1:10) = 1.*group(1:10);
group(11:20) = 2.*group(11:20);
group(21:30) = 3.*group(21:30);

[class, err, logp, coeff] = classify(test, train, group);

correct =0;
for i = 1:length(class)
    if class(i) == group(i)
        correct = correct +1;
    end
end
error_3 = correct/length(group);

figure(3)
plot(group,'b.','MarkerSize',20);
hold on
plot(class,'r.','MarkerSize',10);
hold off
yticks([1,2,3])
legend('True Solution','Classified Data','Location','NorthWest')
title('Test 3')
xlabel('Songs')
ylabel('Labels')
%%
%Plot singular values and energy 
dS = diag(S);
figure(1)
plot(1:length(dS), dS)

%calculate energy
energy = zeros(length(dS),1);
e_total = sum(dS);
for i = 1:length(dS)
   energy(i)= sum(dS(1:i))/e_total;
end
figure(2)
plot(1:length(dS),energy)


%%
% Test code for the spectrogram. delete later if needed
pcolor(tslide_p,ksp,(spec.'/max(max(abs(spec))))) 
shading interp
ylim([0 max(ksp)])
set(gca,'Fontsize',[14]) 
colormap(hot)

    
