clear all; clc; close all;

% Importing the cropped face 
D = 'C:\Users\Chris\Dropbox\School\WIN2020\AMATH 582\Homework\Homework4\yalefaces_cropped\CroppedYale';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
index = 1;
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.pgm')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        A = imread(F,'pgm');
        X_f(:,index) = reshape(A,[32256,1]);
        index = index+1;
    end
end
%%
% Subtract average from the face for cropped
X_f_avg = zeros(192*168, 2414);

for i = 1:192*168
    X_f_avg(i,:) = X_f(i,:) - mean(X_f(i,:));
end
%%
% SVD for cropped images
[Uc,Sc,Vc] = svd(double(X_f),0);
[Uc_a,Sc_a,Vc_a] = svd(X_f_avg,0);

%%
% Plot singular values

dSc = diag(Sc);
dSc_a = diag(Sc_a);

figure(1)
subplot(2,2,1)
plot(1:length(Sc),dSc,'.r');
xlabel('Modes')
ylabel('Value')
subplot(2,2,2)
semilogy(1:length(Sc),dSc,'.r');
xlabel('Modes')
ylabel('Value')

%calculate energy
energy = zeros(length(dSc),1);
e_total = sum(dSc);
for i = 1:length(dSc)
   energy(i)= sum(dSc(1:i))/e_total;
end

subplot(2,2,[3,4])
plot(1:length(dSc),energy)
xlabel('Modes')
ylabel('Cumulative Energy (%)')

%mean subtracted
figure(2)
subplot(2,2,1)
plot(1:length(Sc_a),dSc_a,'.r');
xlabel('Modes')
ylabel('Value')
subplot(2,2,2)
semilogy(1:length(Sc_a),dSc_a,'.r');
xlabel('Modes')
ylabel('Value')

%calculate energy
energy_a = zeros(length(dSc_a),1);
e_a_total = sum(dSc_a);
for i = 1:length(dSc_a)
   energy_a(i)= sum(dSc_a(1:i))/e_a_total;
end

subplot(2,2,[3,4])
plot(1:length(dSc_a),energy_a)
xlabel('Modes')
ylabel('Cumulative Energy (%)')
%%
% Rank approximation to the faces.
rank = [1 5 10 50 100];
figure(3)
subplot(3,2,1)
imshow(reshape(uint8(X_f_avg(:,1)),[192,168]))
title('original')

for k = 1:length(rank)
    Xc_approx = Uc_a(:,1:rank(k))*Sc_a(1:rank(k),1:rank(k))*Vc_a(:,1:rank(k))';
    C = reshape(uint8(Xc_approx(:,1)),[192,168]);
    subplot(3,2,k+1)
    imshow(C);
    title(['Rank ', num2str(rank(k)), ' Approximation'])
end

%%
% Eigenfaces for cropped
j = 4;
for k = 1:j
    maxv = max(Uc_a(:,k));
    minv = min(Uc_a(:,k));
    
    figure(10)
    subplot(1,j,k)
    imshow(reshape(Uc_a(:,k),[192,168]),[minv maxv])
    title(['Eigenface ', num2str(k)])
end

%%
% Importing the uncropped faces
D = 'C:\Users\Chris\Dropbox\School\WIN2020\AMATH 582\Homework\Homework4\yalefaces';
S = dir(fullfile(D,'*'));

C = {S(~[S.isdir]).name}; % files in subfolder.
for jj = 1:numel(C)
    F = fullfile(D,C{jj});
    A = imread(F);
    dim = size(A);
    X_f_full(:,jj) = reshape(A,[dim(1)*dim(2),1]);
end

%%
% Subtract average from the face for uncropped
X_f_full_avg = zeros(77760, 165);

for i = 1:77760
    X_f_full_avg(i,:) = X_f_full(i,:) - mean(X_f_full(i,:));
end
%%
% SVD for uncropped images
[Uf,Sf,Vf] = svd(double(X_f_full),0);
[Uf_a,Sf_a,Vf_a] = svd(X_f_full_avg,0);

%%
% Plot singular values

dSf = diag(Sf);
dSf_a = diag(Sf_a);

figure(4)
subplot(2,2,1)
plot(1:length(Sf),dSf,'.r');
xlabel('Modes')
ylabel('Value')
subplot(2,2,2)
semilogy(1:length(Sf),dSf,'.r');
xlabel('Modes')
ylabel('Value')

%calculate energy
energy = zeros(length(dSf),1);
e_total = sum(dSf);
for i = 1:length(dSf)
   energy(i)= sum(dSf(1:i))/e_total;
end

subplot(2,2,[3,4])
plot(1:length(dSf),energy)
xlabel('Modes')
ylabel('Cumulative Energy (%)')

%mean subtracted
figure(5)
subplot(2,2,1)
plot(1:length(Sf_a),dSf_a,'.r');
xlabel('Modes')
ylabel('Value')
subplot(2,2,2)
semilogy(1:length(Sf_a),dSf_a,'.r');
xlabel('Modes')
ylabel('Value')

%calculate energy
energy_a = zeros(length(dSf_a),1);
e_a_total = sum(dSf_a);
for i = 1:length(dSf_a)
   energy_a(i)= sum(dSf_a(1:i))/e_a_total;
end

subplot(2,2,[3,4])
plot(1:length(dSf_a),energy_a)
xlabel('Modes')
ylabel('Cumulative Energy (%)')
%%
% Rank approximation to the faces (full).
rank = [1 5 10 50 100];
figure(6)
subplot(3,2,1)
imshow(reshape(uint8(X_f_full_avg(:,1)),[243,320]))
title('original')

for k = 1:length(rank)
    Xf_approx = Uf_a(:,1:rank(k))*Sf_a(1:rank(k),1:rank(k))*Vf_a(:,1:rank(k))';
    C = reshape(uint8(Xf_approx(:,1)),[243,320]);
    subplot(3,2,k+1)
    imshow(C);
    title(['Rank ', num2str(rank(k)), ' Approximation'])
end

%%
j = 4; 
for k = 1:j
    maxv = max(Uf_a(:,k));
    minv = min(Uf_a(:,k));
    
    figure(10)
    subplot(1,j,k)
    imshow(reshape(Uf_a(:,k),[243,320]),[minv maxv])
    title(['Eigenface ', num2str(k)])
end