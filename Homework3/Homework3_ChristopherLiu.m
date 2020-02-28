clear all; clc ;close all;

load('cam1_1.mat')
load('cam1_2.mat')
load('cam1_3.mat')
load('cam1_4.mat')
load('cam2_1.mat')
load('cam2_2.mat')
load('cam2_3.mat')
load('cam2_4.mat')  
load('cam3_1.mat')
load('cam3_2.mat')
load('cam3_3.mat')
load('cam3_4.mat')

%%
% Tracking for camera 1 part 1
[x1_1, y1_1] = track_paintcan_gray(vidFrames1_1,[321,228],[20,20,20,20]);
    
%%
% Playing Camera 1 Part 1
for i=1:226
    figure(1)
    imshow(vidFrames1_1(:,:,:,i))
    hold on
    plot(x1_1(i),y1_1(i),'r*')
    hold off
end

%%
% Tracking for camera 2 part 1
[x2_1, y2_1] = track_paintcan_gray(vidFrames2_1,[274,274],[20,20,20,20]);
    
%%
% Play Camera 2 Part 1
for i=1:284
    figure(2)
    imshow(uint8(vidFrames2_1(:,:,:,i)))
    hold on
    plot(x2_1(i),y2_1(i),'r*')
    hold off
end

%%
% Tracking for camera 3 part 1
[x3_1, y3_1] = track_paintcan_gray(vidFrames3_1,[318,271],[15,15,15,15]);

%%
% Play Camera 3 Part 1
for i=1:232
    figure(2)
    imshow(uint8(vidFrames3_1(:,:,:,i)))
    hold on
    plot(x3_1(i),y3_1(i),'r*')
    hold off
end

%%
%Plot position of all 3 cameras, test 1
close all

figure(3)

subplot(3,1,1)
plot(1:226,y1_1,'r')
hold on
plot(1:284,y2_1,'g')
plot(1:232,x3_1,'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')

%truncated
subplot(3,1,2)
plot(1:226,y1_1,'r')
hold on
plot(1:226,y2_1(10:235),'g')
plot(1:226,x3_1(1:226),'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')
title('Truncated')
%%
% SVD part 1
A1(1,:) = x1_1 - mean(x1_1);
A1(2,:) = y1_1 - mean(y1_1);
A1(3,:) = x2_1(10:235) - mean(x2_1(10:235));
A1(4,:) = y2_1(10:235) - mean(y2_1(10:235));
A1(5,:) = x3_1(1:226) - mean(x3_1(1:226));
A1(6,:) = y3_1(1:226) - mean(y3_1(1:226));
[U1,s1,V1] = svd(A1.',0);

sig1 = diag(s1);
energy1_1 = sig1(1)/sum(sig1);
energy2_1 = sum(sig1(1:2))/sum(sig1);
energy3_1 = sum(sig1(1:3))/sum(sig1);

figure(5)
subplot(3,2,1)
plot(1:length(sig1),sig1,'r.','MarkerSize',20)
title('Singular Values')
subplot(3,2,2)
plot(1:length(sig1),sig1,'r.','MarkerSize',20)
set(gca, 'YScale', 'log')
title('Singular Values (semi-log)')

subplot(3,1,2)
plot(1:6,V1(:,1),'.-',1:6,V1(:,2),'.-',1:6,V1(:,3),'.-')
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
xticks([1 2 3 4 5 6])
title('Linear POD Modes')

subplot(3,1,3)
plot(1:226,U1(:,1),1:226,U1(:,2),1:226,U1(:,3))
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
ylabel('Displacement')
xlabel('Frame')
title('Time Evolution Behavior')
%%
% Tracking for camera 1 part 2
[x1_2, y1_2] = track_paintcan(vidFrames1_2,[325,308],[17,17*0.9,15,15*0.9]);

%%
% Tracking Camera 1 Part 2
for i=1:314
    figure(2)
    imshow(uint8(vidFrames1_2(:,:,:,i)))
    hold on
    plot(x1_2(i),y1_2(i),'r*')
    hold off
end

%%
% Tracking for camera 2 part 2
[x2_2, y2_2] = track_paintcan(vidFrames2_2,[314,357],[34,34,34,35]);

%%
% Tracking Camera 2 Part 2
for i=1:356
    figure(2)
    imshow(uint8(vidFrames2_2(:,:,:,i)))
    hold on
    plot(x2_2(i),y2_2(i),'r*')
    hold off
end

%%
% Tracking for camera 3 part 2
[x3_2, y3_2] = track_paintcan_gray(vidFrames3_2,[349,245],[44,45,44,44]);

%%
% Tracking Camera 3 Part 2
for i=1:327
    figure(2)
    imshow(uint8(vidFrames3_2(:,:,:,i)))
    hold on
    plot(x3_2(i),y3_2(i),'r*')
    hold off
end

%%
%Plot position of all 3 cameras part 2
close all

figure(5)

subplot(2,1,1)
plot(1:314,y1_2,'r')
hold on
plot(1:356,y2_2,'g')
plot(1:327,x3_2,'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')

%truncated
subplot(2,1,2)
plot(1:314,y1_2,'r')
hold on
plot(1:314  ,y2_2(15:328),'g')
plot(1:314,x3_2(1:314),'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')
title('Truncated')

%%
% SVD part 2
A2 = zeros(6,314);
A2(1,:) = x1_2 - mean(x1_2);
A2(2,:) = y1_2 - mean(y1_2);
A2(3,:) = x2_2(15:328) - mean(x2_2(15:328));
A2(4,:) = y2_2(15:328) - mean(y2_2(15:328));
A2(5,:) = x3_2(1:314) - mean(x3_2(1:314));
A2(6,:) = y3_2(1:314) - mean(y3_2(1:314));
[U2,s2,V2] = svd(A2.',0);

sig2 = diag(s2);
energy1_2 = sig2(1)/sum(sig2);
energy2_2 = sum(sig2(1:2))/sum(sig2);
energy3_2 = sum(sig2(1:3))/sum(sig2);

figure(6)
subplot(3,2,1)
plot(1:length(sig2),sig2,'r.','MarkerSize',20)
title('Singular Values')
subplot(3,2,2)
plot(1:length(sig2),sig2,'r.','MarkerSize',20)
set(gca, 'YScale', 'log')
title('Singular Values (semi-log)')

subplot(3,1,2)
plot(1:6,V2(:,1),'.-',1:6,V2(:,2),'.-',1:6,V2(:,3),'.-')
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
xticks([1 2 3 4 5 6])
title('Linear POD Modes')

subplot(3,1,3)
plot(1:314,U2(:,1),1:314,U2(:,2),1:314,U2(:,3))
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
ylabel('Displacement')
xlabel('Frame')
title('Time Evolution Behavior')

%%
% Tracking for camera 1 part 3, tracking pink

[x1_3, y1_3] = track_paintcan_pink(vidFrames1_3,[330,290],[20,20,20,20]);
%%
% Tracking Camera 1 Part 3
for i=1:239
    figure(2)
    imshow(uint8(vidFrames1_3(:,:,:,i)))
    hold on
    plot(x1_3(i),y1_3(i),'r*')
    hold off
end

%%
% Tracking for camera 2 part 3, tracking pink

[x2_3, y2_3] = track_paintcan_pink(vidFrames2_3,[252,294],[20,20,20,20]);
%%
% Tracking Camera 2 Part 3
for i=1:281
    figure(2)
    imshow(uint8(vidFrames2_3(:,:,:,i)))
    hold on
    plot(x2_3(i),y2_3(i),'r*')
    hold off
end

%%
% Tracking for camera 3 part 3
[x3_3, y3_3] = track_paintcan_gray(vidFrames3_3,[353,229],[17,17,17,17]);

%%
% Tracking Camera 3 Part 3
for i=1:237
    figure(2)
    imshow(uint8(vidFrames3_3(:,:,:,i)))
    hold on
    plot(x3_3(i),y3_3(i),'r*')
    hold off
end

%%
%Plot position of all 3 cameras part 3
close all

figure(7)

subplot(2,1,1)
plot(1:239,y1_3,'r')
hold on
plot(1:281,y2_3,'g')
plot(1:237,x3_3,'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')

%truncated
subplot(2,1,2)
plot(1:206,y1_3(20:225),'r')
hold on
plot(1:206  ,y2_3(5:210),'g')
plot(1:206,x3_3(10:215),'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')
title('Truncated')

%%
% SVD part 3
A3 = zeros(6,206);
A3(1,:) = x1_3(20:225) - mean(x1_3(20:225));
A3(2,:) = y1_3(20:225) - mean(y1_3(20:225));
A3(3,:) = x2_3(5:210) - mean(x2_3(5:210));
A3(4,:) = y2_3(5:210) - mean(y2_3(5:210));
A3(5,:) = x3_3(10:215) - mean(x3_3(10:215));
A3(6,:) = y3_3(10:215) - mean(y3_3(10:215));
[U3,s3,V3] = svd(A3.',0);
sig3 = diag(s3);

energy1_3 = sig3(1)/sum(sig3);
energy2_3 = sum(sig3(1:2))/sum(sig3);
energy3_3 = sum(sig3(1:3))/sum(sig3);

figure(8)
subplot(3,2,1)
plot(1:length(sig3),sig3,'r.','MarkerSize',20)
title('Singular Values')
subplot(3,2,2)
plot(1:length(sig3),sig3,'r.','MarkerSize',20)
set(gca, 'YScale', 'log')
title('Singular Values (semi-log)')

subplot(3,1,2)
plot(1:6,V3(:,1),'.-',1:6,V3(:,2),'.-',1:6,V3(:,3),'.-')
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
xticks([1 2 3 4 5 6])
title('Linear POD Modes')

subplot(3,1,3)
plot(1:206,U3(:,1),1:206,U3(:,2),1:206,U3(:,3))
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
ylabel('Displacement')
xlabel('Frame')
title('Time Evolution Behavior')


%%
% Tracking for camera 1 part 4
[x1_4, y1_4] = track_paintcan_pink(vidFrames1_4,[402,263],[44,45,44,44]);

%%
% Tracking Camera 1 Part 4
for i=1:392
    figure(2)
    imshow(uint8(vidFrames1_4(:,:,:,i)))
    hold on
    plot(x1_4(i),y1_4(i),'r*')
    hold off
end

%%
% Tracking for camera 2 part 4
[x2_4, y2_4] = track_paintcan_pink(vidFrames2_4,[244,243],[44,45,44,44]);

%%
% Tracking Camera 2 Part 4
for i=1:405
    figure(2)
    imshow(uint8(vidFrames2_4(:,:,:,i)))
    hold on
    plot(x2_4(i),y2_4(i),'r*')
    hold off
end

%%
% Tracking for camera 3 part 4
[x3_4, y3_4] = track_paintcan(vidFrames3_4,[363,244],[40,40,30,30]);

%%
% Tracking Camera 3 Part 4
for i=1:394
    figure(2)
    imshow(uint8(vidFrames3_4(:,:,:,i)))
    hold on
    plot(x3_4(i),y3_4(i),'r*')
    hold off
end

%%
%Plot position of all 3 cameras part 4
close all

figure(9)

subplot(2,1,1)
plot(1:392,y1_4,'r')
hold on
plot(1:405,y2_4,'g')
plot(1:394,x3_4,'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')

%truncated
subplot(2,1,2)
plot(1:392,y1_4,'r')
hold on
plot(1:392  ,y2_4(1:392),'g')
plot(1:392,x3_4(1:392),'b')
hold off
legend('Camera 1', 'Camera 2', 'Camera 3')
title('Truncated')

%%
% SVD part 4
A4 = zeros(6,392);
A4(1,:) = x1_4 - mean(x1_4);
A4(2,:) = y1_4 - mean(y1_4);
A4(3,:) = x2_4(1:392) - mean(x2_4(1:392));
A4(4,:) = y2_4(1:392) - mean(y2_4(1:392));
A4(5,:) = x3_4(1:392) - mean(x3_4(1:392));
A4(6,:) = y3_4(1:392) - mean(y3_4(1:392));
[U4,s4,V4] = svd(A4.',0);

sig4 = diag(s4);
energy1_4 = sig4(1)/sum(sig4);
energy2_4 = sum(sig4(1:2))/sum(sig4);
energy3_4 = sum(sig4(1:3))/sum(sig4);

figure(10)
subplot(3,2,1)
plot(1:length(sig4),sig4,'r.','MarkerSize',20)
title('Singular Values')
subplot(3,2,2)
plot(1:length(sig4),sig4,'r.','MarkerSize',20)
title('Singular Values (semi-log)')
set(gca, 'YScale', 'log')

subplot(3,1,2)
plot(1:6,V4(:,1),'.-',1:6,V4(:,2),'.-',1:6,V4(:,3),'.-')
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
title('Linear POD Modes')

subplot(3,1,3)
plot(1:392,U4(:,1),1:392,U4(:,2),1:392,U4(:,3))
legend('Mode 1','Mode 2','Mode 3','Location','southeast')
ylabel('Displacement')
xlabel('Frame')
title('Time Evolution Behavior')

%%
% Plot all position vectors

figure(11)
plot(1:226,y1_1,'r')
hold on
plot(1:226,y2_1(10:235),'g')
plot(1:226,x3_1(1:226),'b')
hold off
title('Test 1: Ideal Case')
lgd = legend('Camera 1', 'Camera 2', 'Camera 3','Location','southeast');
lgd.FontSize =7;
ylabel('Displacement')
xlabel('Frame')

figure(12)
plot(1:314,y1_2,'r')
hold on
plot(1:314  ,y2_2(15:328),'g')
plot(1:314,x3_2(1:314),'b')
lgd = legend('Camera 1', 'Camera 2', 'Camera 3','Location','southeast');
lgd.FontSize =7;
hold off
title('Test 2: Noisy Case')
ylabel('Displacement')
xlabel('Frame')


figure(13)
plot(1:206,y1_3(20:225),'r')
hold on
plot(1:206  ,y2_3(5:210),'g')
plot(1:206,x3_3(10:215),'b')
hold off
title('Test 3: Horizontal Displacement')
ylabel('Displacement')
xlabel('Frame')
lgd = legend('Camera 1', 'Camera 2', 'Camera 3','Location','southeast');
lgd.FontSize =7;

figure(14)
plot(1:392,y1_4,'r')
hold on
plot(1:392  ,y2_4(1:392),'g')
plot(1:392,x3_4(1:392),'b')
hold off
title('Test 4: Horizontal Displacement and Rotation')
ylabel('Displacement')
xlabel('Frame')
lgd = legend('Camera 1', 'Camera 2', 'Camera 3','Location','southeast');
lgd.FontSize =7;

