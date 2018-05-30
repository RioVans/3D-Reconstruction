%% Read template and images.
clear all
close all
addpath('sift');
addpath('functions');

UPF     = imread('Data/logos/logoUPF.png');
%ED  = imread('Data/logos/UPFbuilding.jpg');
ED  = imread('Data/logos/UPFstand.jpg');
CVC  = imread('Data/logos/logo_master.png');
CVC=imresize(CVC,[122,123]);
upf = sum(double(UPF), 3) / 3 / 255;
ed = sum(double(ED), 3) / 3 / 255;
cvc = sum(double(CVC), 3) / 3 / 255;

%% Compute SIFT keypoints
[points_upf, desc_upf] = sift(upf, 'Threshold', 0.005);
[points_ed, desc_ed] = sift(ed, 'Threshold', 0.01);



figure;
imshow(UPF);%image(imargb)
hold on;
plot(points_upf(1,:), points_upf(2,:),'+g');
figure;
imshow(ED);%image(imbrgb);
hold on;
plot(points_ed(1,:), points_ed(2,:),'+y');


%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_upf, desc_ed);
figure;
plotmatches(upf, ed, points_upf(1:2,:), points_ed(1:2,:), matches_ab, 'Stacking', 'v');
%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_upf(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_ed(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
%[H]=homography2d(xab_a, xab_b,1:1:length(matches_ab));
[H, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000);


%vgg_gui_H(UPF, ED, H);
[M,N]=size(ed);
cvc2=apply_H_v2(double(255*cvc),H,[1,N,1,M]);
A=(cvc2>uint8(255*ed));
CVC2=apply_H_v2(CVC,H,[1,N,1,M]);
for i=1:3
   
    ED2(:,:,i)=CVC2(:,:,i)+uint8(~A).*ED(:,:,i);
end
figure;
imshow([ED2,ED])

