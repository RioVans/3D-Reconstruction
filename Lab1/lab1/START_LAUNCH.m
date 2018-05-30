%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification
display('3D Vision')
display('Master in Computer Vision  2017-2018')
display('------------------------------------------')
display('Juan Felipe Montesinos')
display('Ferran Carrasquer')
display('Yi Xiao')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
display('Similarity: Rotation 40 degrees, translation t(0 40)')
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
H=V2H(angle2R(40,'degrees'),[0;40]);
I2 = apply_H(I, H);
figure('Name',' Original Image','NumberTitle','off'); imshow(I); 
figure('Name','1.1. Similarity','NumberTitle','off'); imshow(uint8(I2));


%% 1.2. Affinities
close all
% ToDo: generate a matrix H which produces an affine transformation
display('Affinity: Rotation 20 degrees, scaling (1,2) at 40 degrees')
H=V2H(Affine(40,40,[1 1],'degrees'),[0;40]);
I2 = apply_H(I, H);
figure('Name','1.2. Affinity','NumberTitle','off'); imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
R=angle2R( 20,'degrees');
Rp=angle2R( 40,'degrees');
T=angle2R( -40,'degrees');
D=[1,0;0,2];
A=R*T*D*Rp;

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above
H=V2H(A,[0;40]);

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I3 = apply_H(I, H);
figure('Name','1.2. Composed Affinity','NumberTitle','off');imshow(uint8(I3));


%% 1.3 Projective transformations (homographies)
display('Projection: Rotation 20 degrees,vanishing (10-3,2 10-3)')
% ToDo: generate a matrix H which produces a projective transformation
H=V2H([0.001,0.002],Affine(20,40,[1 1],'degrees'),[0;60]);
I2 = apply_H(I, H);
figure('Name','1.3. Projection','NumberTitle','off'); imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification  image 0000_s
clear all
% choose the image points
display('-------------------------------------------------------')
display('-------------------------------------------------------')
display('Rectification')
display('Affine Rectification: Image 0000_s')
display('Vanishing line in red')
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 493; %top horizontal
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';

i = 186; %bottom horizontal
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';

i = 48; %left vertical
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';

i = 508; %right vertical
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

l1=cross(p1,p2);%top horizontal
l2=cross(p3,p4);%bottom horizontal
l3=cross(p5,p6);%left vertical
l4=cross(p7,p8);%right vertical




l1=cross(p1,p2);%top horizontal
l2=cross(p3,p4);%bottom horizontal
l3=cross(p5,p6);%left vertical
l4=cross(p7,p8);%right vertical

%I=I(:,1:280,:);
v1=euc2proj(cross(l1,l2));
v2=euc2proj(cross(l3,l4));
v=cross(v1,v2);
v=euc2proj(v);

% show the chosen lines in the imagefigure('Name','0000_s original image','off');imshow(I);
figure('Name','Original image 0000_s','NumberTitle','off');
imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'g');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
plot(t, -(v(1)*t + v(3)) / v(2), 'r');
Hap=V2H([v(1),v(2)]);

I=I(:,1:274,:);
[Io,pedo]=apply_H(I,Hap);
Ql=Hap;
Ql(1:2,3)=pedo';
Lq1=euc2proj(inv(Ql)'*l1);%top horizontal
Lq2=euc2proj(inv(Ql)'*l2);%bottom horizontal
Lq3=euc2proj(inv(Ql)'*l3);%left vertical
Lq4=euc2proj(inv(Ql)'*l4);%right vertical

figure('Name','Affine Rectified image','NumberTitle','off');
imshow(uint8(Io));
hold on;
t=1:0.1:1000;
plot(t, -(Lq1(1)*t + Lq1(3)) / Lq1(2), 'y');
plot(t, -(Lq2(1)*t + Lq2(3)) / Lq2(2), 'y');
plot(t, -(Lq3(1)*t + Lq3(3)) / Lq3(2), 'g');
plot(t, -(Lq4(1)*t + Lq4(3)) / Lq4(2), 'g');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 2. Affine Rectification  image 0001_s

clear all
close all
display('-------------------------------------------------------')
display('-------------------------------------------------------')
display('Rectification')
display('Affine & metric Rectification: Image 0001_s')
% choose the image points
I=imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

% indices of lines
i = 614; %top horizontal
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';

i = 159; %bottom horizontal
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';

i = 645; %left vertical
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';

i = 541; %right vertical
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

l1=cross(p1,p2);%top horizontal
l2=cross(p3,p4);%bottom horizontal
l3=cross(p5,p6);%left vertical
l4=cross(p7,p8);%right vertical
l5=cross(p5,p8);%diagonal 
l6=cross(p6,p7);%diagonal

% show the chosen lines in the image

v1=euc2proj(cross(l1,l2));
v2=euc2proj(cross(l3,l4));
v=cross(v1,v2);
v=euc2proj(v);
figure('Name','Original image 0001_s','NumberTitle','off');imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'g');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'r');

Hap=V2H([v(1),v(2)]);
I=I(:,1:573,:);
[Io,pedo]=apply_H(I,Hap);
Ql=Hap;
Ql(1:2,3)=pedo';
Lq1=euc2proj(inv(Ql)'*l1);%top horizontal
Lq2=euc2proj(inv(Ql)'*l2);%bottom horizontal
Lq3=euc2proj(inv(Ql)'*l3);%left vertical
Lq4=euc2proj(inv(Ql)'*l4);%right vertical
Lq5=euc2proj(inv(Ql)'*l5);%diagonal
Lq6=euc2proj(inv(Ql)'*l6);%diagonal

figure('Name','Affine rectified image 0001_s','NumberTitle','off')
imshow(uint8(Io));
hold on;
t=1:0.1:1000;
plot(t, -(Lq1(1)*t + Lq1(3)) / Lq1(2), 'y');
plot(t, -(Lq2(1)*t + Lq2(3)) / Lq2(2), 'y');
plot(t, -(Lq3(1)*t + Lq3(3)) / Lq3(2), 'g');
plot(t, -(Lq4(1)*t + Lq4(3)) / Lq4(2), 'g');
plot(t, -(Lq5(1)*t + Lq5(3)) / Lq5(2), 'r');
plot(t, -(Lq6(1)*t + Lq6(3)) / Lq6(2), 'r');


%% 3.1 Metric rectification after the affine rectification image 0001_s


L=[Lq1(1)*Lq4(1),Lq1(1)*Lq4(2)+Lq1(2)*Lq4(1),Lq1(2)*Lq4(2);...
 Lq5(1)*Lq6(1),Lq5(1)*Lq6(2)+Lq5(2)*Lq6(1),Lq5(2)*Lq6(2)];
 % Lq2(1)*Lq4(1),Lq2(1)*Lq4(2)+Lq2(2)*Lq4(1),Lq2(2)*Lq4(2);]
[~,~,V]=svd(L);
s=V(:,3);
S=[s(1),s(2);s(2),s(3)];
%S=nearestSPD(S);
K=chol(S,'lower');
K1=inv(K);
H=V2H(K1);
[I2,pedo]=apply_H(Io,H);
Ql=H;
Ql(1:2,3)=pedo';
Lq1=euc2proj(inv(Ql)'*Lq1);%top horizontal
Lq2=euc2proj(inv(Ql)'*Lq2);%bottom horizontal
Lq3=euc2proj(inv(Ql)'*Lq3);%left vertical
Lq4=euc2proj(inv(Ql)'*Lq4);%right vertical
Lq5=euc2proj(inv(Ql)'*Lq5);%diagonal
Lq6=euc2proj(inv(Ql)'*Lq6);%diagonal

figure('Name','Affine & metric rectified image 0001_s','NumberTitle','off')
imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(Lq1(1)*t + Lq1(3)) / Lq1(2), 'y');
plot(t, -(Lq2(1)*t + Lq2(3)) / Lq2(2), 'y');
plot(t, -(Lq3(1)*t + Lq3(3)) / Lq3(2), 'g');
plot(t, -(Lq4(1)*t + Lq4(3)) / Lq4(2), 'g');
plot(t, -(Lq5(1)*t + Lq5(3)) / Lq5(2), 'r');
plot(t, -(Lq6(1)*t + Lq6(3)) / Lq6(2), 'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



