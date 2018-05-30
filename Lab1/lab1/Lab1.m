%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before



%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification  image 0000_s
clear all
% choose the image points
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

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'g');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');

%I=I(:,1:280,:);
v1=euc2proj(cross(l1,l2));
v2=euc2proj(cross(l3,l4));
v=cross(v1,v2);
v=euc2proj(v);
Hap=V2H([v(1),v(2)]);
I=I(:,1:274,:);
[Io,pedo]=apply_H(I,Hap);
Ql=Hap;
Ql(1:2,3)=pedo';
Lq1=euc2proj(inv(Ql)'*l1);%top horizontal
Lq2=euc2proj(inv(Ql)'*l2);%bottom horizontal
Lq3=euc2proj(inv(Ql)'*l3);%left vertical
Lq4=euc2proj(inv(Ql)'*l4);%right vertical

figure;
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
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'g');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'r');

v1=euc2proj(cross(l1,l2));
v2=euc2proj(cross(l3,l4));
v=cross(v1,v2);
v=euc2proj(v);
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

figure;
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

figure;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



