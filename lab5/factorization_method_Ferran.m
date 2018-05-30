function [Pproj,Xproj]=factorization_method(x1,x2)

%{  
  Yi Xiao
  3D VISION
  Master in Computer Vision  2017-2018
  Computer Vision Center (Barcelona, Spain)
  ---------------------------------------------
  ---------------------------------------------
  The entries are x1,x2, where:
  - x1 and x2: homogeneous coordinates of 2D points in two images.

  The output are Pproj and Xproj, where:
  - Pproj: 3*Ncam x 4 matrix containing the camera matrices
  - Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
%}

%normalization
[x1_norm, T1] = normalise2dpts(x1);
[x2_norm, T2] = normalise2dpts(x2);
%camer 1
F1 = fundamental_matrix(x1, x1);
[U1 S1 V1]=svd(F1);
e1 = V1(:,end)/V1(3,3);% here
% camera 2
F2 = fundamental_matrix(x2, x1);
[U1 S1 V1]=svd(F2);
e2 = V1(:,end)/V1(3,3);% here
d_old=0;

%Initialize all alpha1 = 1
A = ones(2,length(x1));
for j=1:length(x1)
    A(1,j) = x1(:,j)'* F1 *(cross(e1,x1(:,j)))/norm(cross(e1,x1(:,j))).^2*A(1,j);
end

for j=1:length(x2)
    A(2,j) = x1(:,j)'* F2 *(cross(e2,x2(:,j)))/norm(cross(e2,x2(:,j))).^2*A(1,j);
end

while 1
    alpha1_j_3rows = [A(1,:);A(1,:);A(1,:)];
    alpha2_j_3rows = [A(2,:);A(2,:);A(2,:)];
    W = [alpha1_j_3rows.*x1_norm;alpha2_j_3rows.*x2_norm];

    % Compute the SVD of the matrix W
    [u d v]=svd(W);
    v_p=v';

    Pproj = u*d(:,1:4);
    Xproj = v_p(1:4,:);

    for i=1:2
        x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
    end
    x_d{1} = euclid(x1_norm);
    x_d{2} = euclid(x2_norm);
    for i = 1:length(x1_norm)
        d1(1,i)=sqrt((x_d{1}(1,i)-x_proj{1}(1,i)).^2+(x_d{1}(2,i)-x_proj{1}(2,i)).^2);
        d2(1,i)=sqrt((x_d{2}(1,i)-x_proj{2}(1,i)).^2+(x_d{2}(2,i)-x_proj{2}(2,i)).^2);
    end
    d = sum(sum(d1.^2+d2.^2));
    if (abs(d - d_old)/d)< 0.1
        break;
    else
        d_old = d;
        temp1 = Pproj(1:3,:)*Xproj;
        A(1,:) = temp1(3,:);
        temp2 = Pproj(4:6,:)*Xproj;
        A(2,:) = temp2(3,:);    
    end
end

Pproj(1:3,:)=inv(T1)*Pproj(1:3,:);
Pproj(4:6,:)=inv(T2)*Pproj(4:6,:);
end

