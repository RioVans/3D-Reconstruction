function [Pproj,Xproj]=factorization_method(x,Ncam)

%{  
  Yi Xiao
  3D VISION
  Master in Computer Vision  2017-2018
  Computer Vision Center (Barcelona, Spain)
  ---------------------------------------------
  ---------------------------------------------
  The entries are x1,x2, where:
  - x: the cell that contains the homogeneous coordinates of 2D points in two images x1 and x2.
  - Ncam: the number of the cameras

  The output are Pproj and Xproj, where:
  - Pproj: 3*Ncam x 4 matrix containing the camera matrices
  - Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
%}

%normalization
for i = 1:Ncam
    [xi_norm, T] = normalise2dpts(x{i});
    x_norm{i} = xi_norm;
    T_m{i} = T;
end

d_old=0;
%Initialize the first row of alpha = 1
lambda = ones(Ncam,length(x_norm{1}));

for i =1:Ncam
    F = fundamental_matrix(x_norm{i},x_norm{1});
    F_m{i} = F;
    [U S V]= svd( F_m{i});
    e= V(:,3)/V(3,3);
    e_m{i} = e;
end

for i = 1:Ncam
    for j = 1: length(x{1})
        lambda(i,j) = x_norm{1}(:,j)'*F_m{i}*cross(e_m{i},x_norm{i}(:,j))/norm(cross(e_m{i},x_norm{i}(:,j))).^2 * lambda(1,j);
    end
end

while 1
    alpha1_j_3rows = [lambda(1,:);lambda(1,:);lambda(1,:)];
    alpha2_j_3rows = [lambda(2,:);lambda(2,:);lambda(2,:)];
    W = [alpha1_j_3rows.* x_norm{1};alpha2_j_3rows.* x_norm{2}];

    % Compute the SVD of the matrix W
    [u d v]=svd(W);

    Pproj = u*d(:,1:4);
    Xproj = v(:,1:4)';

    for i=1:2
        x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
    end
    x_d{1} = euclid( x_norm{1});
    x_d{2} = euclid( x_norm{2});
    for i = 1:length( x_norm{1})
        d1(1,i)=sqrt((x_d{1}(1,i)-x_proj{1}(1,i)).^2+(x_d{1}(2,i)-x_proj{1}(2,i)).^2);
        d2(1,i)=sqrt((x_d{2}(1,i)-x_proj{2}(1,i)).^2+(x_d{2}(2,i)-x_proj{2}(2,i)).^2);
    end
    d = sum(sum(d1.^2+d2.^2));
    if (abs(d - d_old)/d)< 0.1
        break;
    else
        d_old = d;
        temp1 = Pproj(1:3,:)*Xproj;
        lambda(1,:) = temp1(3,:);
        temp2 = Pproj(4:6,:)*Xproj;
        lambda(2,:) = temp2(3,:);    
    end
end

Pproj(1:3,:)=inv(T_m{1})*Pproj(1:3,:);
Pproj(4:6,:)=inv(T_m{2})*Pproj(4:6,:);
end