function X = triangulate(x1, x2, P1, P2, imsize)
%{  
  Yi Xiao
  3D VISION
  Master in Computer Vision  2017-2018
  Computer Vision Center (Barcelona, Spain)
  ---------------------------------------------
  ---------------------------------------------
  The entries are (x1, x2, P1, P2, imsize), where:
  - x1, and x2 are the Euclidean coordinates of two matching points in two different images.
  - P1 and P2 are the two camera matrices
  - imsize is a two-dimensional vector with the image size

  The output:
  - X: the triangulation
%}
x1(3,:)=1;
x2(3,:)=1;

%Scaling and translation so that both pixel coordinates are in the interval [-1, 1].

H=[2/imsize(1) 0 -1;0 2/imsize(2) -1;0 0 1];
x1=H*x1;
x2=H*x2;
P1=H*P1;
P2=H*P2;

% to get the A matrix
A(1,:) = x1(1,:)*P1(3,:)-P1(1,:);
A(2,:) = x1(2,:)*P1(3,:)-P1(2,:);
A(3,:) = x2(1,:)*P2(3,:)-P2(1,:);
A(4,:) = x2(2,:)*P2(3,:)-P2(2,:);

% to do the SVD
[U D V] = svd(A);
X = V(:,4);
X=X./X(end);
end
