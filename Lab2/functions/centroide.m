function [centroid]=centroide(M)
[m,n]=size(M);
X=sum(M(1,:))/n;
Y=sum(M(2,:))/n;
centroid=[X;Y];