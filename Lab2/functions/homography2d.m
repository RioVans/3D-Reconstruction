function [H]=homography2d(x1,x2,points)
i=1;
p=1;
u=length(points);
A=[]; 
r=1;
suma=0;
sa1=x1(:,points);
sa2=x2(:,points);
s1x_avg=sum(sa1(1,:))/u;
s1y_avg=sum(sa1(2,:))/u;
s2x_avg=sum(sa2(1,:))/u;
s2y_avg=sum(sa2(2,:))/u;
for p=1:u
    d=x1(:,points(p));
    pre=((d(1)-s1x_avg)^2+(d(2)-s1y_avg)^2)^0.5;
    suma=suma+pre;
    p=p+1;
end
s1=2^0.5/(1/u*suma);
T1=[s1,0,-s1*s1x_avg;0,s1,-s1*s1y_avg;0,0,1];
x1new=T1*x1;
p=1;
suma=0;
for p=1:u
    d=x2(:,points(p));
    pre=((d(1)-s2x_avg)^2+(d(2)-s2y_avg)^2)^0.5;
    suma=suma+pre;
    p=p+1;
end
s2=2^0.5/(1/u*suma);
T2=[s2,0,-s2*s2x_avg;0,s2,-s2*s2y_avg;0,0,1];
x2new=T2*x2;

for i=1:u
    a=x1new(:,points(i));
    b=x2new(:,points(i));
    A(r,:)=[0,0,0,-a(1),-a(2),-1,b(2)*a(1),b(2)*a(2),b(2)];
    A(r+1,:)=[a(1),a(2),1,0,0,0,-a(1)*b(1),-b(1)*a(2),-b(1)];
    r=r+2;
    i=i+1;
end
    [U,D,V] = svd(A);
    H = reshape(V(:,9),3,3)';
    H=-H;
    H = T2\H*T1;
% H = H/H(3,3);
end
