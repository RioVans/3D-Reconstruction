function [ p ] = euc2proj( x)
[ss,s]=size(x);
p=x;
for i=1:s
    p(:,i)=x(:,i)/x(ss,i);
end

