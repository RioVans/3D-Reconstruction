function [ Io ] = correct( I )
[m,n,o]=size(I);
if o==3
    Ir=rgb2gray(I);
end
Ir=imbinarize(Ir);
i=0;
U=ones(m,1);
D=0;
while((D~=0) & (i<n))
    i=i+1;
    C=Ir(:,i);
    D=C&U;
    
end
j=0;
U=ones(1,n);
D=0;
while((D~=0) & (j<m))
    j=j+1;
    C=Ir(j,:);
    D=C&U;
    
end
Io=I(1:j,1:i,:);
end

