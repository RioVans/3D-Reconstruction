function [ Io ,vargout] = xHx(H,I,varargin)
%{  
    Juan Felipe Montesinos
    3D VISION
    Master in Computer Vision  2017-2018
    Computer Vision Center (Barcelona, Spain)
    ---------------------------------------------
    ---------------------------------------------
    Computes the projection H over image I. To avoid quantization artifacts
    (gaps produced in the projection due to discontinuity of images
    quantized in pixels) it is computed a similarity to modify homography
    in order to obtain the minimum size without artifacts and later
    upsample the image.
    
    varargin is rescaling threshold (0.5 is recommended value). No varargin
    means no fixation of artifacts
%}
if  ismatrix(H)==0 || isequal(size(H),[3,3])==0 || H(3,3) ~=1
    error('H is not a property 2D homography')
end
[m,n,o]=size(I);
A=H(1:2,1:2);
s=det(A);
if nargin==3
if s>varargin{1}
    A=sqrt(varargin{1})*A/sqrt(s);
    H(1:2,1:2)=A;
end
end
ext=floor(euc2proj(H*[m;n;1]));
a=floor(euc2proj(H*[m;0;1]));
b=floor(euc2proj(H*[0;n;1]));
c=floor(euc2proj(H*[1;1;1]));
L=zeros(abs(ext(1)),abs(ext(2)));
comp=[1,1];

if ext(1)<1 || a(1)<1 || b(1)<1 || c(1)<1
    comp=comp-[min([ext(1),a(1),b(1),c(1)]),0];
end
if ext(2)<2 || a(2)<1 || b(2)<1 || c(2)<1
    comp=comp-[0,min([ext(2),a(2),b(2),c(2)])];
end
for k=1:o
for i=1:m
    for j=1:n
        v=euc2proj(H*[i;j;1]);
        v=round(v);
        L(v(1)+comp(1),v(2)+comp(2))=I(i,j,k);
    end

end
    Io(:,:,k)=L;
end
if nargin==3
if s>varargin{1}
    Io=imresize(Io,sqrt(s/varargin{1}));
    vargout=comp*sqrt(s/varargin{1});    
end
end
vargout=comp;
Io=uint8(Io);
