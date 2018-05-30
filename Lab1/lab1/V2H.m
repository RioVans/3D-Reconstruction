function [H] = V2H( varargin )
%{  
    Juan Felipe Montesinos
    3D VISION
    Master in Computer Vision  2017-2018
    Computer Vision Center (Barcelona, Spain)
    ---------------------------------------------
    ---------------------------------------------
    Computes a 2D homography matrix. Affine, translation, similarity or
    projective.
    Inputs are typical parameters of each one. General expression is a 3x3
    matrix with 8 degrees of freedom.

    H=|A11 A12 t1|
      |A21 A22 t2|
      |v1  v2  1 |

    H=|A t|
      |v 1|
       
    A: 2x2 matrix
    v: 1x2 vector
    t: 2x1 vector
    
    Different notations shall provide wrong results
    
%}

%Parameters by default
A=[1, 0; 0, 1];
v=[0 0];
t=[0;0];

for i=1:nargin
    if ismatrix(varargin{i})==1 && isequal(size(varargin{i}),[2,2])==1 
        A=varargin{i};
    elseif ismatrix(varargin{i})==0 && isvector(varargin{i})==0
        error('Wrong input parameters. Your input is not a matrix or a vector')
    end
    if isvector(varargin{i})==1 && isequal(size(varargin{i}),[2,1])==1 
        t=varargin{i};
    end
    if isvector(varargin{i})==1 && isequal(size(varargin{i}),[1,2])==1 
        v=varargin{i};
    end
end
H=[A,t;v,1];
end

