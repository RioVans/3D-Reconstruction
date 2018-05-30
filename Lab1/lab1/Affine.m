function [ A] = Affine( angle1,angle2,D,varargin )
%{  
    Juan Felipe Montesinos
    3D VISION
    Master in Computer Vision  2017-2018
    Computer Vision Center (Barcelona, Spain)
    ---------------------------------------------
    ---------------------------------------------
    For a given angle it computes the 2D rotation matrix in Euclidean
    space considering counterclockwise rotation.

    Inputs:
    -angle: angle in radians or degrees*
    -varargin: "degrees" in case of input angle as degrees.
    
    Output:
    -R:rotation matrix.
%}
D=[D(1),0;0,D(2)];
if nargin==4
    R=angle2R( angle1,varargin);
    Rp=angle2R( angle2,varargin);
    Rn=angle2R( -angle2,varargin);
else
    R=angle2R( angle1);
    Rp=angle2R( angle2);
    Rn=angle2R( -angle2);
end
A=R*Rn*D*Rp;
end

