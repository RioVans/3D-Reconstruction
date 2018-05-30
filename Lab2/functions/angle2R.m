function [ R ] = angle2R( angle,varargin)
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
if nargin == 2
if strcmp(varargin{1},'degrees')==1
    angle= angle*pi()/180;
else
    error('Wrong varargin. Use <<degrees>> or empty varagin (radians by default)')
end
end
R=[cos(angle),-sin(angle);sin(angle),cos(angle)];
end

