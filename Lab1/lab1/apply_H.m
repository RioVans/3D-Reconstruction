function [output,vargout]=apply_H(I,H)
%{  
    Yi Xiao
    3D VISION
    Master in Computer Vision  2017-2018
    Computer Vision Center (Barcelona, Spain)
    ---------------------------------------------
    ---------------------------------------------
    Computes the projection H over image I.
   
    input: H - The homography
           I - The original image

    output: output - the output image after applying the H

%}

[m,n,o]=size(I);
c1=[1;1;1];
c2=[1;m;1];
c3=[n;1;1];
c4=[n;m;1];

c1_H= H*c1;
c2_H= H*c2;
c3_H= H*c3;
c4_H= H*c4;

c1_H = round(c1_H./c1_H(3,:));
c2_H = round(c2_H./c2_H(3,:));
c3_H = round(c3_H./c3_H(3,:));
c4_H = round(c4_H./c4_H(3,:));

x_min_output_image=min([c1_H(1,1),c2_H(1,1),c3_H(1,1),c4_H(1,1)]);
x_max_output_image=max([c1_H(1,1),c2_H(1,1),c3_H(1,1),c4_H(1,1)]);
y_min_output_image=min([c1_H(2,1),c2_H(2,1),c3_H(2,1),c4_H(2,1)]);
y_max_output_image=max([c1_H(2,1),c2_H(2,1),c3_H(2,1),c4_H(2,1)]);

[X,Y] = meshgrid(x_min_output_image:x_max_output_image, y_min_output_image:y_max_output_image);
cols_output_image = x_max_output_image - x_min_output_image + 1;
rows_output_image = y_max_output_image - y_min_output_image + 1;
Z = ones(rows_output_image,cols_output_image);

projected_output_coordinates = inv(H)*[X(:) Y(:) Z(:)]';

temp=projected_output_coordinates(3,:);
projected_output_coordinates_nor = projected_output_coordinates(1,:)./temp;
HX = reshape(projected_output_coordinates_nor, rows_output_image, cols_output_image);
projected_output_coordinates_nor = projected_output_coordinates(2,:)./temp;
HY = reshape(projected_output_coordinates_nor, rows_output_image, cols_output_image);

output_image(:,:,1) = interp2(double(I(:,:,1)), HX, HY, 'linear', 0);
output_image(:,:,2) = interp2(double(I(:,:,2)), HX, HY, 'linear', 0);
output_image(:,:,3) = interp2(double(I(:,:,3)), HX, HY, 'linear', 0);

vargout = [abs(x_min_output_image),abs(y_min_output_image)];
output=output_image;
end