I=imread('penta.png');
%Set of points:
%Red/Sky Blue/Dark Blue/Orange/Purple/Dark Orange/Green
%Original points
x2=[251 407 497 403 97 128 2;4 72 182 473 473 267 182;1 1 1 1 1 1 1];
%transformed points H1
x1=[374 611 743 604 144 190 5; 3 54 137 354 354 201 138; 1 1 1 1 1 1 1];
%{
H1_solution =

    1.5000   -0.0000    0.0000
   -0.0000    0.7500    0.0000
   -0.0000   -0.0000    1.0000

%}
imtool(I)

[H1] = homography2d(x1, x2);
imtool(apply_H(I,H1))