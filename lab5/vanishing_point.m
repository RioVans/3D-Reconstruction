function [v] = vanishing_point(xo1, xf1, xo2, xf2)
%{  
  Yi Xiao
  3D VISION
  Master in Computer Vision  2017-2018
  Computer Vision Center (Barcelona, Spain)
  ---------------------------------------------
  ---------------------------------------------
  Computes the vanishing point formed by the line that joins points xo1 and xf1 and the line 
  that joins points x02 and xf2.

  The entries are xo1, xf1, xo2, xf2
  The output is the vanishing point
%}

l1 = cross(xo1, xf1);
l1 = l1/l1(end);

l2 = cross(xo2, xf2);
l2 = l2/l2(end);

v = cross(l1, l2);

end
