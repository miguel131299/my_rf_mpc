function H = hatMap(a)
% convert from vector in R^3 to skew-symmetric matrix in so(3) 

H=[0 -a(3) a(2);
   a(3) 0 -a(1);
   -a(2) a(1) 0];

end