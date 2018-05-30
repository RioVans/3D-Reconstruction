function [F_es] = fundamental_matrix(x1_test, x2_test)
%function that estimates F using the normalised 8 point algorithm
    
    [x1_norm, t1] = normalise2dpts(x1_test); 
    [x2_norm, t2] = normalise2dpts(x2_test); 
    
    u1 = x1_norm(1,:)';
    v1 = x1_norm(2,:)';  
    u2 = x2_norm(1,:)';
    v2 = x2_norm(2,:)';

    A = [u1.*u2 v1.*u2 u2 u1.*v2 v1.*v2 v2 u1 v1 ones(size(u1,1),1)];

    [u, d, v] = svd(A);

    F_es = v(:,9);
    F_es = [F_es(1) F_es(2) F_es(3);...
            F_es(4) F_es(5) F_es(6);...
            F_es(7) F_es(8) F_es(9)];

    [u, d, v] = svd(F_es);
    d(3,3) = 0;
    F_es = u*d*v';

    F_es = t2'*F_es*t1;
    
end

