function [Y_initial] = gs_errfunction( P0, Xobs )

    H = reshape(P0(1:9), [3,3]);
    n_points = size(Xobs,1) / 2;
    x1 = Xobs(1:n_points); 
    x1 = reshape(x1, [2,size(x1,1)/2]);
    x2 = Xobs(n_points+1:end); 
    x2 = reshape(x2, [2,size(x2,1)/2]);     
    % reprojection error
    xhat = P0(9+1:end);
    xhat = reshape(xhat, [2,size(xhat,1)/2]); 
    xhat = [xhat ; ones(1,size(xhat,2))]; 
    xhatp = H*xhat;
    Y_initial = [disteuc(x1,euclid(xhat));disteuc(x2,euclid(xhatp))];

end

function dist = disteuc(v1,v2)
    dist = sqrt(sum((v1-v2).^2));
end