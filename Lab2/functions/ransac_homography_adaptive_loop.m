function [H, idx_inliers] = ransac_homography_adaptive_loop(x1, x2, th, max_it)

[Ncoords, Npoints] = size(x1);
%x1=xab_a;x2=xab_b;max_it=10000;


% ransac
it = 0;
best_inliers = [];
while it < max_it
    
    points = randomsample(Npoints, 4);
%      H = homography2d_copy(x1(:,points), x2(:,points));
    H = homography2d(x1, x2,points); % ToDo: you have to create this function
    inliers = compute_inliers(H, x1, x2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    p=0.99;
    max_it = log(1-p)/log(pNoOutliers);
    
    it = it + 1;
end

% compute H from all the inliers
% H = homography2d(x1(:,best_inliers), x2(:,best_inliers));
H = homography2d(x1, x2,best_inliers);
idx_inliers = best_inliers;


function idx_inliers = compute_inliers(H, x1, x2, th)
    % Check that H is invertible
    if abs(log(cond(H))) > 15
        idx_inliers = [];
        return
    end
    % compute the symmetric geometric error
    Hx1=H*x1;
    invHx2=H\x2;    
    x1= normalise(x1);
    x2= normalise(x2);     
    Hx1= normalise(Hx1);
    invHx2= normalise(invHx2); 
    
    d2=sum((x1-invHx2).^2)+sum((x2-Hx1).^2);
   
    idx_inliers = find(d2 < th.^2);


function xn = normalise(x)    
    xn = x ./ repmat(x(end,:), size(x,1), 1);

    
function item = randomsample(npts, n)
	a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
	  % Generate random value in the appropriate range 
	  r = ceil((npts-i+1).*rand);
	  item(i) = a(r);       % Select the rth element from the list
	  a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat