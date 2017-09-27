function [x,l] = par2top(p, nx, ny)
    n_cp = length(p)/2; % number of control points
    cp = reshape(p,n_cp,2);
    
    % draw spline
    sp = spline(cp(:,1),cp(:,2));
    
    % convert spline to voxel
    x = zeros(ny, nx);
    X = 1:nx;
    Y = 1:ny;
    mesh = meshgrid(X,Y);
    