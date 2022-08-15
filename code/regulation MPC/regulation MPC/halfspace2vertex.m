function [V, nr] = halfspace2vertex(A, b, x0)
% [V, nr] = halfspace2vertex(A, b, [x0])
%
% Finds extreme points of polyhedron A*x <= b. Note that the polyhedron must
% have a point strictly on its interior.
%
% If provided, x0 must be a point on the interior of the polyhedron. If it is
% not given, one is found by solving a linear program.
%
% V is returned as an N by 2 matrix with each row giving an extreme point.
%
% Second output nr is a list of the non-redundant constraints of the polytope.

% Check inputs.
narginchk(2, 3);
Nc = size(A, 1);
Nx = size(A, 2);
if ~isvector(b) || length(b) ~= Nc
    error('b is the incorrect size!');
end
b = b(:); % Make sure b is a column vector.

% Sort out interior point.
if nargin() < 3
    if all(b > 0)
        % The origin is on the interior. Can rescale rows so that b = 1.
        x0 = zeros(Nx, 1);
        A = bsxfun(@rdivide, A, b);
        b = ones(size(b));
    else
        x0 = findinteriorpoint(A, b);
    end
elseif ~isvector(x0) || length(x0) ~= Nx
    error('Invalid size for x0!');
end
x0 = x0(:); % Make sure x0 is a column vector.

% Get non-redundant constraints from A and b.
[nr, ~, ~, k] = removeredundantcon(A, b, x0);

% Now loop through facets to find vertices.
V = zeros(size(k, 1), Nx);
keep = true(size(k, 1), 1);
for ix = 1:size(k, 1)
    F = A(k(ix,:),:);
    g = b(k(ix,:));
    [keep(ix), V(ix,:)] = fullranksolve(F, g);
end

V = V(keep,:);
[~, u] = unique(round(V*1e6), 'rows');
V = V(u,:);

% If in 2D, sort the vertices.
if Nx == 2
    V = polarsort(V);
end

end%function

function [fullrank, x] = fullranksolve(A, b);
    % Checks whether the system is full rank and if so, solves it. If it is not
    % full rank, a vector of NaNs are returned.
    Nx = size(A, 1);
    [U, S, V] = svd(A);
    s = diag(S);
    tol = Nx*eps(s(1)); % Rank tolerance.
    fullrank = all(s >= tol);
    if fullrank
        x = V*diag(1./s)*U'*b;
    else
        x = NaN(Nx, 1);
    end
end%function

function [p, s] = polarsort(p)
    % [p, s] = polarsort(p)
    %
    % Sorts the [n by 2] matrix p so that the points are counter-clockwise starting
    % at the theta = 0 axis. For ties in theta, sorts on radius.
    x = p(:,1);
    y = p(:,2);
    x = (x - mean(x))/std(x); % Normalize so that the origin is at the center.
    y = (y - mean(y))/std(y);
    [th, r] = cart2pol(x, y);
    [~, s] = sortrows([th, r]); % Sort on theta then r.
    p = p(s,:);
end%function
