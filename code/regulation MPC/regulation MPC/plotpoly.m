function p = plotpoly(p, varargin)
% p = plotpoly(p, ...)
% p = plotpoly({A, b}, ...)
% p = plotpoly(struct('A', A, 'b', b), ...)
%
% Plots the polyhedron with vertices given in the 2 by N matrix p or given by
% the extreme points of A*x <= b. Any additional arguments are passed to plot.
%
% Returns the extreme points matrix p. If you only want this matrix, pass
% false as the only additional argument, and the plot will not be made.
if nargin() < 1
    error('Argument p is required!');
end

% First, find the vertices if {A, b} given.
if iscell(p)
    p = halfspace2vertex(p{1}, p{2})';
elseif isstruct(p)
    p = halfspace2vertex(p.A, p.b)';
end

% Next, sort the vertices.
ptilde = bsxfun(@rdivide, bsxfun(@plus, p, -mean(p, 2)), std(p, 0, 2));
x = ptilde(1,:);
y = ptilde(2,:);
[th, r] = cart2pol(x, y);
thneg = (th < 0);
th(thneg) = th(thneg) + 2*pi(); % Makes theta in [0, 2*pi] instaed of [-pi, pi].
[~, s] = sortrows([th', r']); % Sort on theta then r.
p = p(:,s);

% Duplicate first data point to give closed cycle.
p = p(:,[1:end,1]);

% Now plot.
if isempty(varargin) || ~isequal(varargin{1}, false())
    plot(p(1,:), p(2,:), varargin{:}, 'LineWidth', 2);
end
end%function