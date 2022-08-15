function [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq, ielim)
% [A, b, Aeq, beq] = fouriermotzkin(A, b, [Aeq], [beq], [ielim])
%
% Perform one step of Fourier-Motzkin elimination for the set defined by
%
%   {x \in R^n : A*x <= b, Aeq*x = beq}
%
% Note that the inequality constrants can be degenerate, although degeneracy
% will create significantly more constraints compared to explicit equality
% constraints.
%
% Optional argument ielim decides which column to eliminate. Default is the last
% column.

% Check arguments and get sizes.
narginchk(2, 5);
Nlt = size(A, 1);
if ~isvector(b) || length(b) ~= Nlt
    error('Invalid size for b!');
end
Nx = size(A, 2);
if nargin() < 3 || isempty(Aeq)
    Aeq = zeros(0, Nx);
    beq = zeros(0, 1);
    Neq = 0;
elseif nargin() < 4
    error('beq is required if Aeq is given!');
else
    if size(Aeq, 2) ~= Nx
        error('Aeq has the wrong number of columns!');
    end
    Neq = size(Aeq, 1);
    if ~isvector(beq) || length(beq) ~= Neq
        error('Invalid size for beq!');
    end
end
if nargin() < 5
    ielim = Nx;
elseif ~isscalar(ielim) || ielim <= 0 || ielim > Nx || round(ielim) ~= ielim
    error('ielim must be a scalar positive integer less than Nx!');
end

% Now decide what to do. First, look for an equality constraint with the
% variable of interest.
[pivval, pivrow] = max(abs(Aeq(:,ielim)));
if isempty(pivrow) || pivval == 0
    % No suitable equality constraint found. Need to change inequality
    % constraints.
    [A, b] = fmelim(A, b, ielim);
    Aeq(:,ielim) = [];
else
    % An equality constraint is available. Handle that.
    a = Aeq(pivrow,:);
    c = a(ielim);

    % Make the pivots.
    Aeq = Aeq - Aeq(:,ielim)*a/c;
    beq = beq - Aeq(:,ielim)*beq(pivrow)/c;

    A = A - A(:,ielim)*a/c;
    b = b - A(:,ielim)*beq(pivrow)/c;

    % Get rid of the appropriate rows and columns.
    A(:,ielim) = [];
    Aeq(:,ielim) = [];
    Aeq(pivrow,:) = [];
    beq(pivrow) = [];
end

% Zap any rows that are all zeros.
[A, b] = removezerorows(A, b);
[Aeq, beq] = removezerorows(Aeq, beq);

end%function


% *****************************************************************************
% Helper Functions
% *****************************************************************************

function [Ae, be] = fmelim(A, b, ielim)
    % Performs one step of Fourier-Motzkin elimination for inequality
    % constraints.
    c = A(:,ielim);
    I0 = find(c == 0);
    Ip = find(c > 0);
    Im = find(c < 0);

    Nx = size(A, 2);
    Ne = length(I0) + length(Ip)*length(Im);

    E = [A(:,1:ielim - 1), A(:,ielim + 1:end), b];
    Ee = [E(I0,:); kron(c(Ip), E(Im,:)) - kron(E(Ip,:), c(Im))];

    Ae = Ee(:,1:end-1);
    be = Ee(:,end);
end%function

function [A, b] = removezerorows(A, b)
    % Removes rows of A that are all zeros.
    keeprows = any(A, 2);
    A = A(keeprows,:);
    b = b(keeprows);
end%function