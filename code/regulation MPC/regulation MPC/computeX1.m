function X1 = computeX1(Z, A, B, Xf)
% X1 = computeX1(Z, [A], [B], [Xf])
%
% Computes the feasible set X_1 for the system x^+ = Ax + Bu subject to
% constraints Gx + Hu <= psi and x^+ \in Xf.
%
% Z must be a struct with fields G, H, and psi.
%
% A and B are only necessary if the terminal constraint is given.
%
% Xf can be either a struct with fields A and b to define a polytope, or a
% single vector to define a point constraint. If not provided, it is assumed
% that Xf is the entire space.
%
% X1 is returned as a struct with fields A and b defining a set of inequality
% constraints.

% Check arguments.
narginchk(1, 4);
if ~isstruct(Z) || ~all(isfield(Z, {'G', 'H', 'psi'}))
    error('Invalid input for Z!');
end
Nx = size(Z.G, 2);
Nu = size(Z.H, 2);
if nargin() >= 3
    sys = struct('A', A, 'B', B); % Save these.
end

% Preallocate constraint matrices.
A = [Z.G, Z.H];
b = Z.psi;
Aeq = zeros(0, Nx + Nu);
beq = zeros(0, 1);

% Do some more stuff if there is a terminal constraint.
if nargin() >= 4 && ~isempty(Xf)
    if ~isequal([size(sys.A), size(sys.B)], [Nx, Nx, Nx, Nu])
        error('Incorrect size for A or B!');
    end
    if isstruct(Xf)
        % Polyhedron.
        if ~all(isfield(Xf, {'A', 'b'}))
            error('Struct Xf must have fields A and b!');
        end
        A = [A; Xf.A*[sys.A, sys.B]];
        b = [b; Xf.b];
    elseif isvector(Xf) && length(Xf) == Nx
        % Terminal equality constraint.
        Aeq = [sys.A, sys.B];
        beq = Xf;
    else
        error('Invalid input for Xf!');
    end
end

% Now do the projection step.
for i = 1:Nu
    [A, b, Aeq, beq] = fouriermotzkin(A, b, Aeq, beq);
end
X1 = struct('A', A, 'b', b);

end%function