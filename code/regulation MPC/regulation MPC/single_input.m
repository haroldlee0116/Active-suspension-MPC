function [P, gamma] = single_input(G, H, Psi)

% Define the sets by basing on the i-th column of H
s = size(Psi,1);

I_0 = [];
I_p = [];
I_m = [];

for i = 1:s
    if H(i) == 0
        I_0 = [I_0 i];
    end
    if H(i) > 0
        I_p = [I_p i];
    end
    if H(i) < 0
        I_m = [I_m i];
    end
end

s_0 = numel(I_0);%元素个数
s_p = numel(I_p);
s_m = numel(I_m);

% Set the row of matrix [P gamma]

% Define C
C = [G Psi];

aux = [];

% Define row by row [P gamma]
for i = I_0
    aux = [aux; C(i,:)];
end

for i = I_p
   for j = I_m
      aux = [aux; H(i,:)*C(j,:) - H(j,:)*C(i,:)];
   end
end

% Return the desired matrix/vector
P = aux(:,1:end-1);
gamma = -aux(:,end); % gamma = -aux(:,end);


