function [H,h]=costgen(predmod,weight,dim)

Qbar=blkdiag(kron(eye(dim.N),weight.Q),weight.P);
Rbar=kron(eye(dim.N),weight.R);
H=predmod.S'*Qbar*predmod.S+Rbar;   
h=predmod.S'*Qbar*predmod.T;
 
end