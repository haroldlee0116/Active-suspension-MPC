function predmod=predmodgen_outputMPC(LTI,dim)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

%Prediction matrix from initial state
T=zeros(dim.nx*(dim.N+1),dim.nx);
for k=0:dim.N
    T(k*dim.nx+1:(k+1)*dim.nx,:)=LTI.A^k;
end

%Prediction matrix from input
S=zeros(dim.nx*(dim.N+1),dim.nu*(dim.N));
for k=1:dim.N
    for i=0:k-1
        S(k*dim.nx+1:(k+1)*dim.nx,i*dim.nu+1:(i+1)*dim.nu)=LTI.A^(k-1-i)*LTI.B;
    end
end

predmod.T=T;
% LTI.D=[LTI.D];
% S=S+blkdiag(LTI.D,LTI.D,LTI.D,LTI.D,LTI.D);
predmod.S=S;

% P=[LTI.C];
% % P=eye(dim.nx);
% 
% 
% for k=1:dim.N-1 % k=1:dim.N
%     P=[P;LTI.C*(LTI.A^k)];
% %   P=[P;LTI.A^k];
% end
% 
% % Prediction matrix from input
% S = zeros(dim.ny,dim.nu*dim.N);
% S_prec = S;
% 
% for k = 1:dim.N-1 % k= 1:dim_N
%        S = [S;circshift(S_prec,dim.nu,2)];
% %      S(end-dim.nx+1:end,1:dim.nu) = LTI.A^(k-1)*LTI.B;
%        S(end-dim.ny+1:end,1:dim.nu) = LTI.C*(LTI.A^(k-1))*LTI.B;
%        S_prec = S(end-dim.ny+1:end,:);
% end
% LTI.D=[LTI.D];
% S=S+blkdiag(LTI.D,LTI.D,LTI.D,LTI.D,LTI.D);
% predmod.T=P;
% predmod.S=S;
% end


