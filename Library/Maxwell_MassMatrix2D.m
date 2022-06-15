function M_global = Maxwell_MassMatrix2D

% function M_global = Maxwell_MassMatrix2D
% Computes the global Mass matrix of 2D Maxwell problem
% including Jacobian !! (Diference with Hesthaven, Warburton)
% order of variables: [Hx^1; Hy^1; Ez^1; ... Hx^K; Hy^K; Ez^K] = u
% where Hx^k,Hy^k Ez^k are Np x 1 vectors
% Then: M du/dt = S u


Globals2D;


% ----------------------------------------------------------------------
% DEFINE "MASS MATRIX" M_global ------------------------------------------
% ----------------------------------------------------------------------

M_global=spalloc(3*Np*K, 3*Np*K, K*9*Np^2);
%M_global=zeros(3*Np*K);

Zero=zeros(Np);
M_global_k = [MassMatrix Zero Zero; Zero MassMatrix Zero; Zero Zero MassMatrix];


for k=1:K
    
    M_global( (k-1)*3*Np+1: k*3*Np,  (k-1)*3*Np+1: k*3*Np) = J(1,k)*M_global_k;
 
end

%M_global=sparse(M_global);